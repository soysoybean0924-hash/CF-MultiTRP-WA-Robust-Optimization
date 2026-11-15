"""Prototype layered protection runner. Run only on the Linux server.

No claim of a strict robust counterpart or independent physical coverage.
"""
import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter
import numpy as np
import native_solver as solver
from frozen_input import load
from source_policy import classify_links
from service_policy import protected_service
from campaign import evaluate, covariance_change

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def atomic_json(path,obj):
    tmp=path.with_suffix(path.suffix+'.partial')
    if tmp.exists():raise FileExistsError(tmp)
    tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)),encoding='utf-8')
    tmp.replace(path)

def calibrate(root,config,out):
    scene=load(root)
    cfg=json.loads(Path(config).read_text(encoding='utf-8'))
    dev=set(cfg['split']['development_groups']);test=set(cfg['split']['test_groups'])
    if dev & test or dev|test!=set(range(1,len(scene.resources)+1)):
        raise ValueError('development/test group partition mismatch')
    records=[];coverage=[]
    for g,r in enumerate(scene.resources,1):
        labels,stats=classify_links(r,scene.prior_rx)
        tr=stats['trace'];total=tr.sum(axis=-1)
        # Source fractions are descriptive; full C remains the only physical covariance.
        share=tr/np.maximum(total[...,None],1e-30)
        measured=np.mean(np.abs(r.h_true-r.h_est)**2,axis=(2,3))
        predicted=np.trace(r.error_var,axis1=2,axis2=3).real/4
        normalized=measured/np.maximum(predicted,1e-30)
        for k,ue in enumerate(r.users):
            for l in range(r.h_est.shape[1]):
                records.append(dict(group=g,ue=int(ue),trp=l,split='development' if g in dev else 'test',
                    label=str(labels[k,l]),total_trace=float(total[k,l]),
                    source_share=share[k,l].tolist(),normalized_error=float(normalized[k,l])))
    train=np.array([x['normalized_error'] for x in records if x['split']=='development'])
    held=np.array([x['normalized_error'] for x in records if x['split']=='test'])
    radius=float(np.quantile(train,cfg['coverage_target'],method='higher'))
    result=dict(schema=1,status='model_based_group_holdout_prototype',
      input_manifest_sha256=sha(Path(root)/'manifest.json'),config_sha256=sha(config),
      coverage_target=cfg['coverage_target'],radius_squared=radius,
      uncertainty_set='per-link empirical mean-square ball: mean_tx |H-Hhat|^2 <= radius_squared * trace(C)/Nr',
      covariance_semantics='C is physical posterior; source traces are not summed into C again',
      correlation='source cross terms not stored in frozen input; full C used for total, no independent-source guarantee',
      calibration_scope='one scene, group holdout; not independent UE/channel/SRS seed coverage',
      development_count=len(train),test_count=len(held),
      development_coverage=float(np.mean(train<=radius)),test_coverage=float(np.mean(held<=radius)),
      source_names=['residual','pilot','noise','age'],records=records)
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    atomic_json(out/'layer1.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='records'},ensure_ascii=False,indent=2))

def protection(scene,layer,config,mode,level,ablation='none',trigger='proposed'):
    cfg=json.loads(Path(config).read_text(encoding='utf-8'))
    if not 1<=len(scene.resources)<=17:raise ValueError('resource group count outside 1..17')
    rec={(x['group'],x['ue'],x['trp']):x for x in layer['records']}
    shares=np.array([[rec[(g,int(ue),l)]['source_share'] for l in range(r.h_est.shape[1])]
                     for g,r in enumerate(scene.resources,1) for ue in r.users],dtype=float)
    weights=cfg['mapping'];src=np.array([weights['residual_weight'],weights['pilot_weight'],
                     weights['noise_weight'],weights['age_weight']])
    risk=float(np.mean(shares@src))
    mapped=float(np.clip(weights['floor']+risk,weights['floor'],weights['ceiling']))
    calibration_factor=float(np.sqrt(layer['radius_squared']))
    # Layer 2 uses one calibrated common level; only the per-variable allocation differs.
    # Layer 3 keeps variable factors fixed and compares fixed against mapped level.
    base=(level if mode in ('uniform','variable','fixed') else mapped*level/max(cfg['fixed_level'],1e-12))*calibration_factor
    factors=cfg['variable_factors'] if mode in ('variable','fixed','mapped') else dict(beam=1,service=1,power=1)
    strengths={k:float(np.clip(base*factors[k],0,0.3)) for k in ('beam','service','power')}
    if ablation in strengths:strengths[ablation]=0.0
    old=[r.serve.copy() for r in scene.resources]
    service_diagnostics=[]
    if trigger!='off' and strengths['service']>0:
        margin=cfg['service']['fixed_margin'] if trigger=='fixed' else cfg['service']['trigger_margin']
        selected,service_diagnostics=protected_service(scene.resources,strengths['service'],
          'uniform' if mode=='uniform' else 'source',margin,
          cfg['service']['load_penalty'],cfg['service_source_weights'])
        for r,m in zip(scene.resources,selected):r.serve=m
    changes=int(sum(np.count_nonzero(r.serve!=m) for r,m in zip(scene.resources,old)))
    # Extra P regularizes the statistical WMMSE design; it is not a strict robust counterpart.
    penalty=[strengths['beam']*r.error_var for r in scene.resources]
    physical_cap=float(scene.config.power_watt_per_trp)
    effective_cap=physical_cap*(1-cfg['power_reserve_fraction_max']*strengths['power']/0.3)
    scene.config.power_watt_per_trp=effective_cap
    return penalty,dict(mode=mode,level=level,strengths=strengths,source_risk=risk,
      calibration_factor=calibration_factor,
      mapped_level=mapped,physical_cap_w=physical_cap,effective_cap_w=effective_cap,
      changed_service_links=changes,trigger=trigger,ablation=ablation,
      service_rule='uniform_total_risk' if mode=='uniform' else 'source_weighted_risk',
      service_diagnostics=service_diagnostics,
      beam_penalty_trace=float(sum(np.trace(p,axis1=2,axis2=3).real.sum() for p in penalty)))

def run(root,config,layer_path,mode,out,level=None,ablation='none',trigger='proposed',max_outer=None,power_cap=None,max_groups=None):
    out=Path(out)
    if out.exists():raise FileExistsError(f'use a new result directory: {out}')
    layer=json.loads(Path(layer_path).read_text(encoding='utf-8'))
    if layer['input_manifest_sha256']!=sha(Path(root)/'manifest.json'):
        raise ValueError('layer1 frozen-input identity mismatch')
    if layer.get('config_sha256')!=sha(config):
        raise ValueError('layer1/config identity mismatch; recalibrate in a new directory')
    scene=load(root)
    if max_groups is not None:
        if not 1<=max_groups<=len(scene.resources):raise ValueError('invalid max_groups')
        scene.resources=scene.resources[:max_groups]
    if power_cap is not None:
        if not 0<power_cap<=scene.config.power_watt_per_trp:raise ValueError('power cap outside 0..physical cap')
        scene.config.power_watt_per_trp=float(power_cap)
    cfg=json.loads(Path(config).read_text(encoding='utf-8'))
    level=float(cfg['fixed_level'] if level is None else level)
    penalty,mechanism=protection(scene,layer,config,mode,level,ablation,trigger)
    out.mkdir(parents=True)
    identity=dict(input_manifest_sha256=layer['input_manifest_sha256'],
      layer1_sha256=sha(layer_path),config_sha256=sha(config),runner_sha256=sha(__file__),
      solver_sha256=sha(Path(__file__).with_name('native_solver.py')),
      mode=mode,level=level,ablation=ablation,trigger=trigger,power_cap=power_cap,
      groups=len(scene.resources))
    hs=[r.h_est for r in scene.resources]; masks=[r.serve for r in scene.resources]
    vs=solver._initial(scene.resources,hs,masks,scene.config.rank,scene.config.power_watt_per_trp)
    history=[dict(outer=0,**evaluate(scene,vs))]
    status='limit_unconverged';consecutive=0;wall=0.0;last_trace=[]
    limit=int(max_outer or cfg['max_outer'])
    for outer in range(1,limit+1):
        prior=vs;start=perf_counter()
        result=solver.solve(scene,'matrix_statistical',initial_precoders=vs,policy_covariances=penalty)
        wall+=perf_counter()-start;vs=result.precoders
        trace=result.trajectory;last_trace=trace
        inner=bool(len(trace)>1 and trace[-2].get('admm_converged',False))
        rel=trace[-1].get('outer_relative_objective_change')
        move=covariance_change(prior,vs)
        meets=inner and rel is not None and rel<=scene.config.outer_rel_tol and move<=scene.config.outer_precoder_rel_tol
        consecutive=consecutive+1 if meets else 0
        history.append(dict(outer=outer,designedObjective=result.optimization_objective,
          solverReason=result.reason,admmConverged=inner,outerRelativeChange=rel,
          covarianceRelativeChange=move,**evaluate(scene,vs)))
        if result.reason in ('optimization_objective_decrease','inexact_v_step_no_proxy_improvement'):
            status='failed';break
        if not inner:status='inner_inexact';break
        if consecutive>=cfg['consecutive_outer']:
            status='converged';break
    # Last ADMM trace has per-group CG records and primal/dual thresholds.
    state=dict(identity=identity,mechanism=mechanism,status=status,iter=len(history)-1,
      wallSeconds=wall,history=history,lastSolverTrace=last_trace,
      metric_definition='true channel weighted all-UE SE; violation proxy is SE < threshold')
    np.savez_compressed(out/'checkpoint.npz',state=json.dumps(state,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)),**{f'v{i}':v for i,v in enumerate(vs)})
    atomic_json(out/'state.json',state)
    print(json.dumps(dict(status=status,iter=state['iter'],mechanism=mechanism,
         final_true_sum_se=history[-1]['trueSumSE'],wallSeconds=wall),ensure_ascii=False))

def main():
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='command',required=True)
    c=s.add_parser('calibrate');c.add_argument('--input',required=True);c.add_argument('--config',required=True);c.add_argument('--out',required=True)
    r=s.add_parser('run');r.add_argument('--input',required=True);r.add_argument('--config',required=True)
    r.add_argument('--layer',required=True);r.add_argument('--mode',choices=['uniform','variable','fixed','mapped'],required=True)
    r.add_argument('--out',required=True);r.add_argument('--level',type=float);r.add_argument('--ablation',choices=['none','beam','service','power'],default='none')
    r.add_argument('--trigger',choices=['off','fixed','proposed'],default='proposed');r.add_argument('--max-outer',type=int)
    r.add_argument('--power-cap',type=float);r.add_argument('--max-groups',type=int)
    a=p.parse_args()
    if a.command=='calibrate':calibrate(a.input,a.config,a.out)
    else:run(a.input,a.config,a.layer,a.mode,a.out,a.level,a.ablation,a.trigger,a.max_outer,a.power_cap,a.max_groups)
if __name__=='__main__':main()
