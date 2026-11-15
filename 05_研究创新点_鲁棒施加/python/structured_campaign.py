"""Restartable strategy campaign. Frozen physical C remains immutable."""
import argparse
import json
import os
from pathlib import Path
from time import perf_counter
import numpy as np
from frozen_input import load
import native_solver as solver
from campaign import sha, covariance_change, evaluate
from source_policy import classify_links, policy_covariances
from structured_service import select_services_limited

VARIANTS=('nominal','robust','uniform','differentiated',
          'differentiated_service','differentiated_service_graph')


def prepare(scene, variant, config):
    if variant not in VARIANTS:
        raise ValueError(variant)
    labels=[]; statistics=[]
    for r in scene.resources:
        lab,stat=classify_links(r,scene.prior_rx,**{
            'min_relative':config['classification']['minRelative'],
            'dominance':config['classification']['dominance'],
            'margin':config['classification']['margin'],
            'residual_limit':config['classification']['residualLimit']})
        labels.append(lab);statistics.append(stat)
    original=[r.serve.copy() for r in scene.resources]
    if variant in ('differentiated_service','differentiated_service_graph'):
        c=config['service']
        selected=select_services_limited(scene.resources,pilot_penalty=c['pilotPenalty'],
            load_penalty=c['loadPenalty'],switch_margin=c['switchMargin'])
        for r,mask in zip(scene.resources,selected):
            r.serve=mask
    changed=sum(np.count_nonzero(r.serve!=old) for r,old in zip(scene.resources,original))
    policy=None
    if variant in ('uniform','differentiated','differentiated_service',
                   'differentiated_service_graph'):
        mode='uniform' if variant=='uniform' else 'differentiated'
        design_gap=[]
        for r in scene.resources:
            cold=solver._initial([r],[r.h_est],[r.serve],scene.config.rank,
                                  scene.config.power_watt_per_trp)[0]
            estimated=solver.rates(r.h_est,cold)*r.resource_weight
            floor=float(config.get('designRateReference',0.05))
            design_gap.append(np.clip((floor-estimated)/floor,0,1))
        policy=policy_covariances(scene.resources,statistics,config['strength'],mode,
                                  config['sourceWeights'],design_gap=design_gap)
    native='matrix_nominal' if variant=='nominal' else 'matrix_statistical'
    return native,policy,dict(changedLinks=int(changed),
        labelCounts={name:int(sum(np.count_nonzero(x==name) for x in labels))
                     for name in sorted(set(str(v) for x in labels for v in x.flat))},
        policyTrace=float(sum(np.trace(p,axis1=2,axis2=3).real.sum() for p in policy))
                    if policy is not None else 0.0)


def run(input_root, output, config_path, variant, max_groups=None,pilot_input=None):
    config=json.loads(Path(config_path).read_text(encoding='utf-8'))
    if variant not in config['variants'] or variant not in VARIANTS:
        raise ValueError('variant not enabled')
    base_root=input_root
    if variant=='differentiated_service_graph':
        if pilot_input is None:
            raise ValueError('graph branch requires regenerated --pilot-input')
        input_root=pilot_input
    scene=load(input_root)
    if max_groups is not None:
        if not 1<=max_groups<=len(scene.resources):
            raise ValueError('invalid group count')
        scene.resources=scene.resources[:max_groups]
    if variant=='differentiated_service_graph':
        base=load(base_root)
        base.resources=base.resources[:len(scene.resources)]
        for old,new in zip(base.resources,scene.resources):
            if not np.array_equal(old.users,new.users) or not np.allclose(old.h_true,new.h_true,atol=0,rtol=0):
                raise ValueError('graph candidate does not share true channel and schedule')
    native,policy,diagnostic=prepare(scene,variant,config)
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    checkpoint=output/(variant+'.npz')
    identity=dict(manifestSha256=sha(Path(input_root)/'manifest.json'),
        configSha256=sha(config_path),runnerSha256=sha(__file__),
        solverSha256=sha(Path(__file__).with_name('native_solver.py')),
        policySha256=sha(Path(__file__).with_name('source_policy.py')),
        serviceSha256=sha(Path(__file__).with_name('structured_service.py')),
        inputSha256=scene.manifest['inputSha256'],variant=variant,
        baseManifestSha256=sha(Path(base_root)/'manifest.json'),
        groups=len(scene.resources))
    if checkpoint.exists():
        with np.load(checkpoint,allow_pickle=False) as f:
            state=json.loads(str(f['state']))
            vs=[f[f'v{i}'].copy() for i in range(len(scene.resources))]
        if state['identity']!=identity:
            raise ValueError('checkpoint identity mismatch')
    else:
        hs=[r.h_est for r in scene.resources]
        masks=[r.serve for r in scene.resources]
        vs=solver._initial(scene.resources,hs,masks,scene.config.rank,
                           scene.config.power_watt_per_trp)
        state=dict(identity=identity,status='initial',iter=0,consecutive=0,
                   converged=False,wallSeconds=0.0,diagnostic=diagnostic,
                   history=[dict(outer=0,**evaluate(scene,vs))])
        save(checkpoint,state,vs)
    if state['converged'] or state['status'] in ('failed','inner_inexact','limit_unconverged'):
        return state
    for _ in range(config['batchOuter']):
        if state['iter']>=config['maxOuter']:
            state['status']='limit_unconverged';break
        started=perf_counter()
        old=vs
        result=solver.solve(scene,native,initial_precoders=vs,
                            policy_covariances=policy)
        vs=result.precoders
        duration=perf_counter()-started
        trajectory=result.trajectory
        last=trajectory[-1]
        previous=trajectory[-2] if len(trajectory)>1 else {}
        rel=last.get('outer_relative_objective_change')
        qchange=covariance_change(old,vs)
        inner=bool(previous.get('admm_converged',False))
        meets=(rel is not None and rel<=scene.config.outer_rel_tol and
               qchange<=scene.config.outer_precoder_rel_tol and inner)
        state['consecutive']=state['consecutive']+1 if meets else 0
        state['iter']+=1
        state['wallSeconds']+=duration
        state['converged']=state['consecutive']>=config['consecutiveOuterSteps']
        state['status']=('converged' if state['converged'] else
            'failed' if result.reason in ('optimization_objective_decrease',
                              'inexact_v_step_no_proxy_improvement') else
            'inner_inexact' if not inner else
            'limit_unconverged' if state['iter']>=config['maxOuter'] else 'batch_limit')
        state['history'].append(dict(outer=state['iter'],wallSeconds=duration,
            designedObjective=result.optimization_objective,
            solverReason=result.reason,outerRelativeChange=rel,
            covarianceRelativeChange=qchange,admmConverged=inner,
            admmIterations=previous.get('admm_iterations'),**evaluate(scene,vs)))
        save(checkpoint,state,vs)
        if state['status']!='batch_limit':break
    return state


def save(path,state,vs):
    partial=path.with_suffix('.partial')
    if partial.exists():raise FileExistsError(partial)
    with partial.open('wb') as f:
        np.savez_compressed(f,state=json.dumps(state),
                            **{f'v{i}':v for i,v in enumerate(vs)})
        f.flush();os.fsync(f.fileno())
    os.replace(partial,path)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('frozen_input');p.add_argument('output')
    p.add_argument('--config',required=True)
    p.add_argument('--variant',choices=VARIANTS)
    p.add_argument('--max-groups',type=int)
    p.add_argument('--pilot-input')
    p.add_argument('--until-stop',action='store_true')
    a=p.parse_args()
    config=json.loads(Path(a.config).read_text(encoding='utf-8'))
    for variant in ([a.variant] if a.variant else config['variants']):
        while True:
            state=run(a.frozen_input,a.output,a.config,variant,a.max_groups,a.pilot_input)
            print(variant,state['iter'],state['status'],flush=True)
            if state['status'] in ('failed','inner_inexact'):
                raise RuntimeError(f'{variant}: {state["status"]}')
            if not a.until_stop or state['status']!='batch_limit':break


if __name__=='__main__':main()
