"""Six fixed-input branches; one WMMSE outer update per restartable batch step."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import numpy as np
from frozen_input import load
import native_solver as solver

BRANCHES = {
    'ideal_sumse':'ideal', 'nominal_sumse':'matrix_nominal',
    'robust_sumse':'matrix_statistical', 'ideal_pf':'ideal_pf',
    'nominal_pf':'matrix_pf_nominal', 'robust_pf':'matrix_pf_statistical'}
VERSION = 'srs_joint_python_campaign_v3_cg_certificate'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def identity(root, scene, branch):
    return dict(version=VERSION, modelId=scene.manifest['modelId'],
                configSha256=scene.manifest['configSha256'],
                inputSha256=scene.manifest['inputSha256'], branch=branch,
                solverSha256=sha(Path(__file__).with_name('native_solver.py')),
                runnerSha256=sha(__file__),
                manifestSha256=sha(Path(root)/'manifest.json'))


def evaluate(scene, vs):
    ue=np.zeros(scene.config.users)
    for r,v in zip(scene.resources,vs):
        ue[r.users]+=r.resource_weight*solver.rates(r.h_true,v)
    scheduled=np.unique(np.concatenate([r.users for r in scene.resources]))
    floor=scene.config.pf_rate_floor_bits_per_s_hz
    return dict(perUserTrueSE=ue.tolist(), trueSumSE=float(ue.sum()),
                truePF=float(np.log1p(ue[scheduled]/floor).sum()),
                trueRate5=float(np.percentile(ue,5)),
                powerByTrp=solver.power_by_trp(scene.resources,vs).tolist())


def covariance_change(old, new):
    numerator=denominator=0.0
    for x,y in zip(old,new):
        for l in range(x.shape[0]):
            xflat=x[l].transpose(1,0,2).reshape(x.shape[2],-1)
            yflat=y[l].transpose(1,0,2).reshape(y.shape[2],-1)
            q0=xflat @ xflat.conj().T
            q1=yflat @ yflat.conj().T
            numerator+=float(np.vdot(q1-q0,q1-q0).real)
            denominator+=float(np.vdot(q0,q0).real)
    return np.sqrt(numerator)/max(1.0,np.sqrt(denominator))


def _save(path, state, vs):
    tmp=path.with_suffix('.partial')
    if tmp.exists():
        raise FileExistsError(f'inspect interrupted checkpoint: {tmp}')
    with tmp.open('wb') as fh:
        np.savez_compressed(fh, state=json.dumps(state),
                            **{f'v{i}':v for i,v in enumerate(vs)})
        fh.flush(); os.fsync(fh.fileno())
    os.replace(tmp,path)


def run(root, out, branch, smoke_outer=None, admm_iterations=None):
    scene=load(root)
    if admm_iterations is not None:
        scene.config.admm_iterations=int(admm_iterations)
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    checkpoint=out/(branch+'.npz')
    ident=identity(root,scene,branch)
    ident['effectiveAdmmIterations']=scene.config.admm_iterations
    mode='smoke' if smoke_outer is not None else 'formal'
    if checkpoint.exists():
        with np.load(checkpoint,allow_pickle=False) as f:
            state=json.loads(str(f['state']))
            vs=[f[f'v{i}'].copy() for i in range(len(scene.resources))]
        if state['identity'] != ident or state.get('runMode')!=mode:
            raise ValueError('checkpoint model/config/input/branch/code mismatch')
    else:
        vs=[r.initial.copy() for r in scene.resources]
        state=dict(identity=ident,runMode=mode,iter=0,consecutive=0,converged=False,
                   status='not_run',history=[dict(outer=0, status='initial',
                   **evaluate(scene,vs))])
        _save(checkpoint,state,vs)
    if state['converged']:
        return state
    native=BRANCHES[branch]
    cfg=scene.config
    max_outer=(int(smoke_outer) if smoke_outer is not None else
               int(scene.manifest['config']['maxOuter']))
    target=int(scene.manifest['config']['consecutiveOuterSteps'])
    batch=int(scene.manifest['config']['batchOuter'])
    for _ in range(batch):
        if state['iter'] >= max_outer:
            state['status']='limit_unconverged'; break
        old_vs=vs
        result=solver.solve(scene,native,initial_precoders=vs)
        vs=result.precoders
        qchange=covariance_change(old_vs,vs)
        actual=evaluate(scene,vs)
        last=result.trajectory[-1]
        rel=last.get('outer_relative_objective_change')
        prev=result.trajectory[-2] if len(result.trajectory)>1 else {}
        inner=bool(prev.get('admm_converged',False))
        meets=(rel is not None and
               rel<=cfg.outer_rel_tol and qchange<=cfg.outer_precoder_rel_tol
               and inner)
        state['consecutive']=state['consecutive']+1 if meets else 0
        state['iter']+=1
        state['history'].append(dict(outer=state['iter'],
            designedObjective=result.optimization_objective,
            outerRelativeChange=rel, covarianceRelativeChange=qchange,
            acceptedStep=prev.get('accepted_step'),
            solverReason=result.reason,
            admmConverged=inner, admmTrace=prev.get('admm_trace',[]),
            **actual))
        state['converged']=state['consecutive']>=target
        state['status']=('converged' if state['converged'] else
                         'failed' if result.reason in ('optimization_objective_decrease','inexact_v_step_no_proxy_improvement') else
                         'inner_inexact' if not inner else
                         'limit_unconverged' if state['iter']>=max_outer else
                         'batch_limit')
        _save(checkpoint,state,vs)
        if state['status'] in ('failed','inner_inexact'): break
        if state['converged']: break
    return state


def main():
    p=argparse.ArgumentParser()
    p.add_argument('frozen_input'); p.add_argument('output')
    p.add_argument('--branch',choices=list(BRANCHES)+['all'],default='all')
    p.add_argument('--smoke-outer',type=int)
    p.add_argument('--admm-iterations',type=int)
    p.add_argument('--until-stop',action='store_true')
    a=p.parse_args()
    if a.admm_iterations is not None and a.admm_iterations<1:
        p.error('--admm-iterations must be positive')
    if a.smoke_outer is not None and a.smoke_outer<1:
        p.error('--smoke-outer must be positive')
    for b in (BRANCHES if a.branch=='all' else [a.branch]):
        while True:
            s=run(a.frozen_input,a.output,b,a.smoke_outer,a.admm_iterations)
            print(b,s['iter'],s['status'],flush=True)
            if s['status'] in ('failed','inner_inexact'):
                raise RuntimeError(f'{b}: {s["status"]}; stop before next branch')
            if not a.until_stop or s['status']!='batch_limit': break


if __name__=='__main__': main()
