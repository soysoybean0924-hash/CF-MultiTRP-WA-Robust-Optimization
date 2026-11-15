"""Independent true-channel metrics from frozen precoders; no solver state used."""
import argparse
import json
from pathlib import Path
import numpy as np
from frozen_input import load
from campaign import BRANCHES, identity


def metrics(scene, vs):
    ue=np.zeros(scene.config.users)
    power=np.zeros(scene.resources[0].h_true.shape[1])
    for r,v in zip(scene.resources,vs):
        power+=r.resource_weight*np.sum(abs(v)**2,axis=(1,2,3))
        h=r.h_true
        f=np.einsum('klrn,linq->kirq',h,v,optimize=True)
        for k,uid in enumerate(r.users):
            all_cov=np.eye(h.shape[2],dtype=complex)
            interference=all_cov.copy()
            for i in range(len(r.users)):
                block=f[k,i]@f[k,i].conj().T
                all_cov+=block
                if i!=k: interference+=block
            ue[uid]+=r.resource_weight*max(0.0,
                (np.linalg.slogdet(all_cov)[1]-
                 np.linalg.slogdet(interference)[1])/np.log(2))
    active=np.unique(np.concatenate([r.users for r in scene.resources]))
    floor=scene.config.pf_rate_floor_bits_per_s_hz
    return dict(perUserTrueSE=ue.tolist(),trueSumSE=float(ue.sum()),
                truePF=float(np.log1p(ue[active]/floor).sum()),
                trueRate5=float(np.percentile(ue,5)),powerByTrp=power.tolist(),
                maxPowerViolation=float(max(0,np.max(power)-scene.config.power_watt_per_trp)))


def main():
    p=argparse.ArgumentParser(); p.add_argument('frozen_input'); p.add_argument('checkpoints'); p.add_argument('output')
    p.add_argument('--admm-iterations',type=int); a=p.parse_args(); scene=load(a.frozen_input)
    rows={}; checkpoint_dir=Path(a.checkpoints)
    for branch in BRANCHES:
        path=checkpoint_dir/(branch+'.npz')
        if not path.exists(): rows[branch]=dict(status='missing'); continue
        with np.load(path,allow_pickle=False) as f:
            state=json.loads(str(f['state']))
            if state['identity']!=dict(identity(a.frozen_input,scene,branch), effectiveAdmmIterations=(a.admm_iterations if a.admm_iterations is not None else scene.config.admm_iterations)):
                raise ValueError(f'{branch}: identity mismatch')
            vs=[f[f'v{i}'] for i in range(len(scene.resources))]
        rows[branch]=dict(status=state['status'],outer=state['iter'],
                          final=metrics(scene,vs),history=state['history'])
    target=Path(a.output); target.parent.mkdir(parents=True,exist_ok=True)
    if target.exists(): raise FileExistsError(target)
    target.write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')

if __name__=='__main__': main()
