"""Server-only read-only audit of paired smoke mechanism activity."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from frozen_input import load
from layered import protection

def inspect(state_path,input_root,layer,config):
    state=json.loads(Path(state_path).read_text(encoding='utf-8'))
    ck=Path(state_path).with_name('checkpoint.npz')
    with np.load(ck,allow_pickle=False) as z:
        vs=[z[f'v{i}'].copy() for i in range(state['identity']['groups'])]
    scene=load(input_root)
    scene.resources=scene.resources[:state['identity']['groups']]
    layer_data=json.loads(Path(layer).read_text(encoding='utf-8'))
    old_masks=[r.serve.copy() for r in scene.resources]
    _,mechanism=protection(scene,layer_data,config,state['identity']['mode'],
        state['identity']['level'],state['identity']['ablation'],state['identity']['trigger'])
    masks=[r.serve.copy() for r in scene.resources]
    return state,vs,masks,dict(
      configured=mechanism,
      changed_mask_elements=int(sum(np.count_nonzero(a!=b) for a,b in zip(old_masks,masks))),
      mask_sha256=hashlib.sha256(b''.join(m.astype('u1').tobytes() for m in masks)).hexdigest(),
      final_total_power_w=float(sum(state['history'][-1]['powerByTrp'])),
      final_max_trp_power_w=float(max(state['history'][-1]['powerByTrp'])),
      power_cap_slack_w=float(mechanism['effective_cap_w']-max(state['history'][-1]['powerByTrp'])),
      final_true_sum_se=float(state['history'][-1]['trueSumSE']),
      status=state['status'],iter=state['iter'])

def main():
    p=argparse.ArgumentParser()
    for key in ('input','layer','config','uniform','variable','out'):p.add_argument('--'+key,required=True)
    a=p.parse_args()
    u,uv,um,ur=inspect(a.uniform,a.input,a.layer,a.config)
    v,vv,vm,vr=inspect(a.variable,a.input,a.layer,a.config)
    if len(uv)!=len(vv):raise ValueError('different group count')
    result=dict(uniform=ur,variable=vr,
      same_service_masks=all(np.array_equal(x,y) for x,y in zip(um,vm)),
      changed_mask_elements_between_groups=int(sum(np.count_nonzero(x!=y) for x,y in zip(um,vm))),
      precoder_relative_difference=float(np.sqrt(sum(np.vdot(x-y,x-y).real for x,y in zip(uv,vv)))/
        max(1.0,np.sqrt(sum(np.vdot(x,x).real for x in uv)))),
      user_se_max_absolute_difference=float(np.max(np.abs(
        np.asarray(u['history'][-1]['perUserTrueSE'])-np.asarray(v['history'][-1]['perUserTrueSE'])))))
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    Path(a.out).write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
