"""Server-only one-step reproduction to diagnose old nominal inner accuracy.

Does not edit the historical checkpoint or relax solver tolerances.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from frozen_input import load
import native_solver

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--out',required=True)
    p.add_argument('--admm-iterations',type=int,default=100)
    a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=False)
    scene=load(a.input)
    if a.admm_iterations<1:raise ValueError('admm iterations must be positive')
    scene.config.admm_iterations=a.admm_iterations
    vs=native_solver._initial(scene.resources,[r.h_est for r in scene.resources],
                              [r.serve for r in scene.resources],scene.config.rank,
                              scene.config.power_watt_per_trp)
    result=native_solver.solve(scene,'matrix_nominal',initial_precoders=vs)
    trace=result.trajectory
    rows=[]
    for outer in trace:
        for inner in outer.get('admm_trace',[]):
            rows.append(dict(outer=outer['outer'],inner=inner['inner'],
              primal=inner['primal_residual'],primal_limit=inner['eps_primal'],
              dual=inner['dual_residual'],dual_limit=inner['eps_dual'],
              cg_all_accurate=inner['cg_all_accurate'],
              cg_records=inner['cg_records']))
    (out/'nominal_admm.json').write_text(json.dumps(dict(reason=result.reason,
      outer_iterations=result.outer_iterations,admm_iteration_cap=a.admm_iterations,
      eps_abs=scene.config.eps_abs,eps_rel=scene.config.eps_rel,rows=rows),indent=2,
      default=lambda x:x.item() if isinstance(x,np.generic) else str(x)),encoding='utf-8')
    print(json.dumps(dict(reason=result.reason,outer_iterations=result.outer_iterations,
      last_admm=rows[-1] if rows else None),default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
if __name__=='__main__':main()
