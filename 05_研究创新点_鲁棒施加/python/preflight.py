"""Server-only frozen-input and dependency audit; no optimization."""
import argparse
import json
from pathlib import Path
import numpy as np
from frozen_input import load
def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();cfg=json.loads(Path(a.config).read_text(encoding='utf-8'))
    scene=load(a.input)
    if len(scene.resources)!=17:raise ValueError('expected 17 resource groups')
    if any(r.source_cov is None for r in scene.resources):raise ValueError('source covariance absent')
    if any(r.error_var.shape[-2:]!=(4,4) for r in scene.resources):raise ValueError('physical C shape mismatch')
    if cfg['coverage_target']<=0 or cfg['coverage_target']>=1:raise ValueError('invalid coverage target')
    result=dict(status='preflight_ok',groups=len(scene.resources),users=scene.config.users,
      rank=scene.config.rank,physical_power_cap_w=scene.config.power_watt_per_trp,
      source_components=['residual','pilot','noise','age'],
      limitation='one frozen scene; no independent physical test seed or SRS quality sweep')
    out=Path(a.out);out.mkdir(parents=True,exist_ok=False)
    (out/'preflight.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result))
if __name__=='__main__':main()
