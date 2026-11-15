"""Python server input entrance: import and verify MATLAB's frozen export."""
import argparse
import json
from pathlib import Path
from frozen_input import load


def main():
    p=argparse.ArgumentParser()
    p.add_argument('frozen_binary',help='directory exported by export_srs_joint_binary.m')
    p.add_argument('--report',required=True)
    a=p.parse_args()
    scene=load(a.frozen_binary)
    report=dict(modelId=scene.manifest['modelId'],
                inputSha256=scene.manifest['inputSha256'],
                configSha256=scene.manifest['configSha256'],
                users=scene.config.users,groups=len(scene.resources),
                scheduledOccurrences=sum(len(r.users) for r in scene.resources),
                srsSenders=scene.config.users,
                status='imported_frozen_matlab_input')
    target=Path(a.report)
    if target.exists(): raise FileExistsError(target)
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report))


if __name__=='__main__': main()
