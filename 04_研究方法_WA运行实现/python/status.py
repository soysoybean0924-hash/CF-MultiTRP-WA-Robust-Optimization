"""Read checkpoint completion states without starting optimization."""
import argparse
import json
from pathlib import Path
import numpy as np
from campaign import BRANCHES


def main():
    p=argparse.ArgumentParser();p.add_argument('checkpoints');a=p.parse_args()
    root=Path(a.checkpoints)
    for branch in BRANCHES:
        path=root/(branch+'.npz')
        if not path.exists(): print(branch,'missing');continue
        with np.load(path,allow_pickle=False) as f:
            state=json.loads(str(f['state']))
        print(branch,state['iter'],state['status'],
              f"consecutive={state['consecutive']}")


if __name__=='__main__':main()
