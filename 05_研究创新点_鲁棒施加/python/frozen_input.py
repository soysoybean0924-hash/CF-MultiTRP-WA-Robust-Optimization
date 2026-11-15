"""Reads MATLAB's explicit column-major binary export without MAT struct guesses."""
from dataclasses import dataclass
import json
import hashlib
from pathlib import Path
import numpy as np


@dataclass
class Resource:
    users: np.ndarray
    h_true: np.ndarray
    h_est: np.ndarray
    error_var: np.ndarray
    serve: np.ndarray
    resource_weight: float
    initial: np.ndarray
    source_cov: np.ndarray | None = None


@dataclass
class Config:
    users: int
    rank: int
    power_watt_per_trp: float
    outer_iterations: int
    outer_rel_tol: float
    outer_precoder_rel_tol: float
    admm_iterations: int
    cg_iterations: int
    eps_abs: float
    eps_rel: float
    pf_rate_floor_bits_per_s_hz: float
    solver: str = 'matrix_nominal'
    initialization_mode: str = 'svd'


@dataclass
class Scene:
    config: Config
    resources: list
    manifest: dict
    pilot_id: np.ndarray
    pilot_energy: np.ndarray
    prior_rx: np.ndarray


def _array(root, name, shape):
    re = np.fromfile(root / (name + '.real.bin'), dtype='<f8')
    if re.size != np.prod(shape):
        raise ValueError(f'{name}: size mismatch')
    im_path = root / (name + '.imag.bin')
    im = np.fromfile(im_path, dtype='<f8') if im_path.exists() else None
    if im is not None and im.size != re.size:
        raise ValueError(f'{name}: imaginary size mismatch')
    return (re if im is None else re+1j*im).reshape(shape, order='F')


def load(root):
    root = Path(root)
    m = json.loads((root / 'manifest.json').read_text(encoding='utf-8'))
    if m['modelId'] != 'srs_joint_snapshot_v2':
        raise ValueError('wrong model')
    for entry in m['binaryFiles']:
        p = root / entry['name']
        if hashlib.sha256(p.read_bytes()).hexdigest().lower() != entry['sha256'].lower():
            raise ValueError(f'input hash mismatch: {p.name}')
    c = m['config']
    cfg = Config(21*int(c['usersPerCell']), int(c['rank']),
                 float(c['powerWattPerTrp']), 1,
                 float(c['outerObjectiveRelativeTolerance']),
                 float(c['outerCovarianceRelativeTolerance']),
                 int(c['admmIterations']), int(c['cgIterations']),
                 float(c['admmAbsoluteTolerance']),
                 float(c['admmRelativeTolerance']), float(c['pfFloor']))
    resources = []
    for g, row in enumerate(m['resources'], 1):
        prefix = f'g{g:02d}_'
        k, l, nr, nt, rank = (int(row[x]) for x in
                              ('users', 'trps', 'rx', 'tx', 'rank'))
        shape = (k, l, nr, nt)
        source_path = root / (prefix+'sourceCov.real.bin')
        source_cov = (_array(root,prefix+'sourceCov',(k,l,nr,nr,4))
                      if source_path.exists() else None)
        resources.append(Resource(
            _array(root,prefix+'ueIds',(k,)).astype(int)-1,
            _array(root,prefix+'H_true',shape),
            _array(root,prefix+'H_hat',shape),
            _array(root,prefix+'C',(k,l,nr,nr)),
            _array(root,prefix+'serve',(l,k)).astype(bool),
            float(row['resourceWeightNu']),
            _array(root,prefix+'initial',(l,k,nt,rank)), source_cov))
    users=cfg.users; trps=resources[0].h_true.shape[1]
    return Scene(cfg,resources,m,
                 _array(root,'pilotId',(users,2)).astype(int)-1,
                 _array(root,'pilotEnergy',(users,2)),
                 _array(root,'priorRx',(4,4,users,trps)))
