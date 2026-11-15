"""Exact array-level port of pilot/srs_joint_snapshot.m; arrays are supplied, not seeded."""
import numpy as np

MODEL_ID = "srs_joint_snapshot_v2"


def snapshot(h0, r, pilot_id, energy, noise, z, age_ms, tc_ms, innovation,
             port_rows=(0, 1), save_cross=True, observation_events=None):
    """Finite snapshot; optional four rows require two explicit 2-port events.

    The exponential age factor is an alternative model, not a Jakes trajectory.
    """
    if h0.ndim != 4:
        raise ValueError(f'h0 expected [K,L,4,Nt], actual shape {h0.shape}')
    k_count, l_count, nr, nt = h0.shape
    if min(k_count,l_count,nt) < 1 or nr != 4:
        raise ValueError(f'h0 expected K>=1,L>=1,Nr=4,Nt>=1, actual shape {h0.shape}')
    n_ports=len(port_rows)
    if n_ports not in (2,4) or pilot_id.ndim != 2 or pilot_id.shape != (k_count,n_ports):
        raise ValueError(f'pilot_id expected {(k_count,n_ports)}, actual shape {pilot_id.shape}')
    if n_ports==4:
        if observation_events is None or len(observation_events)!=2:
            raise ValueError('four rows require two explicit 2-port observation events')
        if any(len(event['rows'])!=2 or 'time_ms' not in event or
               'resource' not in event or 'energy' not in event
               for event in observation_events):
            raise ValueError('each switch event needs rows, time_ms, resource, energy')
        if tuple(observation_events[0]['rows']+observation_events[1]['rows'])!=tuple(port_rows):
            raise ValueError('switch event rows must match observation columns')
    if z.ndim != 3:
        raise ValueError(f'z expected [Nt,L,T], actual shape {z.shape}')
    t_count = z.shape[2]
    for name, value, expected in (
            ('r', r, (nr, nr, k_count, l_count)),
            ('energy', energy, (k_count, n_ports)),
            ('z', z, (nt, l_count, t_count)),
            ('innovation', innovation, (k_count, l_count, nr, nt))):
        if value.shape != expected:
            raise ValueError(f'{name} expected {expected}, actual shape {value.shape}')
    if len(set(port_rows)) != n_ports or any(row < 0 or row >= nr for row in port_rows):
        raise ValueError(f'port_rows expected distinct indices in [0,{nr-1}], actual {port_rows}')
    for name, value, positive in (('noise', noise, True),
                                  ('tc_ms', tc_ms, True),
                                  ('age_ms', age_ms, False)):
        if np.ndim(value) != 0 or not np.isrealobj(value) or not np.isfinite(value) or (value <= 0 if positive else value < 0):
            raise ValueError(f'{name} must be a finite real scalar with ' +
                             ('positive' if positive else 'nonnegative') + f' value; actual shape {np.shape(value)}')
    if not np.all(np.isfinite(pilot_id)) or not np.all((pilot_id >= 0) & (pilot_id < t_count)) or not np.all(pilot_id == np.floor(pilot_id)):
        raise ValueError(f'pilot_id must contain integer indices in [0,{t_count-1}]')
    if not np.isrealobj(energy) or not np.all(np.isfinite(energy)) or not np.all(energy >= 0):
        raise ValueError('energy must be real, finite, nonnegative')
    a = np.exp(-age_ms / tc_ms)
    h_true = a * h0 + np.sqrt(1 - a*a) * innovation
    h_hat = np.zeros_like(h0, dtype=complex)
    c = np.zeros((nr, nr, k_count, l_count), complex)
    cross = np.zeros((nr, nr, k_count, k_count, l_count), complex) if save_cross else None
    y_all = np.zeros((t_count, nt, l_count), complex)
    for l in range(l_count):
        m = np.zeros((t_count, k_count*nr), complex)
        q = np.zeros((k_count*nr, k_count*nr), complex)
        x = np.zeros((k_count*nr, nt), complex)
        for k in range(k_count):
            ix = slice(k*nr, (k+1)*nr)
            x[ix] = h0[k, l]
            q[ix, ix] = r[:, :, k, l]
            for p, row in enumerate(port_rows):
                m[pilot_id[k, p], k*nr+row] += np.sqrt(energy[k, p])
        y = m @ x + np.sqrt(noise) * z[:, l, :].T
        s = m @ q @ m.conj().T + noise*np.eye(t_count)
        b = np.linalg.solve(s.T, (q @ m.conj().T).T).T
        xhat = b @ y
        c0 = q - b @ m @ q
        c0 = (c0+c0.conj().T)/2
        y_all[:, :, l] = y
        for k in range(k_count):
            ix = slice(k*nr, (k+1)*nr)
            h_hat[k, l] = a*xhat[ix]
            c[:, :, k, l] = a*a*c0[ix, ix] + (1-a*a)*r[:, :, k, l]
            for j in range(k_count) if save_cross else ():
                jx = slice(j*nr, (j+1)*nr)
                cross[:, :, k, j, l] = a*a*c0[ix, jx]
                if j == k:
                    cross[:, :, k, j, l] += (1-a*a)*r[:, :, k, l]
    return dict(modelId=MODEL_ID, H_true=h_true, H_hat=h_hat, C=c,
                crossUE=cross, Y=y_all, a=a,
                observationEvents=observation_events,
                observationEnergy=(float(np.sum(energy)) if observation_events else None))
