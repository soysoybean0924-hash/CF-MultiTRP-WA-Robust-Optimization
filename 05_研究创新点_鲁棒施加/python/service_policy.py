"""One-replacement service policy with explicit uniform/source risk scores.

Development prototype. No test-channel data enters candidate selection.
"""
import numpy as np

def protected_service(resources, level, mode, margin, load_penalty, source_weights):
    if mode not in ('uniform','source') or level<0 or margin<0:
        raise ValueError('invalid service policy')
    weights=np.asarray(source_weights,float)
    if weights.shape!=(4,) or np.any(weights<0):raise ValueError('four nonnegative weights required')
    masks=[];diagnostics=[]
    for g,r in enumerate(resources,1):
        original=r.serve.copy();l_count,k_count=original.shape
        if np.any(original.sum(axis=0)!=2):raise ValueError('requires exactly two serving TRPs per scheduled UE')
        gain=np.sum(abs(r.h_est)**2,axis=(2,3)).real
        total=np.trace(r.error_var,axis1=2,axis2=3).real
        source=np.trace(r.source_cov,axis1=2,axis2=3).real
        if np.any(source<-1e-7):raise ValueError('negative source trace')
        source=np.maximum(source,0)
        # Risk is dimensionless. Uniform uses the same total-error rule on every link;
        # per-variable service uses pilot/age/noise/residual emphasis from configuration.
        numer=total if mode=='uniform' else source@weights
        risk=numer/np.maximum(gain+total,1e-30)
        signal=gain/np.maximum(gain+total,1e-30)
        load=original.sum(axis=1)/max(k_count,1)
        score=signal-level*risk-load_penalty*load[None,:]
        mask=original.copy();switches=0;candidate_gaps=[]
        for k in range(k_count):
            outgoing=min(np.flatnonzero(original[:,k]),key=lambda l:score[k,l])
            incoming=max(np.flatnonzero(~original[:,k]),key=lambda l:score[k,l])
            gap=float(score[k,incoming]-score[k,outgoing]);candidate_gaps.append(gap)
            if gap>margin:
                mask[outgoing,k]=False;mask[incoming,k]=True;switches+=1
        masks.append(mask)
        diagnostics.append(dict(group=g,switches=switches,
          candidate_gap_max=float(max(candidate_gaps,default=0)),
          candidate_gap_positive=int(sum(x>0 for x in candidate_gaps)),
          candidate_gap_above_margin=int(sum(x>margin for x in candidate_gaps))))
    return masks,diagnostics
