"""Descriptive error labels and bounded design regularization; physical C is immutable."""
import numpy as np

SOURCE_NAMES = ('residual', 'pilot', 'noise', 'age')
LABELS = ('PILOT_DOMINANT', 'NOISE_DOMINANT', 'AGE_DOMINANT')


def classify_links(resource, prior_rx, min_relative=0.05,
                   dominance=0.55, margin=0.15, residual_limit=0.45):
    if resource.source_cov is None:
        raise ValueError('source covariance required')
    source=np.trace(resource.source_cov,axis1=2,axis2=3).real
    if np.min(source)<-1e-7:
        raise ValueError('negative source variance')
    source=np.maximum(source,0)
    total=source.sum(axis=-1)
    prior=np.trace(prior_rx[:,:,resource.users,:],axis1=0,axis2=1).real
    relative=total/np.maximum(prior,1e-30)
    shares=source/np.maximum(total[...,None],1e-30)
    physical=shares[...,1:]
    ordered=np.sort(physical,axis=-1)
    gap=ordered[...,-1]-ordered[...,-2]
    labels=np.full(total.shape,'MIXED',dtype='<U31')
    labels[relative<min_relative]='LOW_ABSOLUTE_ERROR'
    residual=(shares[...,0]>=residual_limit) & (relative>=min_relative)
    labels[residual]='RESIDUAL_DOMINANT_OR_UNCERTAIN'
    dominant=(physical.max(axis=-1)>=dominance)&(gap>=margin)&~residual&(relative>=min_relative)
    mapped=np.asarray(LABELS)[physical.argmax(axis=-1)]
    labels[dominant]=mapped[dominant]
    return labels,dict(trace=source,share=shares,relative=relative,gap=gap,
                       total=total)


def policy_covariances(resources, statistics, strength, mode,
                       source_weights=(0.1,1.0,0.2,0.6), relative_scale=0.1,
                       design_gap=None):
    """PSD penalty in design MSE only. Uniform matches total extra trace budget."""
    if not 0<=strength<=0.3 or mode not in ('uniform','differentiated'):
        raise ValueError('strength or mode outside declared range')
    weights=np.asarray(source_weights,float)
    if weights.shape!=(4,) or np.any(weights<0):
        raise ValueError('nonnegative four-source weights required')
    extra=[]
    for index,(resource,stat) in enumerate(zip(resources,statistics)):
        gap=(np.ones_like(stat['relative']) if design_gap is None else
             np.asarray(design_gap[index],float)[:,None])
        if np.any(~np.isfinite(gap)) or np.any((gap<0)|(gap>1)):
            raise ValueError('design gap must lie in [0,1]')
        confidence=np.clip(1-stat['share'][...,0],0,1)
        gate=np.clip(stat['relative']/relative_scale,0,1)*confidence*(0.5+0.5*gap)
        extra.append(strength*gate[...,None,None]*np.einsum(
            'klrsq,q->klrs',resource.source_cov,weights,optimize=True))
    if mode=='differentiated':
        return extra
    extra_total=sum(float(np.trace(x,axis1=2,axis2=3).real.sum()) for x in extra)
    physical_total=sum(float(np.trace(r.error_var,axis1=2,axis2=3).real.sum())
                       for r in resources)
    factor=extra_total/max(physical_total,1e-30)
    return [factor*r.error_var for r in resources]
