"""One slow-layer service update using estimated channels and source statistics."""
import numpy as np


def select_services(resources, minimum=2, pilot_penalty=0.25):
    if minimum < 1 or not 0 <= pilot_penalty <= 1:
        raise ValueError('invalid service policy')
    selections = []
    for resource in resources:
        if resource.source_cov is None:
            raise ValueError('source covariance is required')
        k_count, l_count = resource.h_est.shape[:2]
        if minimum > l_count:
            raise ValueError('minimum exceeds TRP count')
        gain = np.sum(abs(resource.h_est)**2, axis=(2, 3)).real
        total = np.trace(resource.error_var, axis1=2, axis2=3).real
        pilot = np.trace(resource.source_cov[..., 1], axis1=2, axis2=3).real
        score = (gain-pilot_penalty*pilot)/np.maximum(gain+total, 1e-30)
        mask = np.zeros((l_count, k_count), bool)
        for k in range(k_count):
            selected = np.argsort(-score[k], kind='stable')[:minimum]
            mask[selected, k] = True
        selections.append(mask)
    return selections


def select_services_limited(resources, minimum=2, pilot_penalty=0.25,
                            load_penalty=0.05, switch_margin=0.03):
    """At most one equal-budget replacement per UE, using available statistics."""
    selections=[]
    for resource in resources:
        if resource.source_cov is None:
            raise ValueError('source covariance required')
        old=resource.serve.copy()
        l_count,k_count=old.shape
        if np.any(old.sum(axis=0)<minimum):
            raise ValueError('old service violates minimum coverage')
        gain=np.sum(abs(resource.h_est)**2,axis=(2,3)).real
        total=np.trace(resource.error_var,axis1=2,axis2=3).real
        pilot=np.trace(resource.source_cov[...,1],axis1=2,axis2=3).real
        load=old.sum(axis=1)/max(k_count,1)
        score=(gain-pilot_penalty*pilot)/np.maximum(gain+total,1e-30)
        score-=load_penalty*load[None,:]
        for k in range(k_count):
            if old[:,k].sum()!=minimum:
                continue
            out=min(np.flatnonzero(old[:,k]),key=lambda l:score[k,l])
            incoming=max(np.flatnonzero(~old[:,k]),key=lambda l:score[k,l])
            if score[k,incoming]>score[k,out]+switch_margin:
                old[out,k]=False
                old[incoming,k]=True
        selections.append(old)
    return selections
