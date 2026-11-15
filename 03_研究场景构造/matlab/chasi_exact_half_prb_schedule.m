function audit=chasi_exact_half_prb_schedule(groupRb,slots)
% Integer PRBs: 136/137 alternating gives exactly 50% of 273 over even T.
groupRb=double(groupRb(:));assert(sum(groupRb)==273 && ...
    all(groupRb>0) && all(groupRb==floor(groupRb)));
assert(slots>0 && slots==floor(slots) && mod(slots,2)==0);
G=numel(groupRb);allocation=zeros(G,slots);
for t=1:slots
    target=136+mod(t-1,2);remain=target;
    order=mod((0:G-1)+(t-1),G)+1;
    for g=order
        take=min(groupRb(g),remain);allocation(g,t)=take;
        remain=remain-take;if remain==0,break;end
    end
    assert(remain==0);
end
used=sum(allocation,1);
assert(all(used(1:2:end)==136) && all(used(2:2:end)==137));
audit=struct('requestedFraction',.5,'allocatedRbByGroupTti',allocation, ...
    'usedRbByTti',used,'averageFraction',mean(used)/273, ...
    'perTrpFractionByTti',used/273, ...
    'rule','common all-TRP rotating group order, alternate 136/137 RB');
end
