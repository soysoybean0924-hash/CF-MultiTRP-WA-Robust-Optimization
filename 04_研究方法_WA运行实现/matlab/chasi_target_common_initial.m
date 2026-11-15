function [precoder,power]=chasi_target_common_initial(resources,rank,budget,seed)
%CHASI_TARGET_COMMON_INITIAL CSI-independent feasible native target input.
% Per-TRP power is sum_g nu_g sum_{k,d} |x_{l,k,:,d}|^2.
assert(rank>=1 && rank==floor(rank) && budget>0);
assert(seed>=0 && seed==floor(seed));
G=numel(resources); assert(G>=1);
L=size(resources(1).H_true,2);
Nt=size(resources(1).H_true,4);
stream=RandStream('mt19937ar','Seed',seed);
precoder=cell(1,G); power=zeros(L,1);
for g=1:G
    r=resources(g); K=size(r.H_true,1);
    assert(size(r.H_true,2)==L && size(r.H_true,4)==Nt);
    assert(isequal(size(r.servingMask),[L,K]));
    nu=double(r.resourceWeightNu);
    assert(isfinite(nu) && nu>=0);
    x=(randn(stream,L,K,Nt,rank)+ ...
        1i*randn(stream,L,K,Nt,rank))/sqrt(2);
    x=x.*reshape(r.servingMask,L,K,1,1);
    precoder{g}=x;
    power=power+nu*reshape(sum(abs(x).^2,[2,3,4]),L,1);
end
factor=sqrt((budget/2)./max(power,realmin));
power(:)=0;
for g=1:G
    r=resources(g);
    precoder{g}=precoder{g}.*reshape(factor,L,1,1,1);
    power=power+double(r.resourceWeightNu)* ...
        reshape(sum(abs(precoder{g}).^2,[2,3,4]),L,1);
end
assert(all(power<=budget/2+1e-8));
end
