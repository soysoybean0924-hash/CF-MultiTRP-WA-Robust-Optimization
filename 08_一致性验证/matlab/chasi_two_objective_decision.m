function out=chasi_two_objective_decision(base,ideal,method,opts)
% Separate objective decisions on identical nonmissing UE-second samples.
required={'rateMbps','valid','totalAirPowerW'};
for j=1:numel(required)
 f=required{j};assert(isfield(base,f)&&isfield(ideal,f)&&isfield(method,f));
end
same=isequal(size(base.rateMbps),size(ideal.rateMbps),size(method.rateMbps));
assert(same && isequal(size(base.valid),size(base.rateMbps)) && ...
    isequal(size(ideal.valid),size(base.rateMbps)) && ...
    isequal(size(method.valid),size(base.rateMbps)));
assert(all(isfinite([base.totalAirPowerW ideal.totalAirPowerW ...
    method.totalAirPowerW])));
mask=base.valid & ideal.valid & method.valid;
required=true(size(mask));
if isfield(opts,'requiredMask')
 assert(isequal(size(opts.requiredMask),size(mask)));
 required=logical(opts.requiredMask);
 mask=mask & required;
end
out=struct('commonSampleCount',nnz(mask),'sampleMask',mask, ...
    'missingRequiredCount',nnz(required & ~mask), ...
    'missingBaseCount',nnz(required & ~base.valid), ...
    'missingIdealCount',nnz(required & ~ideal.valid), ...
    'missingMethodCount',nnz(required & ~method.valid), ...
    'edgeStatus','missing_samples','powerStatus','missing_samples', ...
    'edgeRatio',NaN,'meanLoss',NaN,'powerReduction',NaN, ...
    'r5BaseMbps',NaN,'r5IdealMbps',NaN,'r5MethodMbps',NaN, ...
    'r5GainDenominatorMbps',NaN, ...
    'meanBaseMbps',NaN,'meanMethodMbps',NaN, ...
    'baseAirPowerW',base.totalAirPowerW, ...
    'methodAirPowerW',method.totalAirPowerW, ...
    'requiredSampleCount',nnz(required));
if ~any(mask) || any(required(:) & ~mask(:)),return;end
b=base.rateMbps(mask);i=ideal.rateMbps(mask);m=method.rateMbps(mask);
assert(all(isfinite([b;i;m])) && all([b;i;m]>=0));
pb=sort(b);pi=sort(i);pm=sort(m);idx=max(1,ceil(.05*numel(b)));
den=pi(idx)-pb(idx);
out.r5BaseMbps=pb(idx);out.r5IdealMbps=pi(idx);
out.r5MethodMbps=pm(idx);out.r5GainDenominatorMbps=den;
out.meanBaseMbps=mean(b);out.meanMethodMbps=mean(m);
if den>0
 out.edgeRatio=(pm(idx)-pb(idx))/den;
 if out.edgeRatio>=.8,out.edgeStatus='pass';else,out.edgeStatus='fail';end
else,out.edgeStatus='nonpositive_denominator';end
meanB=mean(b);
if meanB<=0 || base.totalAirPowerW<=0
 out.powerStatus='nonpositive_denominator';return;
end
out.meanLoss=(meanB-mean(m))/meanB;
out.powerReduction=(base.totalAirPowerW-method.totalAirPowerW)/base.totalAirPowerW;
if out.meanLoss<.1 && out.powerReduction>=.3
 out.powerStatus='pass';else,out.powerStatus='fail';end
end
