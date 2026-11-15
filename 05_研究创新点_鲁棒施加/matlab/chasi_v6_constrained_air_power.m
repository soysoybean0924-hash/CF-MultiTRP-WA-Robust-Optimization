function out=chasi_v6_constrained_air_power(resources,sourceBeam,independentBeam,users,bandwidthHz,opts,calibration)
% Fixed-direction link-power search: minimize actual weighted air power
% within tested uniform/coordinate reductions subject to a causal proxy.
% The nonconvex problem has no global-optimum claim. Future ACK is absent.
assert(numel(resources)==numel(sourceBeam) && ...
    numel(resources)==numel(independentBeam));
assert(isfield(opts,'maximumMeanLoss') && ...
    opts.maximumMeanLoss>0 && opts.maximumMeanLoss<1);
base=chasi_v6_estimated_experience_proxy( ...
    resources,independentBeam,users,bandwidthHz);
initial=chasi_v6_estimated_experience_proxy( ...
    resources,sourceBeam,users,bandwidthHz);
assert(isequal(base.scheduledMask,initial.scheduledMask));
margin=0;calibrated=false;calibrationSha='';
if ~isempty(calibration)
 assert(isfield(calibration,'validationPassed') && ...
     isfield(calibration,'optimismMargin') && ...
     isfield(calibration,'proxyModelId') && ...
     strcmp(calibration.proxyModelId,base.modelId));
 calibrated=logical(calibration.validationPassed);
 if calibrated,margin=double(calibration.optimismMargin);end
 if isfield(calibration,'fileSha256')
  calibrationSha=calibration.fileSha256;
 end
end
assert(isfinite(margin) && margin>=0 && margin<opts.maximumMeanLoss);
threshold=opts.maximumMeanLoss-margin-1e-6;
assert(threshold>0);
budget=double(opts.powerBudgetW(:));
initialPower=airPower(resources,sourceBeam);
assert(numel(budget)==numel(initialPower) && ...
    all(initialPower<=budget+1e-6));
beam=sourceBeam;power=initialPower;
trajectory=struct('iteration',{},'phase',{},'group',{},'trp',{}, ...
    'ue',{},'fraction',{},'meanProxyMbps',{},'proxyLoss',{}, ...
    'totalAirPowerW',{},'trpPowerW',{},'constraintFeasible',{});
reference=base.meanScheduledMbps;
if ~isfinite(reference) || reference<=0
 reason='independent_proxy_denominator_nonpositive';feasible=false;
elseif ~all(isfinite(initial.perUeMbps))
 reason='nonfinite_initial_proxy';feasible=false;
else
 feasible=loss(initial.meanScheduledMbps,reference)<=threshold;
 if feasible,reason='search_started';
 else,reason='source_beam_proxy_infeasible';end
end
evaluationCount=0;
if feasible
 % A bounded global amplitude ray gives a reproducible feasible start.
 grid=linspace(0,1,21);selected=1;lower=0;
 for j=1:numel(grid)
  q=grid(j);trial=scaleBeams(sourceBeam,sqrt(q));
  metric=chasi_v6_estimated_experience_proxy( ...
      resources,trial,users,bandwidthHz);
  evaluationCount=evaluationCount+1;
  trialPower=airPower(resources,trial);
  okay=loss(metric.meanScheduledMbps,reference)<=threshold;
  trajectory(end+1)=row(numel(trajectory)+1,'uniform_grid',0,0,0, ...
      q,metric.meanScheduledMbps,reference,trialPower,okay); %#ok<AGROW>
  if okay,selected=q;beam=trial;power=trialPower;break;end
  lower=q;
 end
 if selected>lower && selected>0
  for j=1:12
   q=(lower+selected)/2;trial=scaleBeams(sourceBeam,sqrt(q));
   metric=chasi_v6_estimated_experience_proxy( ...
       resources,trial,users,bandwidthHz);
   evaluationCount=evaluationCount+1;
   trialPower=airPower(resources,trial);
   okay=loss(metric.meanScheduledMbps,reference)<=threshold;
   trajectory(end+1)=row(numel(trajectory)+1,'uniform_refine', ...
       0,0,0,q,metric.meanScheduledMbps,reference,trialPower,okay); %#ok<AGROW>
   if okay,selected=q;beam=trial;power=trialPower;
   else,lower=q;end
  end
 end
 % Coordinate refinement preserves the same constraint at every accepted
 % point. Each coordinate is a UE-TRP-frequency-group link power.
 links=[];
 for g=1:numel(resources)
  nu=double(resources(g).resourceWeightNu);
  for l=1:size(beam{g},1),for k=1:size(beam{g},2)
   p=sum(abs(beam{g}(l,k,:,:)).^2,'all');
   if p>0,links(end+1,:)=[nu*p,g,l,k];end %#ok<AGROW>
  end,end
 end
 if ~isempty(links)
  links=sortrows(links,-1);
  links=links(1:min(size(links,1),32),:);
  for j=1:size(links,1)
   g=links(j,2);l=links(j,3);k=links(j,4);
   current=beam{g}(l,k,:,:);
   if sum(abs(current).^2,'all')<=0,continue;end
   for fraction=[0 .25 .5 .75]
    trial=beam;trial{g}(l,k,:,:)=sqrt(fraction)*current;
    metric=chasi_v6_estimated_experience_proxy( ...
        resources,trial,users,bandwidthHz);
    evaluationCount=evaluationCount+1;
    trialPower=airPower(resources,trial);
    okay=loss(metric.meanScheduledMbps,reference)<=threshold;
    trajectory(end+1)=row(numel(trajectory)+1, ...
        'link_coordinate',g,l,k,fraction, ...
        metric.meanScheduledMbps,reference,trialPower,okay); %#ok<AGROW>
    if okay,beam=trial;power=trialPower;break;end
   end
  end
 end
 reason='bounded_feasible_search_complete';
end
final=chasi_v6_estimated_experience_proxy( ...
    resources,beam,users,bandwidthHz);
finalLoss=loss(final.meanScheduledMbps,reference);
feasible=feasible && isfinite(finalLoss) && finalLoss<=threshold && ...
    all(power<=budget+1e-6);
if ~feasible
 beam=sourceBeam;power=initialPower;
 final=initial;finalLoss=loss(final.meanScheduledMbps,reference);
end
out=struct('precoder',{beam},'trajectory',trajectory, ...
    'unconstrainedSourcePrecoder',{sourceBeam}, ...
    'status',reason,'feasible',feasible, ...
    'proxyModelId',base.modelId,'proxyCalibrationSha256',calibrationSha, ...
    'calibrationValidated',calibrated,'optimismMargin',margin, ...
    'maximumMeanLoss',opts.maximumMeanLoss, ...
    'proxyLossLimitUsed',threshold, ...
    'independentProxyMeanMbps',reference, ...
    'sourceInitialProxyMeanMbps',initial.meanScheduledMbps, ...
    'finalProxyMeanMbps',final.meanScheduledMbps, ...
    'finalProxyLoss',finalLoss, ...
    'initialTrpPowerW',initialPower,'finalTrpPowerW',power, ...
    'initialTotalAirPowerW',sum(initialPower), ...
    'finalTotalAirPowerW',sum(power), ...
    'evaluationCount',evaluationCount, ...
    'variable','per group x TRP x UE link power, fixed WA beam direction', ...
    'domain','0 <= p_link <= source link p; per TRP weighted sum <= budget', ...
    'optimalityScope','bounded ray plus one coordinate pass; no global certificate');
end

function p=airPower(resources,beam)
L=size(beam{1},1);p=zeros(L,1);
for g=1:numel(resources)
 nu=double(resources(g).resourceWeightNu);
 for l=1:L,p(l)=p(l)+nu*sum(abs(beam{g}(l,:,:,:)).^2,'all');end
end
end

function value=loss(meanValue,reference)
value=(reference-meanValue)/reference;
end

function v=scaleBeams(original,scale)
v=original;for g=1:numel(v),v{g}=scale*v{g};end
end

function item=row(iteration,phase,g,l,k,fraction,meanValue,reference,power,okay)
item=struct('iteration',iteration,'phase',phase,'group',g,'trp',l, ...
    'ue',k,'fraction',fraction,'meanProxyMbps',meanValue, ...
    'proxyLoss',loss(meanValue,reference), ...
    'totalAirPowerW',sum(power),'trpPowerW',power, ...
    'constraintFeasible',okay);
end
