function [ray,rx,tx,identity]=chasi_v6_tr38901_rebuild_link( ...
    link,siteField,drop,meta,trpXY,tti,timeSecondsOverride)
% Frozen 2-D fields evaluated at a straight-path UE position at absolute TTI.
% TR 38.901 7.6.3.1 / 7.6.3.4 field contract with Procedure B candidate
% Steps 5-7. Full spatial-consistency status remains an explicit gate.
assert(tti>=0 && tti==fix(tti) && isfinite(tti));
elapsed=tti*double(meta.ttiMs)/1000;
if nargin>=7,elapsed=timeSecondsOverride;end
assert(isfinite(elapsed) && elapsed>=0);
velocity=double(link.ray.velocityMps(:).');
xy=double(link.ray.initialXY(:).')+elapsed*velocity(1:2);
trp=double(trpXY(:).');
assert(numel(xy)==2 && numel(trp)==2 && all(isfinite(xy)));
isLos=chasi_v6_tr38901_site_field('los',siteField,xy,trp);
draw=chasi_v6_tr38901_site_field('evaluate',siteField,xy,isLos);
frozen=drop{double(isLos)+1}.ray;
draw.ssp.fixedStrongestOriginalIndices= ...
    drop{double(isLos)+1}.ssp.fixedStrongestOriginalIndices;
% These variables are fixed for one drop by 7.6.3.1. The other ray and
% cluster random quantities are evaluated from the same spatial fields.
draw.ray.centerSign=frozen.centerSign;
draw.ray.couplingScore=frozen.couplingScore;
draw.ray.clusterSign=frozen.clusterSign;
scenario=char(meta.scenario);fc=double(meta.carrierHz)/1e9;
if strcmp(scenario,'UMi'),hbs=10;else,hbs=25;end
d2=max(10,norm(xy-trp));d3=hypot(d2,hbs-1.5);
p=chasi_v6_tr38901_lsp_table(scenario,logical(isLos),fc,d2);
stream=RandStream('mt19937ar','Seed',double(meta.generationSeed)+ ...
    1000000*double(link.ueId)+double(siteField.site));
lsp=chasi_v6_tr38901_lsp(p,stream,draw.lspNormal);
cluster=chasi_v6_tr38901_clusters(p,lsp,stream,draw.ssp);
[ray,rx,tx]=chasi_v6_tr38901_rays(p,lsp,cluster,stream, ...
    trp,xy,velocity,link.bearingRad,meta.carrierHz,draw.ray);
pl=chasi_v6_tr38901_outdoor_pathloss( ...
    scenario,isLos,fc,d2,d3,hbs)-lsp.sfDb;
ray.power=ray.power*10^(-pl/10)/double(meta.channelNormalizationNoiseW);
ray.initialXY=xy;
transitionStart=0;
transitionCount=0;
if elapsed>0 && nargin<7
 % Diffuse Doppler is relative to the frozen drop. The LOS term already
 % contains the current geometric carrier phase from Step 10.
 diffuse=ray.clusterIndex>0;
 % 7.6.3.3 permits a hard LOS/NLOS transition. Restart the phase integral
 % for the newly selected condition; do not splice rays from two states.
 [transitionStart,transitionCount]=chasi_v6_los_transition_events( ...
     siteField,double(link.ray.initialXY),velocity,trp,elapsed);
 assert(chasi_v6_tr38901_site_field('los',siteField, ...
     double(link.ray.initialXY(:).')+elapsed*velocity(1:2),trp)==isLos);
 duration=elapsed-transitionStart;
 steps=max(1,ceil(norm(velocity)*duration/.25));
 dt=duration/steps;phaseIntegral=zeros(nnz(diffuse),1);
 for step=1:steps
  midpoint=transitionStart+(step-.5)*dt;
  [midRay,~,~,midIdentity]=chasi_v6_tr38901_rebuild_link( ...
      link,siteField,drop,meta,trpXY,tti,midpoint);
  assert(midIdentity.isLos==isLos && ...
      nnz(midRay.clusterIndex>0)==nnz(diffuse), ...
      'Unresolved LOS transition within phase integration segment');
  phaseIntegral=phaseIntegral+ ...
      (midRay.arrivalUnit(midRay.clusterIndex>0,:)*velocity.')*dt;
 end
 ray.phasesRad(diffuse,:)=ray.phasesRad(diffuse,:)+ ...
     2*pi*phaseIntegral*meta.carrierHz/299792458;
end
identity=struct('absoluteTti',tti,'positionMeters',xy, ...
    'isLos',isLos,'pathlossDb',pl,'site',siteField.site, ...
    'conditionEpochSeconds',transitionStart, ...
    'conditionTransitionCount',transitionCount);
end
