function [resources,meta]=chasi_layer_protection(resources,layer,cfg,mode,level,ablation,trigger)
% Matches layered.py: physical C_stat is retained; design C+policyCov enters U/W and V.
if nargin<6 || isempty(ablation),ablation='none';end
if nargin<7 || isempty(trigger),trigger='proposed';end
assert(ismember(mode,{'uniform','variable','fixed','mapped'}));
assert(numel(resources)>=1 && numel(resources)<=17);
w=cfg.mapping; weight=[w.residual_weight w.pilot_weight w.noise_weight w.age_weight];
% Python protection computes mapped risk over the selected resource groups.
selected=ismember([layer.records.group],1:numel(resources));
shares=vertcat(layer.records(selected).source_share);
risk=mean(shares*weight(:));
mapped=min(w.ceiling,max(w.floor,w.floor+risk));
factor=sqrt(layer.radius_squared);
if strcmp(mode,'mapped'),base=mapped*level/max(cfg.fixed_level,1e-12)*factor;
else,base=level*factor;end
if strcmp(mode,'uniform'),f=[1 1 1];
else,f=[cfg.variable_factors.beam cfg.variable_factors.service cfg.variable_factors.power];end
strength=min(.3,max(0,base*f));
switch ablation
 case 'beam',strength(1)=0; case 'service',strength(2)=0; case 'power',strength(3)=0;
 case 'none', otherwise,error('invalid ablation');
end
cap=resources(1).powerWattPerTrp;
effective=cap*(1-cfg.power_reserve_fraction_max*strength(3)/.3);
changes=0;
for j=1:numel(resources),resources(j).policyCov=[];end
for g=1:numel(resources)
 r=resources(g);r.policyCov=strength(1)*r.C;
 if ~strcmp(trigger,'off') && strength(2)>0
  old=r.servingMask;
  margin=cfg.service.trigger_margin;
  if strcmp(trigger,'fixed'),margin=cfg.service.fixed_margin;end
  if strcmp(mode,'uniform'),rule='uniform';else,rule='source';end
  r.servingMask=chasi_layer_service(r,strength(2),rule,margin, ...
   cfg.service.load_penalty,cfg.service_source_weights);
  changes=changes+nnz(r.servingMask~=old);
 end
 r.powerWattPerTrp=effective;resources(g)=r;
end
meta=struct('mode',mode,'level',level,'source_risk',risk, ...
 'mapped_level',mapped,'calibration_factor',factor,'strengths',strength, ...
 'physical_cap_w',cap,'effective_cap_w',effective,'changed_service_links',changes);
end


