function [resources,record]=srs_experimental_service(resources,data,serviceTrps)
%SRS_EXPERIMENTAL_SERVICE Measured-statistic-only service ablation.
% Pilot Gram, sender powers and prior gains form a dimensionless interference
% exposure. It is a decision feature, not a fitted source covariance.
arguments
 resources (1,:) struct
 data (1,1) struct
 serviceTrps (1,1) double {mustBeInteger,mustBePositive}=2
end
pilot=double(data.pilotId); energy=double(data.pilotEnergy);
prior=double(data.priorRx); senders=size(pilot,1);
assert(isequal(size(energy),size(pilot)) && size(prior,3)==senders);
record=struct('model','experimental_unvalidated_source_risk', ...
 'changedLinks',0,'totalLinks',0,'score',{cell(1,numel(resources))}, ...
 'interferenceToOwnPrior',{cell(1,numel(resources))}, ...
 'ageMs',zeros(1,numel(resources)));
for g=1:numel(resources)
 r=resources(g); [K,L,~,~]=size(r.H);
 assert(isfield(r,'srs_diagnostic') && serviceTrps<=L);
 gram=double(r.srs_diagnostic.pilotGram);
 assert(size(gram,1)>=max(pilot(:)));
 old=logical(r.servingMask); mask=false(L,K); score=zeros(K,L);
 exposure=zeros(K,L);
 for k=1:K
  id=double(r.ueIds(k));
  for l=1:L
   h=squeeze(double(r.H(k,l,:,:)));
   c=squeeze(double(r.C(k,l,:,:)));
   est=sum(abs(h).^2,'all');
   uncertainty=real(trace(c))*size(h,2);
   noise=double(r.srs_noise_variance);
   % Matched pilot leakage from every SRS sender, including unscheduled UEs.
   own=0; interferers=0;
   for p=1:size(pilot,2)
    own=own+energy(id,p)*real(trace(prior(:,:,id,l)));
    for j=1:senders
     if j==id, continue; end
     coupling=abs(gram(pilot(id,p),pilot(j,p)))^2;
     interferers=interferers+energy(j,p)*coupling* ...
         real(trace(prior(:,:,j,l)));
    end
   end
   exposure(k,l)=interferers/max(own,eps);
   reliability=est/(uncertainty+noise+eps);
   score(k,l)=reliability/(1+exposure(k,l));
  end
  [~,order]=sort(score(k,:),'descend');
  mask(order(1:serviceTrps),k)=true;
 end
 resources(g).servingMask=mask;
 record.changedLinks=record.changedLinks+nnz(old~=mask);
 record.totalLinks=record.totalLinks+numel(mask);
 record.score{g}=score;
 record.interferenceToOwnPrior{g}=exposure;
 record.ageMs(g)=double(r.age_ms);
end
end
