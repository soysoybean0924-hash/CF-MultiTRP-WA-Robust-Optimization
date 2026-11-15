function out=chasi_v6_estimated_experience_proxy(resources,beams,users,bandwidthHz)
% Causal Gaussian log-det rate proxy [Mbps] from H_hat and physical C.
% No truth, future ACK, MCS or policyCov is accepted by this evaluator.
assert(numel(resources)==numel(beams) && users>=1 && bandwidthHz>0);
perUe=zeros(users,1);scheduled=false(users,1);
for g=1:numel(resources)
 r=resources(g);
 assert(~isfield(r,'H_true') && isfield(r,'H') && isfield(r,'C'), ...
     'Proxy must receive decision-time H_hat and physical C only');
 H=double(r.H);C=double(r.C);v=double(beams{g});
 [K,L,Nr,Nt]=size(H);streams=size(v,4);
 assert(isequal([size(v,1),size(v,2),size(v,3),size(v,4)], ...
     [L,K,Nt,streams]) && ...
     isequal([size(C,1),size(C,2),size(C,3),size(C,4)], ...
     [K,L,Nr,Nr]));
 ids=double(r.ueIds(:));assert(numel(ids)==K);
 nu=double(r.resourceWeightNu);
 assert(isscalar(nu) && isfinite(nu) && nu>=0);
 trpPower=zeros(L,1);
 for l=1:L,trpPower(l)=sum(abs(v(l,:,:,:)).^2,'all');end
 for k=1:K
  desired=complex(zeros(Nr,streams));interference=eye(Nr);
  for i=1:K
   effective=complex(zeros(Nr,streams));
   for l=1:L
    effective=effective+reshape(H(k,l,:,:),Nr,Nt)* ...
        reshape(v(l,i,:,:),Nt,streams);
   end
   if i==k,desired=effective;
   else,interference=interference+effective*effective';end
  end
  for l=1:L
   interference=interference+trpPower(l)*reshape(C(k,l,:,:),Nr,Nr);
  end
  interference=(interference+interference')/2;
  allCov=interference+desired*desired';
  [li,flagI]=chol(interference,'lower');
  [la,flagA]=chol(allCov,'lower');
  assert(flagI==0 && flagA==0, ...
      'Estimated proxy covariance is not positive definite');
  spectralEfficiency=max(0,2*(sum(log(real(diag(la))))- ...
      sum(log(real(diag(li)))))/log(2));
  perUe(ids(k))=perUe(ids(k))+nu*spectralEfficiency* ...
      bandwidthHz/1e6;
  scheduled(ids(k))=true;
 end
end
assert(all(isfinite(perUe)) && all(perUe>=0));
out=struct('perUeMbps',perUe,'scheduledMask',scheduled, ...
    'meanScheduledMbps',mean(perUe(scheduled)), ...
    'modelId','causal_Hhat_physicalC_Gaussian_logdet_v1', ...
    'unit','Mbps per scheduled UE over current resource-group weights', ...
    'informationSet','current H_hat, physical C, geometry, schedule and beams');
end
