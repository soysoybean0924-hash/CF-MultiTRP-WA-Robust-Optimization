function out=chasi_independent_power_control(resources,cfg)
%CHASI_INDEPENDENT_POWER_CONTROL Local-link waterfill, no TRP cooperation.
% Each TRP uses only its home-link estimated H. The local surrogate ignores
% inter-TRP interference; true network rates are evaluated afterward.
assert(isfield(cfg,'rank') && isfield(cfg,'budget') && isfield(cfg,'users'));
G=numel(resources); L=size(resources(1).H,2); rank=cfg.rank;
local=resources; beams=cell(1,G); gains=cell(1,G); precoder=cell(1,G);
for g=1:G
 r=resources(g); K=size(r.H,1); Nt=size(r.H,4);
 local(g).servingMask=false(L,K);
 beams{g}=complex(zeros(L,K,Nt,rank)); gains{g}=zeros(K,1);
 for k=1:K
  l=double(r.home_trp(k)); assert(l>=1 && l<=L);
  h=reshape(double(r.H(k,l,:,:)),size(r.H,3),Nt);
  [~,~,v]=svd(h,'econ'); v=v(:,1:rank);
  beams{g}(l,k,:,:)=reshape(v,1,1,Nt,rank);
  local(g).servingMask(l,k)=true;
  gains{g}(k)=sum(abs(h*v).^2,'all')/rank;
 end
end
for g=1:G, precoder{g}=complex(zeros(size(beams{g}))); end
for l=1:L
 allGain=[]; allWeight=[]; allGroup=[]; allUser=[];
 for g=1:G
  users=find(local(g).servingMask(l,:)); users=users(:);
  allGain=[allGain;gains{g}(users)]; %#ok<AGROW>
  allWeight=[allWeight;double(resources(g).resourceWeightNu)*ones(numel(users),1)]; %#ok<AGROW>
  allGroup=[allGroup;g*ones(numel(users),1)]; %#ok<AGROW>
  allUser=[allUser;users(:)]; %#ok<AGROW>
 end
 if isempty(allGain) || max(allGain)==0, continue; end
 low=eps; high=max(allGain);
 for t=1:80
  mid=(low+high)/2;
  p=max(0,1/mid-1./max(allGain,eps));
  if sum(allWeight.*p)>cfg.budget, low=mid;
  else, high=mid; end
 end
 p=max(0,1/high-1./max(allGain,eps));
 for i=1:numel(p)
  g=allGroup(i); k=allUser(i);
  precoder{g}(l,k,:,:)=sqrt(p(i)/rank)*beams{g}(l,k,:,:);
 end
end
stateCfg=cfg; stateCfg.branch='matrix_pf_nominal';
stateCfg.outer_iterations=0;
if ~isfield(stateCfg,'pf_floor'), stateCfg.pf_floor=.01; end
out=chasi_matrix_wmmse(local,stateCfg,precoder);
out.branch='independent_pc_local_waterfill';
out.reason='local_sumlog_waterfill_network_evaluated';
out.localPowerControl='per_TRP_weighted_budget_local_gain_only';
end
