function out=chasi_v3_link_power(resources,precoder,budget,risk,opts)
% Per UE-TRP-group amplitude allocation on fixed WA beam directions.
% Maximizes a declared separable estimated-link utility, not true throughput.
% p_lkg in W; sum_g nu_g sum_k p_lkg <= budget_l. No truth enters decisions.
G=numel(resources);L=size(precoder{1},1);assert(numel(budget)==L);
assert(numel(risk)==G);budget=double(budget(:));
assert(all(isfinite(budget)&budget>0));
if ~isfield(opts,'powerPrice'),opts.powerPrice=0;end
if ~isfield(opts,'iterations'),opts.iterations=30;end
if ~isfield(opts,'step'),opts.step=.1;end
assert(opts.iterations>=1 && opts.step>0 && opts.powerPrice>=0);
Kmax=max(arrayfun(@(r)size(r.H,1),resources));
p=zeros(L,Kmax,G);p0=p;gain=p;limit=p;weight=zeros(G,1);
for g=1:G
 r=resources(g);K=size(r.H,1);x=precoder{g};
 assert(isequal(size(risk{g}),[K,L]));
 assert(all(isfinite(risk{g}(:)) & risk{g}(:)>=0 & ...
     risk{g}(:)<=1));
 if isfield(r,'resourceWeightNu'),weight(g)=double(r.resourceWeightNu);
 else,weight(g)=double(r.rb_count)/273;end
 for k=1:K,for l=1:L
  v=reshape(x(l,k,:,:),size(x,3),size(x,4));
  p0(l,k,g)=sum(abs(v).^2,'all');
  limit(l,k,g)=p0(l,k,g)*(1-risk{g}(k,l));
  h=reshape(double(r.H(k,l,:,:)),size(r.H,3),size(r.H,4));
  gain(l,k,g)=sum(abs(h*v).^2,'all')/max(p0(l,k,g),eps);
 end,end
end
p=limit;
for it=1:opts.iterations
 for g=1:G
  K=size(resources(g).H,1);
  for k=1:K,for l=1:L
   if p0(l,k,g)==0,p(l,k,g)=0;continue;end
   gr=gain(l,k,g)*max(0,1-risk{g}(k,l));
   gradient=gr/(1+gr*p(l,k,g))-opts.powerPrice;
   p(l,k,g)=min(limit(l,k,g), ...
       max(0,p(l,k,g)+opts.step*gradient));
  end,end
 end
 for l=1:L
  total=sum(reshape(p(l,:,:),Kmax,G).*weight','all');
  if total>budget(l)
   p(l,:,:)=p(l,:,:)*(budget(l)/total);
  end
 end
end
vout=precoder;
for g=1:G
 K=size(resources(g).H,1);
 for k=1:K,for l=1:L
  if p0(l,k,g)>0
   vout{g}(l,k,:,:)=vout{g}(l,k,:,:)*sqrt(p(l,k,g)/p0(l,k,g));
  end
 end,end
end
actual=zeros(L,1);
for l=1:L,for g=1:G
 actual(l)=actual(l)+weight(g)*sum(abs(vout{g}(l,:,:,:)).^2,'all');
end,end
assert(all(actual<=budget+1e-6),'TRP power infeasible');
out=struct('precoder',{vout},'linkPowerW',p,'initialLinkPowerW',p0, ...
 'riskLimitedLinkPowerW',limit, ...
 'actualTrpPowerW',actual,'budgetW',budget, ...
 'status','feasible_estimated_link_utility', ...
 'limitation','fixed WA beam directions and separable estimated-link proxy');
end
