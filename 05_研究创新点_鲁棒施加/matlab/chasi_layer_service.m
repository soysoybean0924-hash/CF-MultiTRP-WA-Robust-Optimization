function mask=chasi_layer_service(r,level,mode,margin,loadPenalty,weights)
% At most one TRP replacement per UE, following Python service_policy.py.
mask=logical(r.servingMask);[L,K]=size(mask);assert(all(sum(mask,1)==2));
gain=squeeze(sum(sum(abs(r.H_hat).^2,4),3));
total=zeros(K,L);source=zeros(K,L,4);
for k=1:K,for l=1:L
 total(k,l)=real(trace(squeeze(r.C(k,l,:,:))));
 for q=1:4,source(k,l,q)=real(trace(squeeze(r.sourceCov(k,l,:,:,q))));end
end,end
assert(all(source(:)>=-1e-7));source=max(source,0);
if strcmp(mode,'uniform'),numer=total;
else,numer=sum(source.*reshape(weights,1,1,4),3);end
denom=max(gain+total,1e-30);
risk=numer./denom;signal=gain./denom;
load=sum(mask,2)/K;
score=signal-level*risk-loadPenalty*load';
old=mask;
for k=1:K
 present=find(old(:,k));absent=find(~old(:,k));
 [~,oi]=min(score(k,present));[~,ii]=max(score(k,absent));
 outgoing=present(oi);incoming=absent(ii);
 if score(k,incoming)-score(k,outgoing)>margin
  mask(outgoing,k)=false;mask(incoming,k)=true;
 end
end
end
