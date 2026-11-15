function [perUe,perStreamSinrDb]=chasi_v3_true_rates(resources,precoder,totalUsers)
% True-channel static SE after every final amplitude update; offline only.
perUe=zeros(totalUsers,1);
streams=size(precoder{1},4);perStreamSinrDb=nan(totalUsers,streams);
for g=1:numel(resources)
 r=resources(g);H=double(r.H_true);v=precoder{g};
 [K,L,Nr,Nt]=size(H);
 if isfield(r,'resourceWeightNu'),nu=double(r.resourceWeightNu);
 else,nu=double(r.rb_count)/273;end
 for k=1:K
  allCov=eye(Nr);interf=eye(Nr);
  desired=complex(zeros(Nr,streams));
  for j=1:K
   eff=complex(zeros(Nr,streams));
   for l=1:L
    eff=eff+reshape(H(k,l,:,:),Nr,Nt)*reshape(v(l,j,:,:),Nt,streams);
   end
   allCov=allCov+eff*eff';
   if j~=k,interf=interf+eff*eff';end
   if j==k,desired=eff;end
  end
  rate=max(0,(real(log(det(allCov)))-real(log(det(interf))))/log(2));
  perUe(double(r.ueIds(k)))=perUe(double(r.ueIds(k)))+nu*rate;
  if nargout>1
   for s=1:streams
    h=desired(:,s);R=allCov-h*h';R=(R+R')/2;
    if ~any(h),continue;end
    sinr=max(0,real(h'*(R\h)));
    idx=double(r.ueIds(k));value=10*log10(max(sinr,realmin));
    if isnan(perStreamSinrDb(idx,s)),perStreamSinrDb(idx,s)=value;
    else,perStreamSinrDb(idx,s)=min(perStreamSinrDb(idx,s),value);end
   end
  end
 end
end
end
