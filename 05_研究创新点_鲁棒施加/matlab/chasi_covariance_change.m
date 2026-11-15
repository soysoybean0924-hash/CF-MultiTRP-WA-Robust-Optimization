function value=chasi_covariance_change(old,new)
% Frobenius distance of per-TRP transmit covariance, matching Python campaign.
num=0;den=0;
for g=1:numel(old)
 x=old{g};y=new{g};Nt=size(x,3);
 for l=1:size(x,1)
  a=reshape(permute(x(l,:,:,:),[3,2,4,1]),Nt,[]);
  b=reshape(permute(y(l,:,:,:),[3,2,4,1]),Nt,[]);
  q0=a*a';q1=b*b';num=num+sum(abs(q1-q0).^2,'all');
  den=den+sum(abs(q0).^2,'all');
 end
end
value=sqrt(num)/max(1,sqrt(den));
end
