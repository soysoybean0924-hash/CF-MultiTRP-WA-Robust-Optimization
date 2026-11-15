function value=chasi_v6_tr38901_phase_cross_kernel(rho)
% Cross moment of a uniform phase obtained from a correlated normal field.
% f(z)=exp(2*pi*i*Phi(z)). The nonnegative Hermite spectrum makes the
% finite-order covariance positive semidefinite for any Gaussian field.
% Quadrature/truncation error is deliberately retained behind validation.
assert(all(isfinite(rho(:))) && all(abs(rho(:))<=1+1e-10));
persistent spectrum
if isempty(spectrum)
 nodes=48;order=24;
 off=sqrt(1:nodes-1);
 jacobi=diag(off,1)+diag(off,-1);
 [vectors,diagonal]=eig(jacobi);
 z=real(diag(diagonal));weights=real(vectors(1,:)).^2;
 f=exp(2i*pi*.5*(1+erf(z/sqrt(2))));
 spectrum=zeros(order,1);
 previous=ones(nodes,1);current=z;
 for n=1:order
  if n==1,hermite=current;
  else
   next=(z.*current-sqrt(n-1)*previous)/sqrt(n);
   previous=current;current=next;hermite=current;
  end
  coefficient=sum(weights(:).*f(:).*hermite);
  spectrum(n)=abs(coefficient)^2;
 end
 spectrum=spectrum/sum(spectrum);
end
rho=max(-1,min(1,double(rho)));
value=zeros(size(rho));power=rho;
for n=1:numel(spectrum)
 value=value+spectrum(n)*power;
 power=power.*rho;
end
end
