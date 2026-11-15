function out=srs_arvidson_forward(H0,pilotId,coeff,noiseDraw)
%SRS_ARVIDSON_FORWARD Forward simulation of thesis Eq. 23, per coefficient.
% H0: [all SRS senders,TRP,RX,TX]. No scheduled-user filtering is allowed.
% coeff.betaSamePilot and coeff.alpha/gamma require independent calibration.
% The result is a forward perturbation, not a conditional Gaussian posterior.
% Eq.23's reported scalar standard deviation needs an explicit real/complex
% convention before these unit-complex-variance draws can be calibrated.
% Do not combine this output with srs_overlap_snapshot or srs_joint_snapshot.
[K,L,Nr,Nt]=size(H0);
assert(isequal(size(noiseDraw),size(H0)));
assert(size(pilotId,1)==K && isfield(coeff,'betaSamePilot'));
assert(isscalar(coeff.betaSamePilot) && coeff.betaSamePilot>=0);
out.H_hat=complex(zeros(size(H0)));
out.sigma=zeros(K,L,Nr,Nt);
out.pilot=zeros(K,L,Nr,Nt);
out.noise=zeros(K,L,Nr,Nt);
out.residual=zeros(K,L,Nr,Nt);
for k=1:K
 colliders=find(any(pilotId==pilotId(k,:),2));
 colliders(colliders==k)=[];
 for l=1:L
  for nr=1:Nr
   for nt=1:Nt
    target=abs(H0(k,l,nr,nt));
    interfering=reshape(abs(H0(colliders,l,nr,nt)),1,[]);
    c=struct('alpha',coeff.alpha,'beta',...
      coeff.betaSamePilot*ones(numel(colliders),1),'gamma',coeff.gamma);
    [sig,p]=arvidson_srs_sigma(target,interfering,c);
    out.sigma(k,l,nr,nt)=sig;
    out.pilot(k,l,nr,nt)=p.pilot;
    out.noise(k,l,nr,nt)=p.background;
    out.residual(k,l,nr,nt)=p.truncation;
   end
  end
 end
end
% noiseDraw must be unit-variance complex Gaussian, generated independently.
out.H_hat=H0+out.sigma.*noiseDraw;
out.modelId='arvidson_eq23_forward_requires_calibration';
end
