function out=srs_overlap_snapshot(H0,R,pilotId,energy,noise,z,ageMs,tcMs,innovation,portRows,rb,speedKmh,frequencyGhz,priorMean,jointPrior)
%SRS_OVERLAP_SNAPSHOT Finite-bandwidth matched SRS and exact Gaussian posterior.
% Project extension of Arvidson's received-signal model, not fitted Eq.23.
% noise is the post-matched-filter variance from pilotSnrDb; complex CN
% has E|z|^2=1. Gram preserves cross-label matched-noise correlation.
[K,L,Nr,Nt]=size(H0); T=size(z,3); P=size(pilotId,2);
if nargin<14,priorMean=complex(zeros(size(H0)));end
if nargin<15,jointPrior=[];end
assert(isequal(size(priorMean),size(H0)) && ...
    all(isfinite(real(priorMean(:))) & isfinite(imag(priorMean(:)))));
assert(ismember(rb,[16 17]) && T>=max(pilotId(:)) && P==numel(portRows));
assert(T==16,'Two combs with eight cyclic shifts each are required');
N=6*rb; n=(0:N-1)'; shifts=(0:7);
% The common unit-modulus base sequence cancels in the matched Gram.
% Each port occupies one of two disjoint combs; labels 1:8 and 9:16
% encode cyclic shifts within those combs, not sixteen shifts on one comb.
seq=exp(2i*pi*n*shifts/8);
combGram=(seq'*seq)/N;
Gram=blkdiag(combGram,combGram); Gram=(Gram+Gram')/2;
noiseCov=noise*Gram;
noiseFactor=chol(noiseCov+1e-14*max(noise,eps)*eye(T),'lower');
if nargin>=13
    assert(speedKmh>=0 && frequencyGhz>0);
    dopplerHz=(speedKmh/3.6)*(frequencyGhz*1e9)/299792458;
    a=besselj(0,2*pi*dopplerHz*(ageMs/1000));
    ageModel='jakes_isotropic_single_lag_extension';
else
    a=exp(-ageMs/tcMs);
    ageModel='exponential_project_snapshot';
end
out.H_true=a*H0+sqrt(1-a^2)*innovation;
out.H_hat=complex(zeros(size(H0)));
out.C=complex(zeros(Nr,Nr,K,L));
out.Y=complex(zeros(T,Nt,L));
out.sourceRisk=zeros(K,L,2);
out.receivedTargetPower=zeros(K,L);
out.receivedInterferencePower=zeros(K,L);
out.receivedNoisePower=zeros(K,L);
out.receivedPresinr=zeros(K,L);
out.jointPriorKind='independent_ue_block_diagonal';
if ~isempty(jointPrior)
 assert(size(jointPrior,1)==K*Nr && ...
     size(jointPrior,2)==K*Nr && size(jointPrior,3)==L);
 out.jointPriorKind='joint_ue_covariance';
end
out.crossSourceCovariance=complex(zeros(Nr,Nr,K,L));
out.priorCrossFro=zeros(K,K,L);
out.pilotReuseGroup=cell(K,1);
for k=1:K
 out.pilotReuseGroup{k}=find(any(pilotId==pilotId(k,:),2));
end
out.priorMinimumEigen=zeros(1,L);
out.posteriorPsdProjectionFro=zeros(1,L);
out.observationCovariance=complex(zeros(T,T,L));
for l=1:L
 B=complex(zeros(T,K*Nr)); Q=complex(zeros(K*Nr));
 X=complex(zeros(K*Nr,Nt));
 Mean=complex(zeros(K*Nr,Nt));
 for k=1:K
  ix=(k-1)*Nr+(1:Nr);
  X(ix,:)=reshape(H0(k,l,:,:),Nr,Nt);
  Mean(ix,:)=reshape(priorMean(k,l,:,:),Nr,Nt);
  Q(ix,ix)=R(:,:,k,l);
  for p=1:P
   B(pilotId(k,p),ix(portRows(p)))=sqrt(energy(k,p));
  end
 end
 if ~isempty(jointPrior)
  Q=double(jointPrior(:,:,l));Q=(Q+Q')/2;
  out.priorMinimumEigen(l)=min(real(eig(Q)));
  assert(out.priorMinimumEigen(l)>=-1e-8*max(1,norm(Q,2)), ...
      'Joint SRS prior is not positive semidefinite');
  for k=1:K
   ix=(k-1)*Nr+(1:Nr);
   assert(norm(Q(ix,ix)-R(:,:,k,l),'fro')< ...
       1e-6*max(1,norm(R(:,:,k,l),'fro')), ...
       'Joint prior marginal differs from supplied UE prior');
   for j=1:K
    jx=(j-1)*Nr+(1:Nr);
    out.priorCrossFro(k,j,l)=norm(Q(ix,jx),'fro');
   end
  end
 end
 M=Gram*B;
 received=M*X;
 noiseRealization=noiseFactor*reshape(z(:,l,:),Nt,T).';
 Y=received+noiseRealization;
 S=M*Q*M'+noiseCov;S=(S+S')/2;
 out.observationCovariance(:,:,l)=S;
 assert(min(real(eig(S)))>0,'SRS observation covariance is singular');
 posterior=Q*M'/S;
 mu=Mean+posterior*(Y-M*Mean);
 MQ=M*Q;
 out.Y(:,:,l)=Y;
 for k=1:K
  ix=(k-1)*Nr+(1:Nr);
  own=M(:,ix)*X(ix,:);
  other=received-own;
  targetPower=mean(abs(own).^2,'all');
  interferencePower=mean(abs(other).^2,'all');
  noisePower=mean(abs(noiseRealization).^2,'all');
  out.receivedTargetPower(k,l)=targetPower;
  out.receivedInterferencePower(k,l)=interferencePower;
  out.receivedNoisePower(k,l)=noisePower;
  out.receivedPresinr(k,l)=targetPower/ ...
      max(interferencePower+noisePower,realmin);
  out.H_hat(k,l,:,:)=reshape(a*mu(ix,:),1,1,Nr,Nt);
  marginal=Q(ix,ix)-posterior(ix,:)*MQ(:,ix);
  marginal=(marginal+marginal')/2;
  [basis,values]=eig(marginal);
  eigenvalues=real(diag(values));
  assert(min(eigenvalues)>=-1e-7*max(1,norm(Q(ix,ix),2)), ...
      'Joint SRS marginal posterior is materially non-PSD');
  projected=basis*diag(max(0,eigenvalues))*basis';
  out.posteriorPsdProjectionFro(l)=max( ...
      out.posteriorPsdProjectionFro(l),norm(projected-marginal,'fro'));
  out.C(:,:,k,l)=a^2*(projected+projected')/2+ ...
      (1-a^2)*R(:,:,k,l);
  if ~isempty(jointPrior)
   other=setdiff(1:K*Nr,ix);
   J=M(:,ix);Mothers=M(:,other);
   A=posterior(ix,:);F=eye(Nr)-A*J;
   cross=-F*Q(ix,other)*Mothers'*A'- ...
       A*Mothers*Q(other,ix)*F';
   out.crossSourceCovariance(:,:,k,l)= ...
       a^2*(cross+cross')/2;
  end
  peers=find(any(pilotId==pilotId(k,:),2)); peers(peers==k)=[];
  out.sourceRisk(k,l,1)=sum(arrayfun(@(j)real(trace(R(:,:,j,l))),peers));
  out.sourceRisk(k,l,2)=ageMs;
 end
end
out.a=a;
if nargin>=13, out.modelId='srs_overlap_jakes_v3';
else, out.modelId='srs_overlap_snapshot_v1'; end
out.ageModel=ageModel;
out.noiseCov=noiseCov; out.gram=Gram;
end
