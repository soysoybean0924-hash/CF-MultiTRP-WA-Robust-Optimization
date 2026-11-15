function [estimate,C,source,diagnostic,crossOut]=chasi_v6_temporal_posterior( ...
    estimate,C,source,previous,current,crossIn)
% Linear conditional projection from receiver-side model moments.
% For X_t-m_t=B(X_s-m_s)+innovation, B=Cov(t,s)/Cov(s,s).
% Five generating terms transform by B; only aging gains innovation.
[K,L,Nr,Nt]=size(estimate);
if nargin<6 || isempty(crossIn)
    crossIn=complex(zeros(K,L,Nr,Nr));
end
assert(isequal(size(crossIn),[K,L,Nr,Nr]));
crossOut=crossIn;
assert(isequal(size(source),[K,L,Nr,Nr,5]) && ...
    isequal(size(C),[K,L,Nr,Nr]) && ...
    current.absoluteTti>previous.absoluteTti && ...
    isequal(current.ueIds,previous.ueIds));
worstProjection=0;worstReconstruction=0;worstPriorJitter=0;
for k=1:K,for l=1:L
 oldQ=previous.priorRx(:,:,k,l);
 newQ=current.priorRx(:,:,k,l);
 cross=current.crossRx(:,:,k,l);
 oldQ=(oldQ+oldQ')/2;
 minimum=min(real(eig(oldQ)));
 floorValue=1e-10*max(real(trace(oldQ))/Nr,eps);
 jitter=max(0,floorValue-minimum);
 worstPriorJitter=max(worstPriorJitter,jitter);
 B=cross/(oldQ+jitter*eye(Nr));
 innovation=newQ-B*oldQ*B';
 innovation=(innovation+innovation')/2;
 [basis,eigenvalue]=eig(innovation);
 clipped=max(0,real(diag(eigenvalue)));
 projected=basis*diag(clipped)*basis';
 worstProjection=max(worstProjection, ...
     norm(projected-innovation,'fro')/max(norm(newQ,'fro'),eps));
 previousMean=reshape(double(previous.meanH(k,l,:,:)),Nr,Nt);
 currentMean=reshape(double(current.meanH(k,l,:,:)),Nr,Nt);
 oldEstimate=reshape(double(estimate(k,l,:,:)),Nr,Nt);
 oldC=reshape(double(C(k,l,:,:)),Nr,Nr);
 predicted=B*oldC*B'+projected;
 estimate(k,l,:,:)=reshape(currentMean+B*(oldEstimate-previousMean), ...
     1,1,Nr,Nt);
 for q=1:5
  component=reshape(double(source(k,l,:,:,q)),Nr,Nr);
  component=B*component*B';
  if q==4,component=component+projected;end
  source(k,l,:,:,q)=reshape((component+component')/2, ...
      1,1,Nr,Nr,1);
 end
 crossTerm=reshape(double(crossIn(k,l,:,:)),Nr,Nr);
 crossTerm=B*crossTerm*B';crossTerm=(crossTerm+crossTerm')/2;
 crossOut(k,l,:,:)=reshape(crossTerm,1,1,Nr,Nr);
 reconstructed=reshape(sum(source(k,l,:,:,:),5),Nr,Nr)+crossTerm;
 C(k,l,:,:)=reshape(reconstructed,1,1,Nr,Nr);
 worstReconstruction=max(worstReconstruction, ...
     norm(reconstructed-predicted, ...
     'fro')/max(norm(reconstructed,'fro'),eps));
end,end
diagnostic=struct('previousTti',previous.absoluteTti, ...
    'absoluteTti',current.absoluteTti, ...
    'crossEpochMatched',current.crossEpochMatched, ...
    'conditionTransitionCount',current.conditionTransitionCount, ...
    'previousConditionTransitionCount', ...
        current.previousConditionTransitionCount, ...
    'maximumInnovationPsdProjection',worstProjection, ...
    'maximumPreviousPriorDiagonalJitter',worstPriorJitter, ...
    'maximumSourceReconstructionMismatch',worstReconstruction);
end
