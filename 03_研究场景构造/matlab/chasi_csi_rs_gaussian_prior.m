function [priorMean,priorCov,observation]=chasi_csi_rs_gaussian_prior( ...
    Htrue,Rrx,noiseVariance,Nc)
%CHASI_CSI_RS_GAUSSIAN_PRIOR Hypothetical calibrated four-row CSI-RS input.
% Htrue,Nc: [UE,TRP,DL-Rx,TRP-Tx]; Nc independent unit-variance CN.
% Rrx: [DL-Rx,DL-Rx,UE,TRP], zero-mean prior per independent Tx column.
% noiseVariance: scalar or [UE,TRP] in the SAME normalization as Rrx.
% The output conditions on Htrue+sqrt(noiseVariance)*Nc. Real DL CSI-RS
% feedback, delay, receiver noise, and port calibration are not modeled.
K=size(Htrue,1); L=size(Htrue,2); Nr=size(Htrue,3); Nt=size(Htrue,4);
assert(isequal([size(Nc,1),size(Nc,2),size(Nc,3),size(Nc,4)], ...
    [K,L,Nr,Nt]));
assert(isequal([size(Rrx,1),size(Rrx,2),size(Rrx,3),size(Rrx,4)], ...
    [Nr,Nr,K,L]));
if isscalar(noiseVariance)
    noiseVariance=repmat(noiseVariance,K,L);
end
assert(isequal(size(noiseVariance),[K,L]));
assert(all(isfinite(noiseVariance(:))) && all(noiseVariance(:)>0));
priorMean=complex(zeros(K,L,Nr,Nt));
priorCov=complex(zeros(Nr,Nr,K,L));
observation=complex(zeros(K,L,Nr,Nt));
for l=1:L
    for k=1:K
        R=Rrx(:,:,k,l);
        assert(norm(R-R','fro')<=1e-9*max(1,norm(R,'fro')));
        assert(min(eig((R+R')/2))>=-1e-9*max(1,norm(R,'fro')));
        variance=noiseVariance(k,l);
        y=reshape(Htrue(k,l,:,:),Nr,Nt)+ ...
            sqrt(variance)*reshape(Nc(k,l,:,:),Nr,Nt);
        gain=R/(R+variance*eye(Nr));
        mean=gain*y;
        C=R-gain*R;
        priorMean(k,l,:,:)=reshape(mean,1,1,Nr,Nt);
        priorCov(:,:,k,l)=(C+C')/2;
        observation(k,l,:,:)=reshape(y,1,1,Nr,Nt);
    end
end
end
