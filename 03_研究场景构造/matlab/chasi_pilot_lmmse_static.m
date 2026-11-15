function [Hhat,Ce,crossCe,observation]=chasi_pilot_lmmse_static( ...
    Htrue,Rrx,pilotId,pilotEnergy,noiseVariance,srsRows,Np,priorMean)
%CHASI_PILOT_LMMSE_STATIC One static Gaussian pilot observation per TRP.
% Htrue: [UE,TRP,DL-Rx,TRP-Tx]. Rrx: [DL-Rx,DL-Rx,UE,TRP]
% pilotId/energy: [UE,SRS-port], energy=tau_p * pilot power in the
% same noise-normalized units as Rrx. srsRows maps observed SRS ports to
% DL rows and is an EXPLICIT reciprocal-port hypothesis, not a 2Tx/4Rx
% identity. Np: independent unit-variance CN noise [TRP-Tx,TRP,pilot].
% Ce: marginal error covariance [DL-Rx,DL-Rx,UE,TRP] for each Tx column.
% crossCe: joint error blocks [DL-Rx,DL-Rx,UE,UE,TRP]. Full vec(H)
% covariance uses these blocks and I_Tx (independent Tx-column prior).
% For jointly circular Gaussian H and noise this is MMSE; otherwise LMMSE.
% priorMean may be the per-link conditional mean after a validated CSI-RS
% observation. The default is zero. Rrx then denotes conditional covariance.
if nargin<8, priorMean=complex(zeros(size(Htrue))); end
K=size(Htrue,1); L=size(Htrue,2); Nr=size(Htrue,3); Nt=size(Htrue,4);
P=size(pilotId,2); T=size(Np,3);
assert(isequal([size(Rrx,1),size(Rrx,2),size(Rrx,3),size(Rrx,4)], ...
    [Nr,Nr,K,L]));
assert(isequal(size(pilotId),[K,P]) && isequal(size(pilotEnergy),[K,P]));
assert(isequal([size(Np,1),size(Np,2),size(Np,3)],[Nt,L,T]));
assert(isequal([size(priorMean,1),size(priorMean,2), ...
    size(priorMean,3),size(priorMean,4)],[K,L,Nr,Nt]));
assert(numel(srsRows)==P && all(srsRows>=1) && all(srsRows<=Nr));
assert(all(pilotId(:)==floor(pilotId(:))) && all(pilotId(:)>=1) && ...
    all(pilotId(:)<=T));
assert(all(isfinite(pilotEnergy(:))) && all(pilotEnergy(:)>=0));
assert(isscalar(noiseVariance) && isfinite(noiseVariance) && noiseVariance>0);
Hhat=complex(zeros(K,L,Nr,Nt));
Ce=complex(zeros(Nr,Nr,K,L));
crossCe=complex(zeros(Nr,Nr,K,K,L));
observation=complex(zeros(Nt,L,T));
for l=1:L
    A=complex(zeros(T,K*Nr));
    B=complex(zeros(K*Nr,T));
    Hstack=complex(zeros(K*Nr,Nt));
    Mstack=complex(zeros(K*Nr,Nt));
    for k=1:K
        rows=(k-1)*Nr+(1:Nr);
        Hstack(rows,:)=reshape(Htrue(k,l,:,:),Nr,Nt);
        Mstack(rows,:)=reshape(priorMean(k,l,:,:),Nr,Nt);
        for port=1:P
            t=pilotId(k,port);
            A(t,rows(srsRows(port)))=A(t,rows(srsRows(port)))+ ...
                sqrt(pilotEnergy(k,port));
        end
        R=Rrx(:,:,k,l);
        assert(norm(R-R','fro')<=1e-9*max(1,norm(R,'fro')));
        assert(min(eig((R+R')/2))>=-1e-9*max(1,norm(R,'fro')));
        B(rows,:)=R*A(:,rows)';
    end
    Y=A*Hstack+sqrt(noiseVariance)*reshape(Np(:,l,:),Nt,T).';
    observation(:,l,:)=reshape(Y.',Nt,1,T);
    S=(A*B+noiseVariance*eye(T));
    S=(S+S')/2;
    gain=B/S;
    estimate=Mstack+gain*(Y-A*Mstack);
    for k=1:K
        rows=(k-1)*Nr+(1:Nr);
        Hhat(k,l,:,:)=reshape(estimate(rows,:),1,1,Nr,Nt);
        for i=1:K
            cols=(i-1)*Nr+(1:Nr);
            block=-gain(rows,:)*B(cols,:)';
            if i==k
                block=block+Rrx(:,:,k,l);
                Ce(:,:,k,l)=(block+block')/2;
            end
            crossCe(:,:,k,i,l)=block;
        end
    end
end
end
