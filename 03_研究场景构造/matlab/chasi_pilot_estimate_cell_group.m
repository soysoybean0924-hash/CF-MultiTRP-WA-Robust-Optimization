function [Hhat,Ce,crossCe,observation]=chasi_pilot_estimate_cell_group( ...
    Htrue,Rrx,pilotId,pilotEnergy,noiseVariance,srsRows,Np,priorMean)
%CHASI_PILOT_ESTIMATE_CELL_GROUP Explicit one-RB-group solver adapter.
% The source WMMSE solver takes H{k,l,g} as [DL-Rx,TRP-Tx]. This adapter
% does not set a scene, seed, prior, pilot plan, or physical port mapping.
assert(iscell(Htrue) && size(Htrue,3)==1);
K=size(Htrue,1); L=size(Htrue,2);
assert(K>=1 && L>=1 && ~isempty(Htrue{1,1,1}));
[Nr,Nt]=size(Htrue{1,1,1});
H=complex(zeros(K,L,Nr,Nt));
M=complex(zeros(K,L,Nr,Nt));
if nargin<8, priorMean=[]; end
assert(isempty(priorMean) || (iscell(priorMean) && isequal(size(priorMean),size(Htrue))));
for l=1:L
    for k=1:K
        h=Htrue{k,l,1};
        assert(isequal(size(h),[Nr,Nt]) && all(isfinite(h(:))));
        H(k,l,:,:)=reshape(double(h),1,1,Nr,Nt);
        if ~isempty(priorMean)
            m=priorMean{k,l,1};
            assert(isequal(size(m),[Nr,Nt]) && all(isfinite(m(:))));
            M(k,l,:,:)=reshape(double(m),1,1,Nr,Nt);
        end
    end
end
[dense,Ce,crossCe,observation]=chasi_pilot_lmmse_static( ...
    H,Rrx,pilotId,pilotEnergy,noiseVariance,srsRows,Np,M);
Hhat=cell(K,L,1);
for l=1:L
    for k=1:K
        Hhat{k,l,1}=reshape(dense(k,l,:,:),Nr,Nt);
    end
end
end
