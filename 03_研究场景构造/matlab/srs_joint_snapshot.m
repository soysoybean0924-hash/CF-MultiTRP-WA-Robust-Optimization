function out=srs_joint_snapshot(H0,R,pilotId,energy,noise,z,ageMs,tcMs,innovation,portRows,saveCross)
%SRS_JOINT_SNAPSHOT Shared-observation conditional Gaussian snapshot.
% H0 [K,L,4,Nt], R [4,4,K,L], z [Nt,L,T] unit CN, innovation same
% shape as H0 and independent of H0,z. One pilot index per UE/UL port.
% This is an effective 2-port to 4-row selection hypothesis, not calibration.
K=size(H0,1); L=size(H0,2); Nr=size(H0,3); Nt=size(H0,4);
T=size(z,3); P=size(pilotId,2);
if nargin<11, saveCross=true; end
checkShape('H0',H0,[K,L,Nr,Nt]);
if K<1 || L<1 || Nt<1 || Nr~=4
    error('srs_joint_snapshot:ChannelDimensions', ...
        'H0 must have K>=1, L>=1, Nr=4, Nt>=1; actual size %s.', ...
        mat2str(size(H0)));
end
checkShape('pilotId',pilotId,[K,P]);
if P~=2
    error('srs_joint_snapshot:UnsupportedPortCount', ...
        'pilotId must be K-by-2 for this 2T4R model; expected [%d 2], actual %s.', ...
        K,mat2str(size(pilotId)));
end
checkShape('energy',energy,[K,P]);
checkShape('R',R,[Nr,Nr,K,L]);
checkShape('z',z,[Nt,L,T]);
checkShape('innovation',innovation,[K,L,Nr,Nt]);
if numel(portRows)~=P || any(portRows(:)<1) || ...
        any(portRows(:)>Nr) || any(portRows(:)~=floor(portRows(:))) || ...
        numel(unique(portRows(:)))~=P
    error('srs_joint_snapshot:PortRows', ...
        'portRows must select %d distinct rows in 1..%d; actual %s.', ...
        P,Nr,mat2str(portRows));
end
if any(pilotId(:)<1 | pilotId(:)>T | pilotId(:)~=floor(pilotId(:)))
    error('srs_joint_snapshot:PilotIdRange', ...
        'pilotId must contain integer indices in 1..%d; actual size %s.', ...
        T,mat2str(size(pilotId)));
end
if ~isreal(energy) || any(~isfinite(energy(:))) || any(energy(:)<0)
    error('srs_joint_snapshot:EnergyValue', ...
        'energy must be real, finite, nonnegative; actual size %s.', ...
        mat2str(size(energy)));
end
if ~isscalar(noise) || ~isreal(noise) || ~isfinite(noise) || noise<=0
    error('srs_joint_snapshot:NoiseValue', ...
        'noise must be a positive finite real scalar; actual size %s.',mat2str(size(noise)));
end
if ~isscalar(tcMs) || ~isreal(tcMs) || ~isfinite(tcMs) || tcMs<=0
    error('srs_joint_snapshot:CoherenceTime', ...
        'tcMs must be a positive finite real scalar; actual size %s.',mat2str(size(tcMs)));
end
if ~isscalar(ageMs) || ~isreal(ageMs) || ~isfinite(ageMs) || ageMs<0
    error('srs_joint_snapshot:AgeValue', ...
        'ageMs must be a nonnegative finite real scalar; actual size %s.',mat2str(size(ageMs)));
end
a=exp(-ageMs/tcMs); % project snapshot assumption; not Jakes trajectory
out.H_true=a*H0+sqrt(1-a^2)*innovation;
out.H_hat=complex(zeros(size(H0)));
out.C=complex(zeros(Nr,Nr,K,L));
if saveCross, out.crossUE=complex(zeros(Nr,Nr,K,K,L));
else, out.crossUE=[]; end
out.Y=complex(zeros(T,Nt,L));
for l=1:L
    M=complex(zeros(T,K*Nr)); Q=complex(zeros(K*Nr));
    X=complex(zeros(K*Nr,Nt));
    for k=1:K
        ix=(k-1)*Nr+(1:Nr);
        X(ix,:)=reshape(H0(k,l,:,:),Nr,Nt);
        Q(ix,ix)=R(:,:,k,l);
        for p=1:P
            M(pilotId(k,p),ix(portRows(p)))=sqrt(energy(k,p));
        end
    end
    Y=M*X+sqrt(noise)*reshape(z(:,l,:),Nt,T).';
    S=M*Q*M'+noise*eye(T);
    B=Q*M'/S;
    Xhat=B*Y;
    C0=Q-B*M*Q; C0=(C0+C0')/2;
    out.Y(:,:,l)=Y;
    for k=1:K
        ix=(k-1)*Nr+(1:Nr);
        out.H_hat(k,l,:,:)=reshape(a*Xhat(ix,:),1,1,Nr,Nt);
        out.C(:,:,k,l)=a^2*C0(ix,ix)+(1-a^2)*R(:,:,k,l);
        if saveCross, for j=1:K
            jx=(j-1)*Nr+(1:Nr);
            out.crossUE(:,:,k,j,l)=a^2*C0(ix,jx);
            if j==k
                out.crossUE(:,:,k,j,l)=out.crossUE(:,:,k,j,l)+ ...
                    (1-a^2)*R(:,:,k,l);
            end
        end, end
    end
end
out.a=a;
out.modelId='srs_joint_snapshot_v2';
end

function checkShape(name,x,expected)
% MATLAB drops trailing singleton dimensions from size(x). Per-axis checks
% retain the contract for K=1/L=1/T=1; numel rejects extra non-singleton axes.
actual=zeros(1,numel(expected));
for d=1:numel(expected), actual(d)=size(x,d); end
if ~isequal(actual,expected) || numel(x)~=prod(expected)
    error(['srs_joint_snapshot:' name 'Size'], ...
        '%s expected dimensions %s; actual size %s (per-axis %s).', ...
        name,mat2str(expected),mat2str(size(x)),mat2str(actual));
end
end
