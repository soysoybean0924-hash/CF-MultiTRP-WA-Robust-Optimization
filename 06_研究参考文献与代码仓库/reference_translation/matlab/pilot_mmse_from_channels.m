function Hhat=pilot_mmse_from_channels(H,CorrR,Np,pilotIndex,p,tau_p)
% Translation of main_DNN.py:270-298 with the SAME exported H,R,Np,pilots.
% pilotIndex is Python's zero-based pilot vector; Hhat remains nonideal CSI.
[M,N,K,L]=size(H);
assert(isequal(size(CorrR),[M,M,K,L]));
assert(isequal(size(Np),[M,N,L,tau_p]));
pilot=double(pilotIndex(:))+1;
Hhat=complex(zeros(M,N,K,L));
for l=1:L
    for t=1:tau_p
        users=find(pilot==t);
        yp=sqrt(p*tau_p)*sum(H(:,:,users,l),3)+Np(:,:,l,t);
        Psi=p*tau_p*sum(CorrR(:,:,users,l),3)+eye(M);
        for index=1:numel(users)
            k=users(index);
            RPsi=CorrR(:,:,k,l)/Psi;
            Hhat(:,:,k,l)=sqrt(p*tau_p)*RPsi*yp;
        end
    end
end
end
