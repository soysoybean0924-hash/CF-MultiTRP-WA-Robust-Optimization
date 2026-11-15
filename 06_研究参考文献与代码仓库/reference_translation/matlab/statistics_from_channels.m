function stats=statistics_from_channels(H,Hhat,p)
% MATLAB translation of main_DNN.py:300-353, given identical H and Hhat.
% H/Hhat have shape M x realizations x K x L. Hhat comes from pilot MMSE;
% this function never substitutes true H for Hhat in the original branch.
[M,N,K,L]=size(H);
assert(isequal(size(Hhat),[M,N,K,L]));
stats=struct();
names={'MR','RZF'};
for kind=1:2
    a=zeros(L,K); b=zeros(L,L,K,K);
    meanEffective=zeros(K,K,L); secondMoment=zeros(K,K,L);
    for n=1:N
        w=zeros(M,K,L);
        for ap=1:L
            V=reshape(Hhat(:,n,:,ap),M,K);
            if kind==1
                W=V;
            else
                W=((p*(V*V')+eye(M))\V);
            end
            w(:,:,ap)=W./vecnorm(W,2,1);
        end
        for ap=1:L
            for user=1:K
                h=H(:,n,user,ap);
                effective=h'*w(:,:,ap);
                a(ap,user)=a(ap,user)+effective(user)/N;
                for stream=1:K
                    meanEffective(user,stream,ap)= ...
                        meanEffective(user,stream,ap)+effective(stream)/N;
                    secondMoment(user,stream,ap)= ...
                        secondMoment(user,stream,ap)+abs(effective(stream))^2/N;
                end
            end
        end
    end
    for user=1:K
        for stream=1:K
            z=reshape(meanEffective(user,stream,:),L,1);
            b(:,:,user,stream)=z*z';
        end
    end
    for ap=1:L
        b(ap,ap,:,:)=reshape(secondMoment(:,:,ap),1,1,K,K);
    end
    stats.(['a_' names{kind}])=abs(a);
    stats.(['B_' names{kind}])=real(b);
end
end
