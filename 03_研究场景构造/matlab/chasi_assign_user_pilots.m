function userPilot=chasi_assign_user_pilots(beta,tau,seed)
%CHASI_ASSIGN_USER_PILOTS Author pilot_assignment.py:13-24, MATLAB indices.
% beta [UE,TRP] is long-term mean channel power in one common unit.
% tau is the number of orthogonal USER pilot sequences, not SRS ports.
K=size(beta,1); L=size(beta,2);
assert(K>=1 && L>=1 && tau>=1 && tau<=K && tau==floor(tau));
assert(seed>=0 && seed==floor(seed));
assert(all(isfinite(beta(:))) && all(beta(:)>=0));
stream=RandStream('mt19937ar','Seed',seed);
userPilot=zeros(K,1);
userPilot(1:tau)=randperm(stream,tau).';
for k=tau+1:K
    [~,master]=max(beta(k,:));
    interference=zeros(tau,1);
    for t=1:tau
        interference(t)=sum(beta(userPilot(1:k-1)==t,master));
    end
    [~,userPilot(k)]=min(interference);
end
end
