function test_pilot_assignment()
% First tau users get a random permutation; later users minimize local
% assigned beta at the master TRP, following pilot_assignment.py.
beta=[5 1;1 5;9 1;1 9;8 1]; tau=2; seed=771;
before=rng;
userPilot=chasi_assign_user_pilots(beta,tau,seed);
after=rng;
assert(isequal(before,after));
assert(isequal(sort(userPilot(1:tau)),(1:tau)'));
assert(isequal(userPilot,chasi_assign_user_pilots(beta,tau,seed)));
for k=tau+1:size(beta,1)
    [~,master]=max(beta(k,:));
    interference=zeros(tau,1);
    for t=1:tau
        interference(t)=sum(beta(userPilot(1:k-1)==t,master));
    end
    [~,best]=min(interference);
    assert(userPilot(k)==best);
end
ports=chasi_expand_port_pilots(userPilot,tau,2);
assert(isequal(ports(:,1),userPilot));
assert(isequal(ports(:,2),userPilot+tau));
assert(all(ports(:,1)~=ports(:,2)));
fprintf('test_pilot_assignment PASS\n');
end
