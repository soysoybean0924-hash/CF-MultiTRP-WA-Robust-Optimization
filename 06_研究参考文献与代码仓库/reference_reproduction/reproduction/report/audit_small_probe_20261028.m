function audit_small_probe_20261028()
% Read-only, deterministic, small audit probes. No production inputs are written.
addpath('chasi-matlab-core/chasi_pilot_mmse');
addpath('D:/GitHubProjects/CF-MultiTRP-Optimization/matlab/chasi');
rng(20261028,'twister');
Nr=4; Nt=3; R=0.35.^abs((1:Nr)'-(1:Nr));
H=reshape(chol(R,'lower')*((randn(Nr,Nt)+1i*randn(Nr,Nt))/sqrt(2)),1,1,Nr,Nt);
Nc=(randn(1,1,Nr,Nt)+1i*randn(1,1,Nr,Nt))/sqrt(2);
Np=(randn(Nt,1,2)+1i*randn(Nt,1,2))/sqrt(2);
[hat,Ce,~,detail]=chasi_static_measurement_chain(H,reshape(R,Nr,Nr,1,1),1,Nc,[1 2],[1 1],1,[1 2],Np);
err=H-hat;
fprintf('SAMPLE H=[%d %d] normH=%.9g normHat=%.9g normErr=%.9g sampleNMSE=%.9g CeTrace=%.9g ageMs=0\n',Nr,Nt,norm(H(:)),norm(hat(:)),norm(err(:)),sum(abs(err(:)).^2)/sum(abs(H(:)).^2),real(trace(Ce(:,:,1,1))));
assert(all(isfinite(detail.srsObservation(:))));
M=2000; L=chol(R,'lower'); Z=(randn(Nr,Nt*M)+1i*randn(Nr,Nt*M))/sqrt(2);
Hs=reshape(L*Z,1,1,Nr,Nt*M);
Ncs=(randn(size(Hs))+1i*randn(size(Hs)))/sqrt(2);
Nps=(randn(Nt*M,1,2)+1i*randn(Nt*M,1,2))/sqrt(2);
[hats,Ces]=chasi_static_measurement_chain(Hs,reshape(R,Nr,Nr,1,1),1,Ncs,[1 2],[1 1],1,[1 2],Nps);
es=Hs-hats;
emp=sum(abs(es(:)).^2)/sum(abs(Hs(:)).^2);
theory=real(trace(Ces(:,:,1,1)))/real(trace(R));
fprintf('REPEATED M=%d empiricalNMSE=%.9g theoryNMSE=%.9g absDiff=%.9g\n',M,emp,theory,abs(emp-theory));
assert(abs(emp-theory)<0.035);
for rank=[1 2]
    K=rank; LL=2; G=2; resources=struct([]);
    weights=[0.25 0.75];
    for g=1:G
        h=(randn(K,LL,Nr,Nt)+1i*randn(K,LL,Nr,Nt))/sqrt(2);
        mask=false(LL,K); mask(1,:)=true; if rank==2, mask(2,2)=true; end
        resources(g).H=h;
        resources(g).H_true=h;
        resources(g).C=zeros(K,LL,Nr,Nr);
        resources(g).C_stat=resources(g).C;
        resources(g).servingMask=mask;
        resources(g).ueIds=(1:K)';
        resources(g).resourceWeightNu=weights(g);
    end
    cfg=struct('rank',rank,'budget',0.6,'outer_iterations',2,'admm_iterations',100,'cg_iterations',100,'branch','matrix_nominal','users',K,'eps_abs',1e-7,'eps_rel',1e-6);
    nominal=chasi_matrix_wmmse(resources,cfg);
    cfg.branch='matrix_statistical'; statistical=chasi_matrix_wmmse(resources,cfg);
    power=zeros(LL,1); leak=0;
    for g=1:G
        x=nominal.precoder{g}; nu=resources(g).resourceWeightNu;
        for l=1:LL, power(l)=power(l)+nu*sum(abs(x(l,:,:,:)).^2,'all'); end
        leak=max(leak,max(abs(x.*reshape(~resources(g).servingMask,LL,K,1,1)),[],'all'));
    end
    fprintf('NATIVE rank=%d K=%d maxPower=%.9g leak=%.3g zeroCObjectiveDiff=%.3g outerDelta=%.9g innerConverged=%d\n',rank,K,max(power),leak,abs(nominal.optimization_objective-statistical.optimization_objective),nominal.optimization_objective-nominal.trajectory(1).objective,nominal.trajectory(end).previous_admm_converged);
    assert(max(power)<=cfg.budget+1e-7 && leak==0);
    assert(abs(nominal.optimization_objective-statistical.optimization_objective)<1e-9);
end
end
