function test_static_pilot_lmmse()
% Deterministic two-user pilot-contamination fixture, with explicit noise.
H=complex(zeros(2,1,4,1));
H(1,1,1,1)=2; H(2,1,1,1)=4;
R=complex(zeros(4,4,2,1));
R(:,:,1,1)=eye(4); R(:,:,2,1)=eye(4);
pilots=[1;1]; energy=ones(2,1); noise=1;
Np=complex(zeros(1,1,1));
[estimate,Ce,crossCe,observation]=chasi_pilot_lmmse_static( ...
    H,R,pilots,energy,noise,1,Np);
assert(abs(observation(1,1,1)-6)<1e-12);
assert(abs(estimate(1,1,1,1)-2)<1e-12);
assert(abs(estimate(2,1,1,1)-2)<1e-12);
assert(abs(Ce(1,1,1,1)-2/3)<1e-12);
assert(abs(crossCe(1,1,1,2,1)+1/3)<1e-12);
assert(abs(estimate(1,1,3,1))<1e-12);
assert(abs(Ce(3,3,1,1)-1)<1e-12);
pilots=[1;2]; Np=complex(zeros(1,1,2));
[estimate,Ce,crossCe]=chasi_pilot_lmmse_static( ...
    H,R,pilots,energy,noise,1,Np);
assert(abs(estimate(1,1,1,1)-1)<1e-12);
assert(abs(estimate(2,1,1,1)-2)<1e-12);
assert(abs(Ce(1,1,1,1)-0.5)<1e-12);
assert(abs(crossCe(1,1,1,2,1))<1e-12);
energy=4*ones(2,1);
[~,Ce]=chasi_pilot_lmmse_static(H,R,pilots,energy,noise,1,Np);
assert(abs(Ce(1,1,1,1)-0.2)<1e-12);
% Many independent Tx columns form an empirical posterior covariance test.
rng(42,'twister'); n=12000;
R0=0.35.^abs((1:4)'-(1:4)); L0=chol(R0,'lower');
H=complex(zeros(2,1,4,n));
for k=1:2
    z=(randn(4,n)+1i*randn(4,n))/sqrt(2);
    H(k,1,:,:)=reshape(L0*z,1,1,4,n);
    R(:,:,k,1)=R0;
end
Np=(randn(n,1,1)+1i*randn(n,1,1))/sqrt(2);
[estimate,Ce,crossCe]=chasi_pilot_lmmse_static( ...
    H,R,[1;1],ones(2,1),1,1,Np);
errorMatrix=complex(zeros(8,n)); analytic=complex(zeros(8));
for k=1:2
    rows=(k-1)*4+(1:4);
    errorMatrix(rows,:)=reshape(H(k,1,:,:)-estimate(k,1,:,:),4,n);
    for i=1:2
        cols=(i-1)*4+(1:4);
        analytic(rows,cols)=crossCe(:,:,k,i,1);
    end
end
empirical=errorMatrix*errorMatrix'/n;
relative=norm(empirical-analytic,'fro')/norm(analytic,'fro');
assert(relative<0.055);
assert(norm(analytic-analytic','fro')<1e-11);
assert(min(eig(analytic))>-1e-11);
globalNmse=sum(abs(errorMatrix(:)).^2)/sum(abs(H(:)).^2);
theoryNmse=real(trace(analytic))/8;
assert(abs(globalNmse-theoryNmse)<0.035);
for k=1:2
    rows=(k-1)*4+(1:4);
    truthLink=reshape(H(k,1,:,:),4,n);
    empiricalLinkNmse=sum(abs(errorMatrix(rows,:)).^2,'all')/ ...
        sum(abs(truthLink).^2,'all');
    theoryLinkNmse=real(trace(Ce(:,:,k,1)))/real(trace(R0));
    assert(abs(empiricalLinkNmse-theoryLinkNmse)<0.04);
end
assert(abs(Ce(1,1,1,1)-2/3)<1e-12);
assert(real(crossCe(1,1,1,2,1))<-0.3);
% Two observed UE ports on distinct pilots must be solved jointly when
% their prior receive rows are correlated.
H2=complex(zeros(1,1,4,1)); H2(1,1,1,1)=2; H2(1,1,2,1)=4;
R2=reshape(R0,4,4,1,1);
N2=complex(zeros(1,1,2));
[estimate2,Ce2]=chasi_pilot_lmmse_static( ...
    H2,R2,[1 2],[1 1],1,[1 2],N2);
expected2=R0(:,1:2)/(R0(1:2,1:2)+eye(2))*[2;4];
expectedCe2=R0-R0(:,1:2)/(R0(1:2,1:2)+eye(2))*R0(1:2,:);
assert(norm(reshape(estimate2,4,1)-expected2)<1e-12);
assert(norm(reshape(Ce2,4,4)-expectedCe2,'fro')<1e-12);
% A CSI-RS-conditioned prior has a nonzero mean. If truth equals that
% conditional mean and pilot noise is zero, the mean must be preserved.
priorMean=H2;
[estimateMean,CeMean]=chasi_pilot_lmmse_static( ...
    H2,R2,[1 2],[1 1],1,[1 2],N2,priorMean);
assert(norm(estimateMean(:)-priorMean(:))<1e-12);
assert(norm(CeMean(:)-Ce2(:))<1e-12);
fprintf('posterior relative covariance error %.5g, NMSE %.5g vs %.5g\n', ...
    relative,globalNmse,theoryNmse);
fprintf('test_static_pilot_lmmse PASS\n');
end
