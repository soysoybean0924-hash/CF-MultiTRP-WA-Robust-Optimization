function audit_pilot_statistics_20261029(outputDir)
% Independent, read-only static-estimator audit. No solver or input generator is called.
if nargin<1
    outputDir=fullfile('reproduction','chasi_pilot_mmse','analysis','pilot_statistics_20261029');
end
addpath(fullfile('chasi-matlab-core','chasi_pilot_mmse'));
if ~isfolder(outputDir), mkdir(outputDir); end
fprintf('AUDIT_START seed=20261029 K=2 L=2 Nr=4 TxColumns=4000 checkpoints=[500 2000 4000]\n');
test_pilot_assignment;
test_csi_rs_gaussian_prior;
test_static_measurement_chain;
test_static_pilot_lmmse;

seed=20261029; rng(seed,'twister');
K=2; L=2; Nr=4; N=4000; rows=[1 2];
Rbase=0.35.^abs((1:Nr)'-(1:Nr));
Z=cn([K,L,Nr,N]); Nc=cn([K,L,Nr,N]); Np=cn([N,L,4]);
baseGain=[1 .4; .2 .08];
scenarios={ ...
    struct('name','reuse_E1','gain',baseGain,'pilotId',[1 3;1 3],'energy',1,'srsNoise',10), ...
    struct('name','reuse_E4','gain',baseGain,'pilotId',[1 3;1 3],'energy',4,'srsNoise',10), ...
    struct('name','reuse_noise2p5','gain',baseGain,'pilotId',[1 3;1 3],'energy',1,'srsNoise',2.5), ...
    struct('name','orth_E1','gain',baseGain,'pilotId',[1 3;2 4],'energy',1,'srsNoise',10), ...
    struct('name','changed_links','gain',[1 .15;.6 .04],'pilotId',[1 3;1 3],'energy',1,'srsNoise',10)};
records=cell(1,numel(scenarios)*3); recordIndex=0; indices=[500 2000 N];
for j=1:numel(scenarios)
    sc=scenarios{j}; H=complex(zeros(K,L,Nr,N)); Rrx=complex(zeros(Nr,Nr,K,L));
    csiVariance=zeros(K,L);
    for l=1:L
        for k=1:K
            Rrx(:,:,k,l)=sc.gain(k,l)*Rbase;
            H(k,l,:,:)=reshape(sqrt(sc.gain(k,l))*chol(Rbase,'lower')* ...
                reshape(Z(k,l,:,:),Nr,N),1,1,Nr,N);
            csiVariance(k,l)=sc.gain(k,l)/10^(-3/10);
        end
    end
    energy=sc.energy*ones(K,2);
    [hat,Ce,crossCe,detail]=chasi_static_measurement_chain( ...
        H,Rrx,csiVariance,Nc,sc.pilotId,energy,sc.srsNoise,rows,Np);
    E=H-hat;
    assert(isequal([size(hat,1),size(hat,2),size(hat,3),size(hat,4)],[K,L,Nr,N]));
    assert(isequal(size(Ce),[Nr,Nr,K,L]) && ...
        isequal(size(crossCe),[Nr,Nr,K,K,L]));
    assert(all(isfinite(hat(:))) && all(isfinite(Ce(:))) && ...
        all(isfinite(detail.srsObservation(:))));
    theoryErr=0; theoryTrue=0; hermitianMax=0; minEigen=Inf;
    linkTheoryErr=zeros(K,L); linkTheoryNmse=zeros(K,L);
    for l=1:L
        for k=1:K
            C=Ce(:,:,k,l); R=Rrx(:,:,k,l);
            hermitianMax=max(hermitianMax,norm(C-C','fro'));
            minEigen=min(minEigen,min(real(eig((C+C')/2))));
            linkTheoryErr(k,l)=real(trace(C));
            linkTheoryNmse(k,l)=real(trace(C))/real(trace(R));
            theoryErr=theoryErr+real(trace(C));
            theoryTrue=theoryTrue+real(trace(R));
        end
    end
    assert(hermitianMax<1e-10 && minEigen>-1e-10);
    if strcmp(sc.name,'reuse_E1')
        crossTheory=crossCe(:,:,1,2,1);
        e1=reshape(E(1,1,:,:),Nr,N); e2=reshape(E(2,1,:,:),Nr,N);
        crossEmp=e1*e2'/N;
        fprintf('CROSS reused_UE1_UE2_TRP1 theory11=%.9g empirical11=%.9g+%.9gi froDiff=%.9g\n', ...
            real(crossTheory(1,1)),real(crossEmp(1,1)),imag(crossEmp(1,1)), ...
            norm(crossEmp-crossTheory,'fro'));
    end
    for t=1:numel(indices)
        n=indices(t); e=E(:,:,:,1:n); h=H(:,:,:,1:n); hh=hat(:,:,:,1:n);
        empiricalErr=sum(abs(e(:)).^2)/n;
        empiricalTrue=sum(abs(h(:)).^2)/n;
        empiricalNmse=empiricalErr/empiricalTrue;
        errorMean=mean(reshape(e,[],n),2);
        meanNorm=norm(errorMean);
        correlation=sum(e(:).*conj(hh(:)))/sqrt(sum(abs(e(:)).^2)*sum(abs(hh(:)).^2));
        linkEmpErr=zeros(K,L); linkEmpNmse=zeros(K,L); linkCovRel=zeros(K,L);
        for l=1:L
            for k=1:K
                ek=reshape(e(k,l,:,:),Nr,n); hk=reshape(h(k,l,:,:),Nr,n);
                linkEmpErr(k,l)=sum(abs(ek(:)).^2)/n;
                linkEmpNmse(k,l)=sum(abs(ek(:)).^2)/sum(abs(hk(:)).^2);
                linkCovRel(k,l)=norm(ek*ek'/n-Ce(:,:,k,l),'fro')/norm(Ce(:,:,k,l),'fro');
            end
        end
        rec=struct('scenario',sc.name,'samples',n,'seed',seed, ...
            'pilotEnergy',sc.energy,'srsNoiseVariance',sc.srsNoise, ...
            'csiRsSnrDb',-3,'gains',sc.gain,'pilotId',sc.pilotId, ...
            'theoryMsePerTxColumn',theoryErr,'empiricalMsePerTxColumn',empiricalErr, ...
            'theoryTrueEnergyPerTxColumn',theoryTrue, ...
            'empiricalTrueEnergyPerTxColumn',empiricalTrue, ...
            'theoryNmse',theoryErr/theoryTrue,'empiricalNmse',empiricalNmse, ...
            'errorMeanNorm',meanNorm,'errorHatCorrelationReal',real(correlation), ...
            'errorHatCorrelationImag',imag(correlation), ...
            'linkTheoryMse',linkTheoryErr,'linkEmpiricalMse',linkEmpErr, ...
            'linkTheoryNmse',linkTheoryNmse,'linkEmpiricalNmse',linkEmpNmse, ...
            'linkCovarianceRelativeError',linkCovRel, ...
            'covarianceHermitianError',hermitianMax,'covarianceMinEigenvalue',minEigen);
        recordIndex=recordIndex+1; records{recordIndex}=rec;
        fprintf(['CASE %-19s n=%d E=%.9g/%.9g NMSE=%.9g/%.9g ', ...
            'meanNorm=%.4g corr=%.4g%+.4gi covRelMax=%.4g\n'], ...
            sc.name,n,empiricalErr,theoryErr,empiricalNmse,theoryErr/theoryTrue, ...
            meanNorm,real(correlation),imag(correlation),max(linkCovRel(:)));
    end
    if strcmp(sc.name,'reuse_E1')
        fprintf('NOISE unitNc=%.9g unitNp=%.9g expected=1\n', ...
            mean(abs(Nc(:)).^2),mean(abs(Np(:)).^2));
        fprintf('LINK base UE1_TRP1 MSE=%.9g/%.9g NMSE=%.9g/%.9g; UE2_TRP1 MSE=%.9g/%.9g NMSE=%.9g/%.9g\n', ...
            linkEmpErr(1,1),linkTheoryErr(1,1),linkEmpNmse(1,1),linkTheoryNmse(1,1), ...
            linkEmpErr(2,1),linkTheoryErr(2,1),linkEmpNmse(2,1),linkTheoryNmse(2,1));
    end
end

% Frozen input is read, never generated or overwritten. Inspect three groups.
inputPath=fullfile('reproduction','chasi_pilot_mmse','inputs', ...
    'target_static_gaussian_v1_seed19.mat');
loaded=load(inputPath,'data'); frozen=loaded.data;
assert(strcmp(frozen.modelId,'target_gaussian_static_pilot_lmmse_hypothesis_v1'));
assert(numel(frozen.resources)==17 && frozen.sceneConfig.users==210);
assert(frozen.config.rank==2 && frozen.config.pilotSequences==8 && ...
    frozen.config.pilotMatchedEnergy==1 && frozen.config.pilotSnrDb==-10 && ...
    frozen.config.csiRsSnrDb==-3 && isequal(frozen.config.srsRows(:)',[1 2]));
frozenGroups=[1 9 17]; frozenRecords=cell(1,numel(frozenGroups));
for j=1:numel(frozenGroups)
    g=frozenGroups(j); r=frozen.resources(g); [Kf,Lf,Nrf,Ntf]=size(r.H_true);
    assert(isequal([Lf,Nrf,Ntf],[21 4 64]));
    assert(isequal(size(r.H),size(r.H_true)) && ...
        isequal(size(r.C),[Kf,Lf,Nrf,Nrf]) && isequal(r.C,r.C_stat));
    assert(all(isfinite(r.H(:))) && all(isfinite(r.H_true(:))) && all(isfinite(r.C(:))));
    errorEnergy=sum(abs(double(r.H_true(:))-double(r.H(:))).^2);
    trueEnergy=sum(abs(double(r.H_true(:))).^2);
    savedNmse=frozen.diagnostics(g).nmse;
    recomputed=errorEnergy/trueEnergy;
    assert(abs(recomputed-savedNmse)<1e-6);
    block=reshape(double(r.C(1,1,:,:)),Nrf,Nrf);
    herm=norm(block-block','fro'); mineig=min(real(eig((block+block')/2)));
    assert(herm<1e-6 && mineig>-1e-6);
    frozenRecords{j}=struct('group',g,'scheduledUsers',Kf, ...
        'trps',Lf,'rxPorts',Nrf,'txPorts',Ntf,'ageMs',r.age_ms, ...
        'rbCount',r.rb_count,'resourceWeightNu',r.resourceWeightNu, ...
        'savedNmse',savedNmse,'recomputedNmse',recomputed, ...
        'covarianceHermitianErrorFirstLink',herm, ...
        'covarianceMinEigenvalueFirstLink',mineig);
    fprintf('FROZEN g=%d K=%d H=[%d %d %d %d] ageMs=%.9g nu=%.9g NMSE=%.9g/%.9g CeMinEig=%.4g\n', ...
        g,Kf,Kf,Lf,Nrf,Ntf,r.age_ms,r.resourceWeightNu,recomputed,savedNmse,mineig);
end
result=struct('model','static calibrated-port Gaussian hypothesis', ...
    'seed',seed,'rxCorrelationRho',0.35,'rxPorts',Nr,'txColumns',N, ...
    'checkpoints',indices,'cases',[records{:}],'frozenGroups',[frozenRecords{:}], ...
    'frozenInputPath',inputPath,'frozenInputSha256', ...
    strtrim(fileread([inputPath '.sha256'])));
path=fullfile(outputDir,'summary.json'); fid=fopen(path,'wt'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(result,'PrettyPrint',true)); fclose(fid);
fprintf('AUDIT_DONE cases=%d frozenGroups=%d output=%s\n',numel(records),numel(frozenRecords),path);
end

function z=cn(sz)
z=(randn(sz)+1i*randn(sz))/sqrt(2);
end
