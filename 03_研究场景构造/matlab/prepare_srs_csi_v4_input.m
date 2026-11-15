function info=prepare_srs_csi_v4_input(root,configPath,inputPath)
% Independent v4 input: finite SRS posterior, chronological CSI-RS updates.
assert(~isfile(inputPath) && isfolder(root));
addpath(fullfile(root,'03_研究场景构造','matlab'),fullfile(root,'04_研究方法_WA运行实现','matlab'));
c=jsondecode(fileread(configPath));
assert(strcmp(c.modelId,'srs_csi_joint_v4'));
assert(c.prbUtilization==.5 && strcmp(c.scheduleMode,'legacy_random_full_groups'));
assert(c.measurementBeforeDecisionAtSameTti && ...
    strcmp(c.srsSenderPolicy,'all_ues_at_each_groups_hop_time'));
if ~isfield(c,'speedKmh'), c.speedKmh=3; end
cfg=chasi_config(struct('seed',c.sceneSeed,'usersPerCell',c.usersPerCell, ...
    'rank',c.rank,'powerWattPerTrp',c.powerWattPerTrp, ...
    'prbUtilization',c.prbUtilization,'scheduleMode',c.scheduleMode, ...
    'csiRsSnrDb',c.csiRsSnrDb,'csiRsPeriodTti',c.csiRsPeriodTti, ...
    'speedKmh',c.speedKmh, ...
    'channelMode',c.channelMode,'frequencyGhz',c.frequencyGhz, ...
    'noiseWattPerRb',c.noiseWattPerRb, ...
    'serviceTrps',2,'srsMode','pilot_reuse_gaussian', ...
    'srsPilotSequences',c.uePilotLabels,'decisionTti',c.decisionTti, ...
    'durationTti',c.decisionTti,'ttiMs',c.ttiMs));
assert(cfg.srsHopIntervalTti==c.hopIntervalTti && ...
    cfg.srsPeriodTti==c.revisitTti && cfg.srsHops==c.resourceGroups && ...
    c.pilotSequenceResources==2*c.uePilotLabels);
% The legacy estimator runs inside scene generation, but none of its H/C,
% age/noise diagnostics are transferred into this input.
scene=chasi_generate_scene(cfg);
K=cfg.users; L=cfg.trps; Nr=4; Nt=cfg.txPorts; G=cfg.srsHops;
if isfield(c,'validationGroups')
    assert(c.validationGroups>=1 && c.validationGroups<=G);
    G=c.validationGroups;
end
d=zeros(K,L);
for k=1:K
    d(k,:)=sqrt(sum((scene.trp_xy-scene.ue_xy(k,:)).^2,2));
end
% Source gain includes its sector pattern, path loss, optional shadow and
% noise normalization. Distance-only is a separate declared ablation.
switch char(c.propagationPrior)
    case 'source_large_scale_gain'
        assert(c.rxCorrelationRho==0.35, ...
            'Source scene receive correlation is fixed at 0.35.');
        gain=double(scene.largeScaleGain);
        anchor=double(scene.largeScaleHomeMedian);
    case 'distance_power_law_ablation'
        anchor=double(scene.largeScaleHomeMedian);
        homeDistance=d(sub2ind(size(d),(1:K)',double(scene.home(:))));
        referenceDistance=median(homeDistance);
        gain=anchor*(referenceDistance./max(d,25)).^c.distanceExponent;
    otherwise
        error('Unknown propagation prior');
end
assert(isequal(size(gain),[K,L]) && all(isfinite(gain(:))) && ...
    all(gain(:)>0) && isfinite(anchor) && anchor>0);
Rbase=c.rxCorrelationRho.^abs((1:Nr)'-(1:Nr)); cholRx=chol(Rbase,'lower');
R=complex(zeros(Nr,Nr,K,L));
for l=1:L, for k=1:K, R(:,:,k,l)=gain(k,l)*Rbase; end, end
pilotId=[mod((0:K-1)',c.uePilotLabels)+1, ...
    mod((0:K-1)',c.uePilotLabels)+1+c.uePilotLabels];
energy=ones(K,2)*c.pilotMatchedEnergy;
noise=anchor/10^(c.pilotSnrDb/10);
stream=RandStream('mt19937ar','Seed',c.measurementSeed);
csiStream=RandStream('mt19937ar','Seed',c.measurementSeed+77);
resources=scene.resources(1:G);
for g=1:G
    H0=complex(zeros(K,L,Nr,Nt)); innovation=H0;
    for l=1:L, for k=1:K
        H0(k,l,:,:)=reshape(sqrt(gain(k,l))*cholRx*cn(stream,Nr,Nt),1,1,Nr,Nt);
        innovation(k,l,:,:)=reshape(sqrt(gain(k,l))*cholRx*cn(stream,Nr,Nt),1,1,Nr,Nt);
    end, end
    z=cn(stream,Nt,L,2*c.uePilotLabels);
    firstTti=(g-1)*cfg.srsHopIntervalTti;
    measurementTti=firstTti+floor((c.decisionTti-firstTti)/cfg.srsPeriodTti)*cfg.srsPeriodTti;
    ageMs=(c.decisionTti-measurementTti)*c.ttiMs;
    assert(ageMs>=0);
    if strcmp(c.modelId,'srs_csi_joint_v4')
        snap=srs_overlap_snapshot(H0,R,pilotId,energy,noise,z,0, ...
            c.ageCoherenceMs,innovation,c.ulPortToDlRowsOneBased(:)', ...
            resources(g).srs_rb_count,c.speedKmh,c.frequencyGhz);
    elseif strcmp(c.modelId,'srs_overlap_snapshot_v1')
        snap=srs_overlap_snapshot(H0,R,pilotId,energy,noise,z,ageMs, ...
            c.ageCoherenceMs,innovation,c.ulPortToDlRowsOneBased(:)', ...
            resources(g).srs_rb_count);
    else
        snap=srs_joint_snapshot(H0,R,pilotId,energy,noise,z,ageMs, ...
            c.ageCoherenceMs,innovation,c.ulPortToDlRowsOneBased(:)',false);
    end
    ids=double(resources(g).ueIds(:));
    for local=1:numel(ids)
        if nnz(resources(g).servingMask(:,local))<2
            [~,nearest]=sort(d(ids(local),:));
            second=nearest(find(nearest~=scene.home(ids(local)),1));
            resources(g).servingMask(second,local)=true;
        end
    end
    resources(g).H_true=single(snap.H_true(ids,:,:,:));
    resources(g).H=single(snap.H_hat(ids,:,:,:));
    resources(g).C=single(permute(snap.C(:,:,ids,:),[3,4,1,2]));
    resources(g).C_stat=resources(g).C;
    resources(g).prior_gain=single(gain(ids,:));
    resources(g).omitted_gain=zeros(numel(ids),L,'single');
    resources(g).age_ms=ageMs;
    resources(g).pilot_ids=pilotId(ids,:);
    resources(g).co_pilot_indices=arrayfun( ...
        @(local)find(pilotId(:,1)==pilotId(ids(local),1)), ...
        (1:numel(ids))','UniformOutput',false);
    resources(g).srs_noise_variance=noise;
    if isfield(snap,'sourceRisk')
        resources(g).srs_source_risk=single(snap.sourceRisk(ids,:,:));
    end
    resources(g).srs_rows=c.ulPortToDlRowsOneBased(:)';
    resources(g).srs_diagnostic=struct('modelId',c.modelId, ...
        'measurementTti',measurementTti,'decisionTti',c.decisionTti, ...
        'measurementTimeMs',measurementTti*c.ttiMs, ...
        'decisionTimeMs',c.decisionTti*c.ttiMs, ...
        'ageMs',ageMs,'a',snap.a,'noiseVariance',noise, ...
        'srsSenderIds',(1:K)','scheduledIds',ids, ...
        'crossRepresentation','priorRx + pilotId + energy + noise; reconstruct Q-QM''/S*M*Q');
    if isfield(snap,'gram')
        resources(g).srs_diagnostic.pilotGram=snap.gram;
        resources(g).srs_diagnostic.noiseCov=snap.noiseCov;
        resources(g).srs_diagnostic.ageModel=snap.ageModel;
        resources(g).srs_diagnostic.speedKmh=c.speedKmh;
        resources(g).srs_diagnostic.frequencyGhz=c.frequencyGhz;
    end
    % Attribute the SRS posterior before any CSI-RS update. The v3 source
    % function is used only at this exact SRS-time seam (a=1).
    resources(g).srs_diagnostic.a=1;
    priorData=struct('modelId','srs_overlap_jakes_v3', ...
        'resources',resources(g),'pilotId',pilotId,'pilotEnergy',energy, ...
        'priorRx',R,'pilotNoiseVariance',noise);
    source=srs_overlap_source_statistics(priorData,1);
    sourceCov=double(source.covariance);
    currentTti=measurementTti;
    csiEvents=(floor(measurementTti/c.csiRsPeriodTti)+1)*c.csiRsPeriodTti: ...
        c.csiRsPeriodTti:c.decisionTti;
    est=double(resources(g).H);truth=double(H0(ids,:,:,:));
    cov=double(resources(g).C);
    for eventTti=csiEvents
        [est,truth,cov,sourceCov]=advance(est,truth,cov,sourceCov, ...
            R,ids,eventTti-currentTti,c.ttiMs,c.speedKmh,c.frequencyGhz,csiStream);
        for local=1:numel(ids),for l=1:L
            P=reshape(cov(local,l,:,:),Nr,Nr);
            obsVar=gain(ids(local),l)/10^(c.csiRsSnrDb/10);
            J=P/(P+obsVar*eye(Nr));F=eye(Nr)-J;
            y=reshape(truth(local,l,:,:),Nr,Nt)+sqrt(obsVar)*cn(csiStream,Nr,Nt);
            mu=reshape(est(local,l,:,:),Nr,Nt);
            est(local,l,:,:)=reshape(mu+J*(y-mu),1,1,Nr,Nt);
            cov(local,l,:,:)=reshape(F*P*F'+J*(obsVar*eye(Nr))*J',1,1,Nr,Nr);
            for q=1:4
                Q=reshape(sourceCov(local,l,:,:,q),Nr,Nr);
                sourceCov(local,l,:,:,q)=reshape(F*Q*F',1,1,Nr,Nr,1);
            end
            Q=reshape(sourceCov(local,l,:,:,3),Nr,Nr)+J*(obsVar*eye(Nr))*J';
            sourceCov(local,l,:,:,3)=reshape(Q,1,1,Nr,Nr,1);
        end,end
        currentTti=eventTti;
    end
    [est,truth,cov,sourceCov]=advance(est,truth,cov,sourceCov, ...
        R,ids,c.decisionTti-currentTti,c.ttiMs,c.speedKmh,c.frequencyGhz,csiStream);
    resources(g).H=single(est);resources(g).H_true=single(truth);
    resources(g).C=single(cov);resources(g).C_stat=resources(g).C;
    resources(g).sourceCov=single(sourceCov);
    resources(g).csi_diagnostic=struct('eventTti',csiEvents, ...
        'srsMeasurementTti',measurementTti,'decisionTti',c.decisionTti, ...
        'lastCsiTti',currentTti,'csiAgeMs',(c.decisionTti-currentTti)*c.ttiMs);
end
% Native v4 runner calls the same deterministic feasible initialization
% for every branch; the historical common-initial helper assumes the old
% resource occupancy and fails its half-budget assertion at 50% PRB.
initialPrecoder=cell(1,G);initialPower=zeros(L,1);
if ismember(c.modelId,{'srs_overlap_snapshot_v1','srs_overlap_jakes_v3','srs_csi_joint_v4'})
    snapshotFile='srs_overlap_snapshot.m';
else
    snapshotFile='srs_joint_snapshot.m';
end
data=struct('modelId',c.modelId,'config',c, ...
    'configPath',configPath, ...
    'parameterSource','finite_sequence_matched_filter_analytic_project_extension', ...
    'sceneSeed',c.sceneSeed,'measurementSeed',c.measurementSeed, ...
    'trafficSeed',c.sceneSeed,'initialSeed',c.initialSeed, ...
    'configSha256',chasi_sha256_file(configPath), ...
    'sceneConfig',cfg,'resources',resources, ...
    'schemaVersion',4,'generatedGroups',G,'fullScaleInput',G==cfg.srsHops, ...
    'prbRequested',c.prbUtilization, ...
    'prbEffectiveByTrp',double(scene.active_prb_ratio), ...
    'prbEffective',mean(double(scene.active_prb_ratio)), ...
    'propagationPrior',c.propagationPrior, ...
    'srsSenderIds',(1:K)','pilotId',pilotId,'pilotEnergy',energy, ...
    'priorRx',R,'pilotNoiseVariance',noise, ...
    'initialPrecoder',{initialPrecoder},'initialPower',initialPower, ...
    'sourceCodeSha256',chasi_sha256_file([mfilename('fullpath') '.m']), ...
    'snapshotCodeSha256',chasi_sha256_file(fullfile(root,'03_研究场景构造','matlab',snapshotFile)));
data.sceneGeneratorSha256=chasi_sha256_file(fullfile(root,'03_研究场景构造','matlab','chasi_generate_scene.m'));
data.sceneConfigSha256=chasi_sha256_file(fullfile(root,'02_研究背景搭建','matlab','chasi_config.m'));
parent=fileparts(inputPath); if ~isfolder(parent), mkdir(parent); end
tmp=[inputPath '.partial']; assert(~isfile(tmp));
save(tmp,'data','-v7.3'); assert(movefile(tmp,inputPath));
hash=chasi_sha256_file(inputPath);
fid=fopen([inputPath '.sha256'],'wt'); assert(fid>=0);
fprintf(fid,'%s\n',hash); fclose(fid);
info=struct('inputPath',inputPath,'inputSha256',hash);
end

function [est,truth,cov,source]=advance(est,truth,cov,source,R,ids,dt,ttiMs,speed,fc,stream)
if dt==0,return;end
assert(dt>0,'Future observation would be used');
a=besselj(0,2*pi*(speed/3.6)*(fc*1e9)/299792458*(dt*ttiMs/1000));
[K,L,Nr,Nt]=size(est);
for k=1:K,for l=1:L
 prior=R(:,:,ids(k),l);noise=chol(prior,'lower')*cn(stream,Nr,Nt);
 old=reshape(truth(k,l,:,:),Nr,Nt);
 truth(k,l,:,:)=reshape(a*old+sqrt(1-a^2)*noise,1,1,Nr,Nt);
 est(k,l,:,:)=a*est(k,l,:,:);
 for q=1:4
  Q=reshape(source(k,l,:,:,q),Nr,Nr)*a^2;
  if q==4,Q=Q+(1-a^2)*prior;end
  source(k,l,:,:,q)=reshape(Q,1,1,Nr,Nr,1);
 end
 cov(k,l,:,:)=reshape(sum(source(k,l,:,:,:),5),1,1,Nr,Nr);
end,end
end

function x=cn(stream,varargin)
x=(randn(stream,varargin{:})+1i*randn(stream,varargin{:}))/sqrt(2);
end
