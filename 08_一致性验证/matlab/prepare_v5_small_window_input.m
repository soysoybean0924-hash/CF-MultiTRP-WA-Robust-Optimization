function path=prepare_v5_small_window_input(root)
% New reduced-cost identity for a 1002-TTI causal/ACK audit only.
startup;
source=fullfile(root,'03_研究场景构造','inputs', ...
    'srs_csi_joint_v5_prb50_attribution_smoke','input.mat');
out=fullfile(root,'08_一致性验证','v5_window_fixture_20261113');
if ~isfolder(out),mkdir(out);end
path=fullfile(out,'input.mat');
if isfile(path)
    assert(strcmp(chasi_sha256_file(path),strtrim(fileread([path '.sha256']))));
    return
end
x=load(source,'data');d=x.data;
r=d.resources(1);ids=double(r.ueIds(1:2));
d.priorRx=d.priorRx(:,:,ids,1:2);
d.pilotId=d.pilotId(ids,:);d.pilotEnergy=d.pilotEnergy(ids,:);
d.sceneConfig.users=2;d.sceneConfig.trps=2;
r.ueIds=int32((1:2)');r.home_trp=int32((1:2)');
r.H=r.H(1:2,1:2,:,1:4);r.H_true=r.H_true(1:2,1:2,:,1:4);
r.C=r.C(1:2,1:2,:,:);r.C_stat=r.C;
r.sourceCov5=r.sourceCov5(1:2,1:2,:,:,:);
r.sourceCov=r.sourceCov(1:2,1:2,:,:,:);
r.servingMask=true(2,2);r.prior_gain=r.prior_gain(1:2,1:2);
r.omitted_gain=zeros(2,2,'single');r.pilot_ids=d.pilotId;
r.co_pilot_indices={1;2};r.srs_diagnostic.srsSenderIds=(1:2)';
r.srs_diagnostic.scheduledIds=(1:2)';
r.resourceWeightNu=.5;r.tx_rb_count=136.5;
d.resources=r;d.prbSchedule=chasi_exact_half_prb_schedule(273,1002);
d.prbEffective=.5;d.fullScaleInput=false;
d.fixtureIdentity='v5_window_two_ue_two_trp_one_group_1002tti';
d.measurementSeed=d.measurementSeed+333;d.trafficSeed=d.trafficSeed+333;
data=d;save(path,'data','-v7.3');
fid=fopen([path '.sha256'],'w');fprintf(fid,'%s\n',chasi_sha256_file(path));fclose(fid);
end
