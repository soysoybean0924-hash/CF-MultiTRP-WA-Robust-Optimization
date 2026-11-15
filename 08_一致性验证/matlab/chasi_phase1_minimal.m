function report=chasi_phase1_minimal(root,outDir)
% Reproducible derived 2-group, 4-UE paired MATLAB smoke experiment.
% One frozen H_true/H_hat/C source for every branch; no regenerated channel.
assert(~isfolder(outDir),'Output directory already exists');
addpath(fullfile(root,'03_研究场景构造','matlab'));
addpath(fullfile(root,'04_研究方法_WA运行实现','matlab'));
s=chasi_load_frozen_binary(fullfile(root,'03_研究场景构造','inputs','structured_main_v1'),2);
assert(numel(s.resources)==2);
r=s.resources;
for g=1:2
    r(g).H_true=r(g).H_true(1:2,:,:,:);
    r(g).H=r(g).H(1:2,:,:,:);
    r(g).H_hat=r(g).H;
    r(g).C=r(g).C(1:2,:,:,:);
    r(g).C_stat=r(g).C;
    r(g).sourceCov=r(g).sourceCov(1:2,:,:,:,:);
    r(g).servingMask=r(g).servingMask(:,1:2);
    r(g).initial=r(g).initial(:,1:2,:,:);
    r(g).ueIds=(2*g-1:2*g)';
    r(g).omitted_gain=zeros(2,size(r(g).H,2));
    score=reshape(sum(abs(r(g).H).^2,[3,4]),2,[]);
    score(~r(g).servingMask')=-inf;
    [~,r(g).home_trp]=max(score,[],2);
    assert(max(abs(sum(r(g).sourceCov,5)-r(g).C),[],'all')<1e-5);
end
cfg=s.config;
solver=struct('rank',cfg.rank,'budget',cfg.powerWattPerTrp,...
    'outer_iterations',2,'outer_rel_tol',cfg.outerObjectiveRelativeTolerance,...
    'outer_precoder_tol',cfg.outerCovarianceRelativeTolerance,...
    'admm_iterations',15,'cg_iterations',10,...
    'eps_abs',cfg.admmAbsoluteTolerance,'eps_rel',cfg.admmRelativeTolerance,...
    'pf_floor',cfg.pfFloor,'users',4);
scene=struct('resources',r,'active_prb_ratio',32/273,...
    'config',struct('users',4,'durationTti',1000,'seed',cfg.sceneSeed,...
    'trafficModel','target_busy_fraction','trafficLoad',0.30,...
    'scsHz',30000,'ttiMs',1,'rbTotal',273));
[busy,tail,terminal]=chasi_traffic(4,1000,cfg.sceneSeed,'target_busy_fraction',0.30);
traffic=struct('busy',busy,'tail',tail,'terminal',terminal);
mkdir(outDir);
save(fullfile(outDir,'input.mat'),'r','scene','traffic','solver','-v7.3');
names={'no_cooperation','ideal_pf','matrix_pf_nominal','matrix_pf_statistical'};
report=struct('protocol','derived_two_group_pf_smoke_v1',...
    'manifest_sha256',s.input_manifest_sha256,'scene_seed',cfg.sceneSeed,...
    'measurement_seed',cfg.measurementSeed,'initial_seed',cfg.initialSeed,...
    'rb_fraction',32/273,'formal_acceptance',false,'branches',struct());
initial={r.initial};
for j=1:numel(names)
    solver.branch=names{j};
    started=tic;
    if j==1, sol=chasi_no_cooperation(r,solver);
    else, sol=chasi_matrix_wmmse(r,solver,initial); end
    met=chasi_evaluate(scene,sol,traffic);
    assert(sol.max_power_violation<1e-6);
    report.branches.(names{j})=struct('reason',sol.reason,...
        'converged',sol.converged,'runtime_s',toc(started),...
        'pf_objective',sol.optimization_objective,...
        'r5_mbps',met.experience_bottom5_mbps_valid_ues,...
        'mean_mbps',met.experience_mean_mbps_valid_ues,...
        'total_power_w',met.total_radiated_power_watt,...
        'max_power_violation_w',sol.max_power_violation);
    save(fullfile(outDir,[names{j} '.mat']),'sol','met','-v7.3');
    fprintf('%s: R5 %.9g Mbps, mean %.9g Mbps, power %.9g W, %s\n',...
        names{j},met.experience_bottom5_mbps_valid_ues,...
        met.experience_mean_mbps_valid_ues,met.total_radiated_power_watt,sol.reason);
end
base=report.branches.no_cooperation.r5_mbps;
ideal=report.branches.ideal_pf.r5_mbps;
report.r5_gain_retention=struct('denominator_mbps',ideal-base,...
    'nominal',(report.branches.matrix_pf_nominal.r5_mbps-base)/(ideal-base),...
    'total_c',(report.branches.matrix_pf_statistical.r5_mbps-base)/(ideal-base));
fid=fopen(fullfile(outDir,'summary.json'),'w','n','UTF-8');assert(fid>0);
fwrite(fid,jsonencode(report,PrettyPrint=true));fclose(fid);
end
