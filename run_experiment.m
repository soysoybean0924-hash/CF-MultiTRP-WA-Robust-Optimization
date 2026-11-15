function result=run_experiment(mode,configPath)
% Local research entry. Generate an independent frozen input, never overwrite.
if nargin<1 || isempty(mode),mode='basic';end
root=fileparts(mfilename('fullpath'));startup();
if nargin<2 || isempty(configPath)
 configPath=fullfile(root,'03_研究场景构造','config','srs_joint_snapshot_v2.json');
end
if strcmp(mode,'layered_diagnostic')
 inputRoot=fullfile(root,'03_研究场景构造','inputs','structured_main_v1');
 study=fullfile(root,'05_研究创新点_鲁棒施加','config','study.json');
 out=fullfile(root,'07_实验结果与进展','local_runs', ...
  char(datetime('now','Format','yyyyMMdd_HHmmss')));
 result=chasi_run_layered_small(inputRoot,study,'fixed',1,1,out);
 return
end
cfg=jsondecode(fileread(configPath));
out=fullfile(root,'07_实验结果与进展','local_runs',char(datetime('now','Format','yyyyMMdd_HHmmss')));
assert(~isfolder(out));mkdir(out);
inputPath=fullfile(out,'input.mat');
result=struct('mode',mode,'input',inputPath,'status','input_created');
prepare_srs_joint_input(root,configPath,inputPath);
if strcmp(mode,'input'),return;end
assert(strcmp(mode,'basic'),'Use basic or input; layered frozen adapter is not validated');
x=load(inputPath,'data');d=x.data;
solverCfg=struct('rank',cfg.rank,'budget',cfg.powerWattPerTrp, ...
 'outer_iterations',1,'admm_iterations',cfg.admmIterations, ...
 'cg_iterations',cfg.cgIterations,'eps_abs',cfg.admmAbsoluteTolerance, ...
 'eps_rel',cfg.admmRelativeTolerance,'pf_floor',cfg.pfFloor, ...
 'users',d.sceneConfig.users,'branch','matrix_statistical');
sol=chasi_matrix_wmmse(d.resources,solverCfg,d.initialPrecoder);
result.status=sol.reason;result.true_per_ue_rate=sol.true_per_ue_rate;
result.power=sol.power;result.input_sha256=chasi_sha256_file(inputPath);
save(fullfile(out,'small_result.mat'),'result','-v7.3');
end

