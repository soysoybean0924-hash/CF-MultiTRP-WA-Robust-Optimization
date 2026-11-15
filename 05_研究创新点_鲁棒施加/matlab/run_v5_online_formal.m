function summary=run_v5_online_formal(root,waOuter,targetTti,runTag,solverOverride)
% New v5 only. Full 17-group, 1000-TTI run; not executed in smoke checks.
startup;
if nargin<3,targetTti=1000;end
if nargin<4,runTag=sprintf('outer%d',waOuter);end
if nargin<5,solverOverride=struct();end
assert(targetTti>=1 && targetTti<=1000 && ...
    ischar(runTag) && ~isempty(regexp(runTag,'^[A-Za-z0-9_]+$','once')));
config=fullfile(root,'03_研究场景构造','config','srs_csi_joint_v5_prb50.json');
input=fullfile(root,'03_研究场景构造','inputs', ...
    'srs_csi_joint_v5_prb50','input.mat');
if ~isfile(input),prepare_srs_csi_v5_input(root,config,input);end
x=load(input,'data');d=x.data;
assert(d.fullScaleInput && numel(d.resources)==17 && ...
    strcmp(d.modelId,'srs_csi_joint_v5') && d.prbEffective==.5);
cfgDir=fullfile(root,'05_研究创新点_鲁棒施加','config','v5_prb50');
if ~isfile(fullfile(cfgDir,'source_all.json'))
    chasi_v3_make_policy_configs(root,cfgDir,input,1:8);
end
out=fullfile(root,'07_实验结果与进展','local_runs', ...
    'srs_csi_joint_v5_prb50',runTag);
run_v5_online_six(root,input,cfgDir,fullfile(out,'online'), ...
    targetTti,waOuter,solverOverride);
summary=chasi_v5_protocol_finalize(root,input, ...
    fullfile(out,'online','latest.mat'),fullfile(out,'branches'));
end
