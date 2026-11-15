function summaries=run_srs_frozen_experience(root,outDir,targetOuter,inputPath)
%RUN_SRS_FROZEN_EXPERIENCE Same-traffic frozen-beam 1-second diagnostic.
% This uses the existing mobile J0 evaluator with beams frozen at the
% snapshot. Its later pilot posterior updates do not enter decisions.
addpath(fullfile(root,'02_研究背景搭建','matlab'), ...
 fullfile(root,'03_研究场景构造','matlab'), ...
 fullfile(root,'04_研究方法_WA运行实现','matlab'));
if nargin<3, targetOuter=20; end
assert(targetOuter>=1 && targetOuter==floor(targetOuter));
if targetOuter==20, suffix='_frozen_experience.mat';
else, suffix=sprintf('_outer%d_frozen_experience.mat',targetOuter); end
if nargin<4
 inputPath=fullfile(root,'03_研究场景构造','inputs', ...
  'srs_overlap_jakes_v3_main','input.mat');
end
inputHash=chasi_sha256_file(inputPath);
if targetOuter==20, summaryName='frozen_experience_summary.json';
else, summaryName=sprintf('frozen_experience_outer%d_summary.json',targetOuter); end
summaryPath=fullfile(outDir,summaryName);
if isfile(summaryPath)
 summaries=jsondecode(fileread(summaryPath));
 assert(all(strcmp({summaries.inputSha256},inputHash)));
 return
end
s=load(inputPath,'data'); data=s.data;
cfg=data.sceneConfig;
cfg.trafficModel='target_busy_fraction'; cfg.trafficLoad=0.30;
scene=chasi_generate_scene(cfg);
scene.resources=data.resources; scene.config=cfg;
trafficSeed=data.trafficSeed;
[busy,tail,terminal]=chasi_traffic(cfg.users,1000,trafficSeed, ...
    cfg.trafficModel,cfg.trafficLoad);
traffic=struct('busy',busy,'tail',tail,'terminal',terminal);
options=struct('windowTti',1000,'offeredMbps',cfg.offeredMbps, ...
 'fastPeriodTti',cfg.fastPeriodTti,'slowPeriodTti',cfg.slowPeriodTti, ...
 'beamPolicy','frozen','trafficSeed',trafficSeed);
names={'ideal_pf','nominal_pf','robust_pf'};
summaries=struct('branch',names,'inputSha256',repmat({''},1,3), ...
 'outerSteps',num2cell(targetOuter*ones(1,3)), ...
 'experienceBottom5Mbps',num2cell(nan(1,3)), ...
 'experienceMeanMbps',num2cell(nan(1,3)), ...
 'completedBursts',num2cell(zeros(1,3)), ...
 'totalBursts',num2cell(zeros(1,3)), ...
 'maxPowerViolationWatt',num2cell(nan(1,3)), ...
 'runtimeSeconds',num2cell(nan(1,3)));
for i=1:3
 loaded=load(fullfile(outDir,names{i},'latest.mat'),'state'); state=loaded.state;
 assert(state.iter==targetOuter && strcmp(state.inputSha256,inputHash));
 outputPath=fullfile(outDir,[names{i} suffix]);
 if isfile(outputPath)
     prior=load(outputPath,'evaluation','inputHash');
     assert(strcmp(prior.inputHash,inputHash));
     evaluation=prior.evaluation;
 else
     solution=struct('precoder',{state.precoder},'branch',names{i});
     evaluation=chasi_time_evolution(scene,solution,options,traffic);
     save(outputPath,'evaluation','trafficSeed','inputHash','-v7.3');
 end
 summaries(i).inputSha256=inputHash;
 summaries(i).experienceBottom5Mbps=evaluation.experience_bottom5_mbps_valid;
 summaries(i).experienceMeanMbps=evaluation.experience_mean_mbps_valid;
 summaries(i).completedBursts=evaluation.completed_bursts;
 summaries(i).totalBursts=evaluation.total_bursts;
 summaries(i).maxPowerViolationWatt=evaluation.max_power_violation_watt;
 summaries(i).runtimeSeconds=evaluation.runtime_seconds;
 fprintf('EXPERIENCE %s R5=%.8g MEAN=%.8g completed=%d/%d time=%.3f\n', ...
  names{i},summaries(i).experienceBottom5Mbps, ...
  summaries(i).experienceMeanMbps, ...
  summaries(i).completedBursts,summaries(i).totalBursts, ...
  summaries(i).runtimeSeconds);
end
assert(~isfile(summaryPath),'Existing summary must be preserved');
fid=fopen(summaryPath,'wt'); assert(fid>=0);
fprintf(fid,'%s',jsonencode(summaries,PrettyPrint=true)); fclose(fid);
end
