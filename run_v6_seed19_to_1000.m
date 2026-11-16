function final=run_v6_seed19_to_1000(root,runTag,curvePath)
% Run the unmodified full-scale seed-19 v6 campaign in cumulative batches.
% A failed batch stops here; existing checkpoints, raw MAT files and logs
% are never promoted to an acceptance result by this wrapper.
if nargin<1 || isempty(root),root=fileparts(mfilename('fullpath'));end
if nargin<2 || isempty(runTag),runTag='seed19_full';end
if nargin<3,curvePath='';end
assert(isfolder(root) && ischar(runTag) && ...
    ~isempty(regexp(runTag,'^[A-Za-z0-9_]+$','once')));
assert(isempty(curvePath) || isfile(curvePath));
addpath(root);startup;
parentConfig=jsondecode(fileread(fullfile(root,'03_研究场景构造', ...
    'config','srs_csi_joint_v5_prb50.json')));
assert(parentConfig.sceneSeed==19 && parentConfig.usersPerCell==10 && ...
    parentConfig.prbUtilization==.5 && parentConfig.resourceGroups==17 && ...
    parentConfig.maxOuter==2000 && parentConfig.admmIterations==100 && ...
    parentConfig.cgIterations==100 && parentConfig.protocolTti==1000, ...
    'The requested seed-19 full-scale solver configuration is absent');
out=fullfile(root,'07_实验结果与进展','local_runs', ...
    'srs_csi_joint_v6_formal','source_all_strength1000',runTag);
if ~isfolder(out),mkdir(out);end
checkpoint=fullfile(out,'online','latest.mat');
completed=0;
if isfile(checkpoint)
    loaded=load(checkpoint,'state');
    completed=loaded.state.iter;
    assert(strcmp(loaded.state.identity.inputSha256, ...
        chasi_sha256_file(fullfile(root,'03_研究场景构造','inputs', ...
        'srs_csi_joint_v6_formal','input.mat'))), ...
        'Existing checkpoint belongs to a different frozen input');
    clear loaded
end
assert(completed>=0 && completed<=1000);
firstTarget=max(20,20*ceil((completed+1)/20));
if completed==1000,firstTarget=1000;end
progressPath=fullfile(out,'campaign_progress.log');
summary=[];
for target=firstTarget:20:1000
    stamp=datestr(now,'yyyymmdd_HHMMSS_FFF');
    logPath=fullfile(out,sprintf('batch_tti_%04d_%s.log',target,stamp));
    assert(~isfile(logPath),'Refusing to overwrite an existing batch log');
    started=tic;
    diary(logPath);
    try
        fprintf('FORMAL_BATCH_START target=%d previous=%d log=%s\n', ...
            target,completed,logPath);
        summary=run_v6_formal_batch(root,runTag,target,curvePath);
        assert(summary.ttiCount==target && isfile(checkpoint), ...
            'Batch checkpoint did not reach the requested cumulative TTI');
        checkpointHash=chasi_sha256_file(checkpoint);
        appendProgress(progressPath,sprintf( ...
            'tti=%d elapsed_s=%.3f checkpoint_sha256=%s status=complete', ...
            target,toc(started),checkpointHash));
        fprintf('FORMAL_BATCH_COMPLETE target=%d elapsed_s=%.3f\n', ...
            target,toc(started));
        diary off
        completed=target;
        if target<1000,clear summary;end
    catch ME
        appendProgress(progressPath,sprintf( ...
            'tti=%d elapsed_s=%.3f status=failed identifier=%s message=%s', ...
            target,toc(started),ME.identifier,ME.message));
        diary off
        rethrow(ME)
    end
end
assert(completed==1000 && ~isempty(summary), ...
    'The full 1000-TTI formal window has not completed');
assert(isfield(summary,'fourBranchPaired') && ...
    numel(summary.fourBranchPaired.raw)==4, ...
    'Final evaluation lacks four-branch paired raw results');
paired=summary.fourBranchPaired;
final=struct('inputSha256',summary.inputSha256, ...
    'checkpointSha256',chasi_sha256_file(checkpoint), ...
    'ttiCount',summary.ttiCount,'edgeStatus',summary.edgeStatus, ...
    'powerStatus',summary.powerStatus, ...
    'blockingReasons',{summary.blockingReasons}, ...
    'edgeRatio',paired.edgeRatio,'meanLoss',paired.meanLoss, ...
    'powerReduction',paired.powerReduction, ...
    'commonSampleCount',paired.commonSampleCount, ...
    'campaignProgressPath',progressPath);
finalPath=fullfile(out,'campaign_final_status.mat');
if isfile(finalPath)
    previous=load(finalPath,'final');
    assert(isequaln(previous.final,final), ...
        'Existing final campaign status differs');
    return
end
csvPath=fullfile(out,'campaign_four_branch_metrics.csv');
assert(~isfile(csvPath),'Refusing to overwrite branch metrics CSV');
csvPartial=[csvPath '.partial'];
assert(~isfile(csvPartial),'Inspect incomplete metrics CSV first');
fid=fopen(csvPartial,'w');assert(fid>0);
fprintf(fid,'branch,edge5_mbps,mean_mbps,mean_air_power_w,wa_converged,power_feasible,service_feasible\n');
for j=1:4
    raw=paired.raw(j);logs=raw.waLog;
    fprintf(fid,'%s,%.17g,%.17g,%.17g,%d,%d,%d\n', ...
        raw.name,raw.edge5Mbps,raw.meanMbps,raw.airPowerW, ...
        ~isempty(logs) && all([logs.converged]), ...
        ~isempty(logs) && all([logs.powerFeasible]), ...
        ~isempty(logs) && all([logs.serviceFeasible]));
end
fclose(fid);
[ok,msg]=movefile(csvPartial,csvPath);assert(ok,msg);
partial=[finalPath '.partial'];
assert(~isfile(partial),'Inspect incomplete campaign final status first');
save(partial,'final','-v7.3');
[ok,msg]=movefile(partial,finalPath);assert(ok,msg);
end

function appendProgress(path,line)
fid=fopen(path,'a');assert(fid>0,'Cannot append campaign progress');
fprintf(fid,'%s %s\n',datestr(now,31),line);
fclose(fid);
end
