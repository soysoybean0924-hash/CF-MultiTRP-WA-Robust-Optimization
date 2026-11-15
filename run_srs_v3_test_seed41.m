% Predeclared independent seed test; same model and design parameters.
startup;
root=pwd;
configPath=fullfile(root,'03_研究场景构造','config', ...
    'srs_overlap_jakes_v3_test_seed41.json');
inputDir=fullfile(root,'03_研究场景构造','inputs', ...
    'srs_overlap_jakes_v3_test_seed41');
if ~isfolder(inputDir), mkdir(inputDir); end
inputPath=fullfile(inputDir,'input.mat');
if ~isfile(inputPath)
    info=prepare_srs_joint_input(root,configPath,inputPath);
    fprintf('TEST_INPUT %s\n',info.inputSha256);
end
loadedInput=load(inputPath,'data'); resources=loadedInput.data.resources;
base=fullfile(root,'07_实验结果与进展','local_runs', ...
    'srs_overlap_jakes_v3_test_seed41');
if ~isfolder(base), mkdir(base); end
names={'ideal_pf','nominal_pf','robust_pf'}; targetOuter=20;
summary=struct('branch',names,'inputSha256',repmat({''},1,3), ...
    'outerSteps',num2cell(zeros(1,3)), ...
    'complete',num2cell(false(1,3)), ...
    'stopReason',repmat({''},1,3), ...
    'trueSumSE',num2cell(nan(1,3)), ...
    'trueR5',num2cell(nan(1,3)), ...
    'truePF',num2cell(nan(1,3)), ...
    'totalPower',num2cell(nan(1,3)), ...
    'maxPowerViolation',num2cell(nan(1,3)), ...
    'innerSolved',num2cell(false(1,3)), ...
    'innerCount',num2cell(zeros(1,3)), ...
    'relativeObjectiveChange',num2cell(nan(1,3)), ...
    'precoderChange',num2cell(nan(1,3)));
for i=1:3
    branchDir=fullfile(base,names{i});
    if ~isfolder(branchDir), mkdir(branchDir); end
    checkpoint=fullfile(branchDir,'latest.mat'); already=0;
    if isfile(checkpoint)
        x=load(checkpoint,'state'); already=x.state.iter;
    end
    if already<targetOuter
        state=run_srs_joint_batch(root,inputPath,names{i}, ...
            targetOuter-already,checkpoint);
    else
        x=load(checkpoint,'state'); state=x.state;
    end
    e=state.history(end);
    summary(i).inputSha256=state.inputSha256;
    summary(i).outerSteps=state.iter;
    summary(i).complete=state.complete;
    summary(i).stopReason=state.stopReason;
    summary(i).trueSumSE=e.trueSumSE;
    summary(i).trueR5=e.trueRate5;
    summary(i).truePF=e.truePF;
    power=zeros(size(resources(1).H,2),1);
    for g=1:numel(resources)
        power=power+double(resources(g).resourceWeightNu)* ...
            reshape(sum(abs(state.precoder{g}).^2,[2 3 4]),[],1);
    end
    summary(i).totalPower=sum(power);
    summary(i).maxPowerViolation=max(0,max(power)-state.nativeCfg.budget);
    summary(i).innerSolved=e.innerSolved;
    summary(i).innerCount=e.innerCount;
    summary(i).relativeObjectiveChange=e.relativeObjectiveChange;
    summary(i).precoderChange=e.precoderChange;
    fprintf('INDEPENDENT_TEST %s outer=%d complete=%d SE=%.7g R5=%.7g\n', ...
        names{i},state.iter,state.complete,e.trueSumSE,e.trueRate5);
end
assert(numel(unique({summary.inputSha256}))==1);
fid=fopen(fullfile(base,'paired_twenty_step_detailed_summary.json'),'wt'); assert(fid>=0);
fprintf(fid,'%s',jsonencode(summary,PrettyPrint=true)); fclose(fid);
