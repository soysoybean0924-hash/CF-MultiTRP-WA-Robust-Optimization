function result=run_srs_overlap_paired(root,targetOuter)
%RUN_SRS_OVERLAP_PAIRED Serial same-input reference/nominal/total-C WA.
% Checkpoints carry input/config/solver hashes and preserve nonconvergence.
arguments
 root (1,:) char
 targetOuter (1,1) double {mustBeInteger,mustBePositive}=1
end
inputPath=fullfile(root,'03_研究场景构造','inputs', ...
    'srs_overlap_jakes_v3_main','input.mat');
assert(isfile(inputPath),'Generate input with prepare_srs_joint_input first');
branches={'ideal_pf','nominal_pf','robust_pf'};
outRoot=fullfile(root,'07_实验结果与进展','local_runs', ...
    'srs_overlap_jakes_v3_main');
result=struct('inputSha256',chasi_sha256_file(inputPath), ...
    'branches',{branches},'states',{cell(1,numel(branches))});
for i=1:numel(branches)
    checkpoint=fullfile(outRoot,branches{i},'latest.mat');
    already=0;
    if isfile(checkpoint)
        saved=load(checkpoint,'state'); already=saved.state.iter;
    end
    if already<targetOuter
        result.states{i}=run_srs_joint_batch(root,inputPath, ...
            branches{i},targetOuter-already,checkpoint);
    else
        saved=load(checkpoint,'state'); result.states{i}=saved.state;
    end
end
end
