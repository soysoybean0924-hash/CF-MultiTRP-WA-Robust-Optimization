% Same one-second frozen-beam diagnostic on predeclared independent seed.
startup;
addpath(fullfile(pwd,'08_一致性验证','matlab'));
outDir=fullfile(pwd,'07_实验结果与进展','local_runs', ...
    'srs_overlap_jakes_v3_test_seed41');
inputPath=fullfile(pwd,'03_研究场景构造','inputs', ...
    'srs_overlap_jakes_v3_test_seed41','input.mat');
run_srs_frozen_experience(pwd,outDir,20,inputPath);
if ~isfile(fullfile(outDir,'independent_pc.mat'))
    run_srs_independent_pc(pwd,fullfile(outDir,'independent_pc.mat'),inputPath);
end
if ~isfile(fullfile(outDir,'independent_pc_frozen_experience.mat'))
    run_srs_independent_experience(pwd,outDir,inputPath);
end
export_frozen_protocol_candidates(pwd,20,outDir);
