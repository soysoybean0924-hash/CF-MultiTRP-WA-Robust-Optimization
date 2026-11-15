% Run after the three fixed-service branches reach the requested cap.
startup;
addpath(fullfile(pwd,'08_一致性验证','matlab'));
base=fullfile(pwd,'07_实验结果与进展','local_runs', ...
    'srs_overlap_jakes_v3_main');
validate_srs_source_components_small(fullfile(pwd,'08_一致性验证', ...
    'srs_source_components_small.mat'));
validate_srs_source_attribution(pwd,fullfile(pwd,'08_一致性验证', ...
    'srs_v3_source_attribution.mat'));
summary=summarize_srs_overlap_paired(pwd);
fid=fopen(fullfile(base,'paired_current_summary.json'),'wt'); assert(fid>=0);
fprintf(fid,'%s',jsonencode(summary,PrettyPrint=true)); fclose(fid);
export_srs_overlap_trajectory(pwd,fullfile(base,'paired_current_trajectory.csv'));
tail=export_srs_overlap_rates(pwd,fullfile(base,'per_ue_current_static_se.csv'));
fid=fopen(fullfile(base,'current_tail5_definition.json'),'wt'); assert(fid>=0);
fprintf(fid,'%s',jsonencode(tail,PrettyPrint=true)); fclose(fid);
export_srs_layer_activation(pwd,fullfile(base,'layer_activation_current.mat'));
export_frozen_protocol_candidates(pwd);
for i=1:numel(summary)
    fprintf('POSTCHECK %s outer=%d complete=%d reason=%s SE=%.7g R5=%.7g\n', ...
        summary(i).branch,summary(i).outerSteps,summary(i).complete, ...
        summary(i).stopReason,summary(i).trueSumSE,summary(i).trueR5);
end
