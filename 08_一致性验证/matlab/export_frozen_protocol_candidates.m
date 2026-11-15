function summary=export_frozen_protocol_candidates(root,targetOuter,base)
%EXPORT_FROZEN_PROTOCOL_CANDIDATES Preserve quantile conventions and gaps.
if nargin<3
    base=fullfile(root,'07_实验结果与进展','local_runs', ...
        'srs_overlap_jakes_v3_main');
end
if nargin<2, targetOuter=20; end
if targetOuter==20
    names={'independent_pc','ideal_pf','nominal_pf','robust_pf'};
    suffix='_frozen_experience.mat'; outName='frozen_protocol_candidates';
else
    names={'ideal_pf','nominal_pf','robust_pf'};
    suffix=sprintf('_outer%d_frozen_experience.mat',targetOuter);
    outName=sprintf('frozen_protocol_candidates_outer%d',targetOuter);
end
n=numel(names);
summary=struct('branch',names,'outerSteps',num2cell(targetOuter*ones(1,n)), ...
    'r5Valid',num2cell(nan(1,n)), ...
    'r5MissingZero',num2cell(nan(1,n)), ...
    'tail5Valid',num2cell(nan(1,n)), ...
    'meanValid',num2cell(nan(1,n)), ...
    'validUeSeconds',num2cell(zeros(1,n)), ...
    'completedBursts',num2cell(zeros(1,n)), ...
    'uncompletedTerminalBursts',num2cell(zeros(1,n)), ...
    'runtimeSeconds',num2cell(nan(1,n)));
if targetOuter==20, summary(1).outerSteps=0; end
for i=1:numel(names)
    x=load(fullfile(base,[names{i} suffix]),'evaluation');
    e=x.evaluation; candidate=e.completed_burst_ue_second_candidate;
    assert(e.online_wmmse_updates==0 && e.users_with_valid_samples==210);
    c=candidate.cdf_candidates;
    summary(i).r5Valid=c.valid_only.linear_percentile_mbps;
    summary(i).r5MissingZero=c.missing_as_zero.linear_percentile_mbps;
    summary(i).tail5Valid=c.valid_only.lower_tail_mean_mbps;
    summary(i).meanValid=candidate.mean_mbps_valid_ue_seconds;
    summary(i).validUeSeconds=candidate.valid_ue_second_samples;
    summary(i).completedBursts=candidate.completed_bursts_counted;
    summary(i).uncompletedTerminalBursts= ...
        candidate.terminal_bursts_not_queue_completed;
    summary(i).runtimeSeconds=e.runtime_seconds;
end
fid=fopen(fullfile(base,[outName '.json']),'wt');
assert(fid>=0); fprintf(fid,'%s',jsonencode(summary,PrettyPrint=true)); fclose(fid);
writetable(struct2table(summary),fullfile(base,[outName '.csv']));
end
