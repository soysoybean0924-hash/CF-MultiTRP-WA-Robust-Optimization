% Formal-scale one-step paired diagnostic. Run from this repository root.
try
    startup;
    fprintf('RUN_START srs_overlap_jakes_v3_main %s\n', ...
        datestr(now,'yyyy-mm-dd HH:MM:SS'));
    result=run_srs_overlap_paired(pwd,20);
    outDir=fullfile(pwd,'07_实验结果与进展','local_runs', ...
        'srs_overlap_jakes_v3_main');
    save(fullfile(outDir,'paired_twenty_step.mat'),'result','-v7.3');
    for i=1:numel(result.states)
        s=result.states{i};
        fprintf('PAIRED_BRANCH %s iter=%d complete=%d reason=%s\n', ...
            result.branches{i},s.iter,s.complete,s.stopReason);
    end
catch ME
    fprintf(2,'%s\n',getReport(ME,'extended'));
    exit(1);
end
