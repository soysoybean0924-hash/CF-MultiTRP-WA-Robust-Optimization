function plot_srs_joint_zh(inputPath,checkpointDir,outDir)
% Independent final-point CDF; missing branches are written to manifest.
assert(~isfolder(outDir),'Plot directory exists; choose a new output directory.');
mkdir(outDir);
goals={'sumse','pf'}; methods={'ideal','nominal','robust'};
labels={'理想 CSI','名义 CSI','误差感知'};
missing={}; data=struct();
for q=1:2
    fig=figure('Visible','off'); hold on;
    for m=1:3
        branch=[methods{m} '_' goals{q}];
        path=fullfile(checkpointDir,[branch '.mat']);
        if ~isfile(path), missing{end+1}=branch; continue; end %#ok<AGROW>
        result=evaluate_srs_joint(inputPath,path);
        x=sort(result.perUserTrueSE(:)); y=(1:numel(x))'/numel(x);
        stairs(x,y,'DisplayName',labels{m}); data.(branch)=result;
    end
    xlabel('每用户真实 SE (bit/s/Hz)'); ylabel('经验 CDF');
    ylim([0 1]); title([upper(goals{q}) ' 优化目标：真实信道评价']);
    legend('show'); grid on;
    base=fullfile(outDir,['true_se_cdf_' goals{q}]);
    saveas(fig,[base '.png']); saveas(fig,[base '.svg']);
    saveas(fig,[base '.pdf']); close(fig);
end
save(fullfile(outDir,'plot_data.mat'),'data','missing');
fid=fopen(fullfile(outDir,'missing_branches.json'),'wt');
fprintf(fid,'%s',jsonencode(missing)); fclose(fid);
end
