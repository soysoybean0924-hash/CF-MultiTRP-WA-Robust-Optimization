function chasi_plot_campaign_cost(progressCsv, figureDir)
%CHASI_PLOT_CAMPAIGN_COST Plot measured wall time of completed resume batches.
% Original checkpoint rounds are excluded. This is not solver-only CPU time.
t = readtable(progressCsv,'TextType','string','Delimiter',',');
assert(all(t.exitCode == 0),'Manifest contains an incomplete MATLAB batch.');
assert(all(t.status == "batch_limit" | t.status == "target_native_tolerance"));
ids = ["ideal_sumse","pilot_sumse","ideal_pf","pilot_pf"];
colors = [0.00 0.38 0.68; 0.82 0.34 0.12; ...
    0.10 0.55 0.31; 0.58 0.27 0.62];
assert(all(ismember(t.branch,ids)));
if ~isfolder(figureDir), mkdir(figureDir); end
fig = figure('Visible','off','Color','white','Position',[80 80 1200 470]);
layout = tiledlayout(fig,1,2,'TileSpacing','compact');
ax1 = nexttile(layout); hold(ax1,'on'); grid(ax1,'on');
ax2 = nexttile(layout); hold(ax2,'on'); grid(ax2,'on');
for b = 1:numel(ids)
    part = sortrows(t(t.branch == ids(b),:),'outer');
    assert(height(part)>0 && numel(unique(part.outer)) == height(part));
    minutes = part.elapsedSeconds/60;
    plot(ax1,part.outer,minutes,'o-','Color',colors(b,:), ...
        'LineWidth',1.6,'MarkerSize',4,'DisplayName',char(ids(b)));
    plot(ax2,part.outer,cumsum(minutes)/60,'o-','Color',colors(b,:), ...
        'LineWidth',1.6,'MarkerSize',4,'DisplayName',char(ids(b)));
end
xlabel(ax1,'Completed outer iteration');
ylabel(ax1,'MATLAB batch wall time (min)');
xlabel(ax2,'Completed outer iteration');
ylabel(ax2,'Cumulative recorded batch time per branch (h)');
legend(ax1,'Location','best','Interpreter','none');
title(layout,'Static pilot LMMSE v1 | seed 19 | serial resume; startup included, original rounds excluded');
savefig(fig,fullfile(figureDir,'batch_cost.fig'));
exportgraphics(fig,fullfile(figureDir,'batch_cost.png'),'Resolution',300);
exportgraphics(fig,fullfile(figureDir,'batch_cost.pdf'),'ContentType','vector');
close(fig);
end
