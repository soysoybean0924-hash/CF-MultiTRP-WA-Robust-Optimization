function chasi_plot_night_diagnostics(checkpointDir, figureDir)
%CHASI_PLOT_NIGHT_DIAGNOSTICS Plot completed outer records of one pilot run.
% Progress may have unequal branch lengths. No missing iteration is filled.
ids = {'ideal_sumse','pilot_sumse','ideal_pf','pilot_pf'};
labels = {'Ideal CSI / sum-SE','Pilot LMMSE / sum-SE', ...
    'Ideal CSI / PF','Pilot LMMSE / PF'};
colors = [0.00 0.38 0.68; 0.82 0.34 0.12; ...
    0.10 0.55 0.31; 0.58 0.27 0.62];
states = cell(1,4);
for b = 1:4
    path = fullfile(checkpointDir,[ids{b} '.mat']);
    assert(isfile(path),'Missing checkpoint: %s',path);
    item = load(path,'state'); states{b} = item.state;
    assert(states{b}.iter == numel(states{b}.history));
    assert(states{b}.iter > 0);
end
for b = 2:4
    assert(strcmp(states{1}.modelId,states{b}.modelId));
    assert(strcmp(states{1}.inputSha256,states{b}.inputSha256));
    assert(strcmp(states{1}.configSha256,states{b}.configSha256));
end
if ~isfolder(figureDir), mkdir(figureDir); end
fig = figure('Visible','off','Color','white','Position',[80 80 1400 850]);
t = tiledlayout(fig,2,3,'TileSpacing','compact');
spec = { ...
    'trueSumSE','True sum-SE (bit/s/Hz)',false,0; ...
    'truePF','True PF utility (natural log)',false,0; ...
    'relativeObjectiveChange','Relative outer objective change',true,1e-4; ...
    'precoderChange','Relative precoder change',true,5e-3; ...
    'innerCount','ADMM inner iterations',false,0; ...
    'residualRatio','Last ADMM residual / threshold',true,1};
for p = 1:size(spec,1)
    ax = nexttile(t); hold(ax,'on'); grid(ax,'on');
    for b = 1:4
        h = states{b}.history;
        if strcmp(spec{p,1},'residualRatio')
            primal = arrayfun(@(x)x.admmPrimal/x.epsPrimal,h);
            dual = arrayfun(@(x)x.admmDual/x.epsDual,h);
            y = max(primal,dual);
        else
            y = arrayfun(@(x)x.(spec{p,1}),h);
        end
        if spec{p,3}, y = max(y,realmin); end
        plot(ax,1:numel(h),y,'o-','Color',colors(b,:), ...
            'LineWidth',1.5,'MarkerSize',4, ...
            'DisplayName',sprintf('%s (%d)',labels{b},states{b}.iter));
    end
    if spec{p,4}>0
        yline(ax,spec{p,4},'k--','Threshold', ...
            'HandleVisibility','off');
    end
    if spec{p,3}, set(ax,'YScale','log'); end
    xlabel(ax,'Completed outer iteration'); ylabel(ax,spec{p,2});
    xlim(ax,[0.5,max(cellfun(@(s)s.iter,states))+0.5]);
    if p==1 || p==2, legend(ax,'Location','best'); end
end
title(t,sprintf('%s | seed 19 | progress, not convergence',states{1}.modelId), ...
    'Interpreter','none');
savefig(fig,fullfile(figureDir,'pilot_progress.fig'));
exportgraphics(fig,fullfile(figureDir,'pilot_progress.png'),'Resolution',250);
exportgraphics(fig,fullfile(figureDir,'pilot_progress.pdf'),'ContentType','vector');
close(fig);

fig = figure('Visible','off','Color','white','Position',[80 80 1200 460]);
t = tiledlayout(fig,1,2,'TileSpacing','compact');
for p = 1:2
    ax = nexttile(t); hold(ax,'on'); grid(ax,'on');
    for b = (2*p-1):(2*p)
        h = states{b}.history;
        y = arrayfun(@(x)x.estimatedObjective,h);
        plot(ax,1:numel(h),y,'o-','Color',colors(b,:), ...
            'LineWidth',1.7,'DisplayName',labels{b});
    end
    xlabel(ax,'Completed outer iteration');
    if p==1, ylabel(ax,'Design sum-SE (bit/s/Hz)');
    else, ylabel(ax,'Design PF utility (natural log)'); end
    legend(ax,'Location','best');
end
title(t,'Design-channel objectives; distinct from true metrics');
savefig(fig,fullfile(figureDir,'pilot_design_objective.fig'));
exportgraphics(fig,fullfile(figureDir,'pilot_design_objective.png'), ...
    'Resolution',250);
exportgraphics(fig,fullfile(figureDir,'pilot_design_objective.pdf'), ...
    'ContentType','vector');
close(fig);

file = fullfile(figureDir,'pilot_outer_records.csv');
fid = fopen(file,'w'); assert(fid>0);
cleaner = onCleanup(@()fclose(fid));
fprintf(fid,['branch,outer,trueSumSE_bit_s_Hz,truePF_natural_log,', ...
    'trueRate5_bit_s_Hz,designObjective,relativeObjectiveChange,', ...
    'precoderChange,innerCount,admmPrimal,epsPrimal,admmDual,epsDual,', ...
    'stopReason\n']);
for b = 1:4
    for j = 1:states{b}.iter
        h = states{b}.history(j);
        fprintf(fid,['%s,%d,%.17g,%.17g,%.17g,%.17g,%.17g,', ...
            '%.17g,%d,%.17g,%.17g,%.17g,%.17g,%s\n'], ...
            ids{b},h.outer,h.trueSumSE,h.truePF,h.trueRate5, ...
            h.estimatedObjective,h.relativeObjectiveChange, ...
            h.precoderChange,h.innerCount,h.admmPrimal,h.epsPrimal, ...
            h.admmDual,h.epsDual,states{b}.stopReason);
    end
end
end
