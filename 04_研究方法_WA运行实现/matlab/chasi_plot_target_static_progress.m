function chasi_plot_target_static_progress(checkpointDir,figureDir)
%CHASI_PLOT_TARGET_STATIC_PROGRESS Plot actual completed outer history only.
% Unequal-length curves are explicitly progress diagnostics, not a fair
% same-iteration performance comparison.
ids={'ideal_sumse','pilot_sumse','ideal_pf','pilot_pf'};
names={'理想·总SE','导频·总SE','理想·PF','导频·PF'};
colors=[0.00 0.38 0.68;0.82 0.34 0.12; ...
    0.10 0.55 0.31;0.58 0.27 0.62];
states=cell(size(ids));
for i=1:numel(ids)
    path=fullfile(checkpointDir,[ids{i} '.mat']);
    if isfile(path)
        saved=load(path,'state');
        if saved.state.iter>0, states{i}=saved.state; end
    end
end
assert(any(~cellfun(@isempty,states)),'No completed outer step to plot');
if ~isfolder(figureDir), mkdir(figureDir); end
font='Microsoft YaHei';
fig=figure('Visible','off','Color','white','Position',[100 100 1300 760]);
tile=tiledlayout(fig,2,2,'TileSpacing','compact');
metrics={'trueSumSE','truePF','relativeObjectiveChange','precoderChange'};
titles={'共同 Htrue：总 SE','共同 Htrue：PF 效用', ...
    '外层相对目标变化（阈值 10^{-4}）', ...
    '预编码相对变化（阈值 5×10^{-3}）'};
for q=1:4
    ax=nexttile(tile); hold(ax,'on'); grid(ax,'on');
    for i=1:numel(ids)
        if isempty(states{i}), continue; end
        h=states{i}.history;
        y=arrayfun(@(row)row.(metrics{q}),h);
        plot(ax,1:numel(h),y,'o-','LineWidth',1.7,'Color',colors(i,:), ...
            'DisplayName',names{i});
    end
    if q==3, yline(ax,1e-4,'k--','HandleVisibility','off'); end
    if q==4, yline(ax,5e-3,'k--','HandleVisibility','off'); end
    if q>=3
        set(ax,'YScale','log');
        bounds=ylim(ax);
        if q==3, threshold=1e-4; else, threshold=5e-3; end
        ylim(ax,[min(bounds(1),threshold/2),max(bounds(2),threshold*2)]);
    end
    title(ax,titles{q},'FontName',font);
    xlabel(ax,'已完成外层轮次','FontName',font);
    legend(ax,'Location','best','FontName',font);
end
title(tile,['静态高斯导频项目假设 v1｜210 UE／17 组｜', ...
    '各分支当前实际长度；非公平同轮次比较'],'FontName',font);
exportgraphics(fig,fullfile(figureDir,'static_progress.png'),'Resolution',250);
exportgraphics(fig,fullfile(figureDir,'static_progress.pdf'),'ContentType','vector');
close(fig);
fig=figure('Visible','off','Color','white','Position',[100 100 1300 480]);
tile=tiledlayout(fig,1,2,'TileSpacing','compact');
fields={'admmPrimal','admmDual'};
limits={'epsPrimal','epsDual'};
labels={'ADMM 原始残差／阈值','ADMM 对偶残差／阈值'};
for q=1:2
    ax=nexttile(tile); hold(ax,'on'); grid(ax,'on');
    for i=1:numel(ids)
        if isempty(states{i}), continue; end
        h=states{i}.history;
        ratio=arrayfun(@(row)row.(fields{q})/row.(limits{q}),h);
        semilogy(ax,1:numel(h),max(ratio,realmin),'o-','LineWidth',1.7, ...
            'Color',colors(i,:),'DisplayName',names{i});
    end
    yline(ax,1,'k--','HandleVisibility','off');
    bounds=ylim(ax);
    ylim(ax,[min(bounds(1),0.1),max(bounds(2),2)]);
    xlabel(ax,'已完成外层轮次','FontName',font);
    ylabel(ax,labels{q},'FontName',font);
    legend(ax,'Location','best','FontName',font);
end
title(tile,'静态高斯导频项目假设 v1｜ADMM 末次内层残差；虚线为阈值 1', ...
    'FontName',font);
exportgraphics(fig,fullfile(figureDir,'static_admm.png'),'Resolution',250);
exportgraphics(fig,fullfile(figureDir,'static_admm.pdf'),'ContentType','vector');
close(fig);
end
