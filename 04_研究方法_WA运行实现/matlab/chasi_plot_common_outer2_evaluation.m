function chasi_plot_common_outer2_evaluation(evaluationDir,figureDir,expectedOuter)
%CHASI_PLOT_COMMON_OUTER2_EVALUATION Draw four equal-round saved results.
% Each curve uses 210 users on one common true channel. It is not a
% cross-scenario network-sum CDF or a converged performance comparison.
summary=jsondecode(fileread(fullfile(evaluationDir,'summary.json')));
ids={'ideal_sumse','pilot_sumse','ideal_pf','pilot_pf'};
if nargin<3, expectedOuter=2; end
for b=1:4, assert(summary.(ids{b}).outer==expectedOuter); end
rates=readtable(fullfile(evaluationDir,'per_user_se.csv'), ...
    'TextType','string','Delimiter',',');
assert(height(rates)==210 && isequal(rates.userId,(1:210)'));
rateNames={'idealSumSE','pilotSumSE','idealPF','pilotPF'};
for b=1:4, assert(all(isfinite(rates.(rateNames{b})))); end
if ~isfolder(figureDir), mkdir(figureDir); end
colors=[0.00 0.38 0.68;0.82 0.34 0.12];

fig=figure('Visible','off','Color','white','Position',[80 80 1200 470]);
t=tiledlayout(fig,1,2,'TileSpacing','compact');
for p=1:2
    ax=nexttile(t); hold(ax,'on'); grid(ax,'on');
    for m=1:2
        b=2*p-2+m;
        values=sort(max(rates.(rateNames{b}),0));
        xstep=zeros(421,1); ystep=zeros(421,1);
        xstep(2:2:end)=values; xstep(3:2:end)=values;
        ystep(2:2:end)=(0:209)'/210;
        ystep(3:2:end)=(1:210)'/210;
        plot(ax,xstep,ystep, ...
            'Color',colors(m,:),'LineWidth',1.7, ...
            'DisplayName',strrep(ids{b},'_',' '));
    end
    if p==1, title(ax,'sum-SE design'); else, title(ax,'PF design'); end
    xlabel(ax,'True user SE (bit/s/Hz)');
    ylabel(ax,'Cumulative proportion of 210 UEs');
    ylim(ax,[0 1]); legend(ax,'Location','southeast');
end
title(t,sprintf('Static pilot LMMSE v1 | seed 19 | common Htrue | completed outer %d',expectedOuter));
saveall(fig,figureDir,'user_cdf');

fig=figure('Visible','off','Color','white','Position',[80 80 1200 760]);
t=tiledlayout(fig,2,2,'TileSpacing','compact');
for b=1:4
    ax=nexttile(t); hold(ax,'on'); grid(ax,'on');
    scatter(ax,rates.sampleNmse,rates.(rateNames{b}),14, ...
        'filled','MarkerFaceAlpha',0.65);
    xlabel(ax,'Per UE sample NMSE (dimensionless)');
    ylabel(ax,'True user SE (bit/s/Hz)');
    title(ax,strrep(ids{b},'_',' '));
end
title(t,sprintf('Static pilot LMMSE v1 | seed 19 | common completed outer %d',expectedOuter));
saveall(fig,figureDir,'user_nmse_vs_true_se');

fig=figure('Visible','off','Color','white','Position',[80 80 1200 760]);
t=tiledlayout(fig,2,2,'TileSpacing','compact');
for b=1:4
    power=readtable(fullfile(evaluationDir,[ids{b} '_trp_power.csv']), ...
        'Delimiter',',');
    assert(height(power)==21 && isequal(power.trpId,(1:21)'));
    assert(all(power.powerW<=power.budgetW+1e-7));
    ax=nexttile(t); hold(ax,'on'); grid(ax,'on');
    bar(ax,power.trpId,power.powerW,'FaceColor',[0.00 0.38 0.68]);
    yline(ax,20,'--','20 W budget','Color',[0.82 0.34 0.12]);
    xlabel(ax,'TRP index'); ylabel(ax,'Weighted transmit power (W)');
    ylim(ax,[0 22]); title(ax,strrep(ids{b},'_',' '));
end
title(t,sprintf('Static pilot LMMSE v1 | seed 19 | common completed outer %d | 0 violations',expectedOuter));
saveall(fig,figureDir,'trp_power');
end

function saveall(fig,folder,name)
savefig(fig,fullfile(folder,[name '.fig']));
exportgraphics(fig,fullfile(folder,[name '.png']),'Resolution',300);
exportgraphics(fig,fullfile(folder,[name '.pdf']),'ContentType','vector');
close(fig);
end
