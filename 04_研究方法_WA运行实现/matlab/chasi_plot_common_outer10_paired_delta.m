function chasi_plot_common_outer10_paired_delta(csvPath,figureDir)
%CHASI_PLOT_COMMON_OUTER10_PAIRED_DELTA Editable paired-UE diagnostic.
% Input must be the independently audited 210-UE common-outer-10 CSV.
t=readtable(csvPath,'Delimiter',',');
assert(height(t)==210 && isequal(t.userId,(1:210)'));
nmse=t.sampleNmse;
delta=[t.idealMinusPilotSumSE,t.idealMinusPilotPF];
assert(all(isfinite(nmse)) && all(nmse>=0) && all(isfinite(delta(:))));
if ~isfolder(figureDir), mkdir(figureDir); end
labels={'sum-SE objective','PF objective'};
colors=[0.07 0.43 0.68;0.83 0.34 0.16];
fig=figure('Visible','off','Color','white','Position',[80 80 1300 530]);
tiles=tiledlayout(fig,1,2,'TileSpacing','compact');
ax=nexttile(tiles); hold(ax,'on'); grid(ax,'on');
for j=1:2
    x=sort(delta(:,j)); n=numel(x);
    stepX=zeros(2*n+1,1); stepY=stepX;
    stepX(1)=x(1);
    for i=1:n
        stepX(2*i)=x(i); stepY(2*i)=i/n;
        stepX(2*i+1)=x(min(i+1,n)); stepY(2*i+1)=i/n;
    end
    plot(ax,stepX,stepY,'Color',colors(j,:),'LineWidth',1.7, ...
        'DisplayName',sprintf('%s (210 UEs)',labels{j}));
end
xline(ax,0,'k--','HandleVisibility','off');
xlabel(ax,'Ideal minus pilot UE SE (bit/s/Hz)');
ylabel(ax,'Cumulative proportion'); ylim(ax,[0 1]);
legend(ax,'Location','southeast');
ax=nexttile(tiles); hold(ax,'on'); grid(ax,'on');
for j=1:2
    scatter(ax,nmse,delta(:,j),15,colors(j,:),'filled', ...
        'MarkerFaceAlpha',0.65,'DisplayName',labels{j});
end
yline(ax,0,'k--','HandleVisibility','off');
xlabel(ax,'Per UE sample NMSE (dimensionless)');
ylabel(ax,'Ideal minus pilot UE SE (bit/s/Hz)');
legend(ax,'Location','northeast');
title(tiles,'Static pilot LMMSE v1 | seed 19 | common Htrue | completed outer 10 (unconverged)', ...
    'Interpreter','none');
savefig(fig,fullfile(figureDir,'paired_user_delta.fig'));
exportgraphics(fig,fullfile(figureDir,'paired_user_delta.png'),'Resolution',300);
exportgraphics(fig,fullfile(figureDir,'paired_user_delta.pdf'),'ContentType','vector');
close(fig);
end
