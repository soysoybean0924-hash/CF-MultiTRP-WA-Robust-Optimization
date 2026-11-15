function chasi_plot_input_nmse(analysisDir,figureDir)
%CHASI_PLOT_INPUT_NMSE Visualize frozen input estimation quality only.
% This figure contains no optimizer result. Each user/link is counted once.
summary=jsondecode(fileread(fullfile(analysisDir,'summary.json')));
assert(strcmp(summary.modelId,'target_gaussian_static_pilot_lmmse_hypothesis_v1'));
assert(strcmp(summary.inputSha256, ...
    'a55c8c16646d16ecc1150304f5676ff5c857d63d9c4c1d8c34ccd041dfc601da'));
users=readtable(fullfile(analysisDir,'per_user.csv'),'Delimiter',',');
links=readtable(fullfile(analysisDir,'per_link.csv'),'Delimiter',',');
assert(height(users)==210 && height(links)==4410);
assert(sum(links.serving)==488);
assert(all(isfinite(users.sampleNmse)) && all(isfinite(links.sampleNmse)));
if ~isfolder(figureDir), mkdir(figureDir); end
fig=figure('Visible','off','Color','white','Position',[80 80 1200 470]);
t=tiledlayout(fig,1,2,'TileSpacing','compact');
ax=nexttile(t); hold(ax,'on'); grid(ax,'on');
drawcdf(ax,users.sampleNmse,[0.00 0.38 0.68],'210 UE sample NMSE');
drawcdf(ax,users.theoryNmse,[0.82 0.34 0.12],'210 UE theory NMSE');
xlabel(ax,'Per UE NMSE (dimensionless)');
ylabel(ax,'Cumulative proportion'); ylim(ax,[0 1]);
legend(ax,'Location','southeast'); title(ax,'Per user across saved groups and 21 TRPs');
ax=nexttile(t); hold(ax,'on'); grid(ax,'on');
drawcdf(ax,links.sampleNmse,[0.00 0.38 0.68],'All 4410 UE-TRP links');
drawcdf(ax,links.sampleNmse(logical(links.serving)), ...
    [0.82 0.34 0.12],'488 ever-serving links');
xlabel(ax,'Per link sample NMSE (dimensionless)');
ylabel(ax,'Cumulative proportion'); ylim(ax,[0 1]);
legend(ax,'Location','southeast'); title(ax,'Per link subsets');
title(t,'Static pilot LMMSE v1 | seed 19 | frozen input; no solver result');
savefig(fig,fullfile(figureDir,'nmse_distribution.fig'));
exportgraphics(fig,fullfile(figureDir,'nmse_distribution.png'),'Resolution',300);
exportgraphics(fig,fullfile(figureDir,'nmse_distribution.pdf'),'ContentType','vector');
close(fig);
end

function drawcdf(ax,raw,color,label)
values=sort(max(raw(:),0)); n=numel(values);
xstep=zeros(2*n+1,1); ystep=zeros(2*n+1,1);
xstep(2:2:end)=values; xstep(3:2:end)=values;
ystep(2:2:end)=(0:n-1)'/n;
ystep(3:2:end)=(1:n)'/n;
plot(ax,xstep,ystep,'Color',color,'LineWidth',1.6,'DisplayName',label);
end
