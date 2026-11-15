function plot_srs_overlap_paired(root,csvPath,outPath)
%PLOT_SRS_OVERLAP_PAIRED Descriptive static SE CDF, no convergence claim.
base=fullfile(root,'07_实验结果与进展','local_runs','srs_overlap_jakes_v3_main');
if nargin<2, csvPath=fullfile(base,'per_ue_static_se.csv'); end
if nargin<3, outPath=fullfile(base,'static_se_cdf.png'); end
t=readtable(csvPath);
value={t.idealSE,t.nominalSE,t.totalCSE};
names={'Ideal CSI','Nominal CSI','Total-C robust'};
fig=figure('Visible','off','Color','w'); hold on;
for i=1:3
 x=sort(value{i}); y=(1:numel(x))'/numel(x);
 stairs(x,y,'LineWidth',1.6,'DisplayName',names{i});
end
xlabel('Per-user static spectral efficiency (bit/s/Hz)');
ylabel('Empirical CDF'); grid on; legend('Location','southeast');
title('Same-input WA diagnostic; branches not yet converged');
exportgraphics(fig,outPath,'Resolution',180);
close(fig);
end
