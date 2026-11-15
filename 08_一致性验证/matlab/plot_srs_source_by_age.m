function plot_srs_source_by_age(root,outPath)
%PLOT_SRS_SOURCE_BY_AGE Internal analytic attribution, not field truth.
file=fullfile(root,'08_一致性验证','srs_v3_source_attribution.mat');
x=load(file,'report'); r=x.report;
inputPath=fullfile(root,'03_研究场景构造','inputs', ...
    'srs_overlap_jakes_v3_main','input.mat');
y=load(inputPath,'data'); resources=y.data.resources;
age=arrayfun(@(g)double(resources(g).age_ms),r.groups);
[sortedAge,order]=sort(age);
trace=r.meanTraceByGroup(order,:);
share=trace./sum(trace,2);
fig=figure('Visible','off','Color','w');
bar(sortedAge,share,'stacked','BarWidth',.9);
xlabel('SRS measurement age (ms)');
ylabel('Mean analytic posterior covariance trace share');
legend({'Target residual','Pilot pollution','Observation noise', ...
    'Aging innovation'},'Location','eastoutside');
ylim([0 1]); grid on;
title('Project Gaussian model attribution; not external source measurement');
exportgraphics(fig,outPath,'Resolution',180); close(fig);
end
