function plot_fullscale_bounded(aDir,beDir)
% Plots actual MATLAB outputs; all curves are explicitly unconverged diagnostics.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
names={'A','B','C','D','E'};
dirs={aDir,beDir,beDir,beDir,beDir};
f=figure('Visible','off','Color','w'); hold on;
for j=1:5
 d=load(fullfile(dirs{j},[names{j} '.mat']),'sol');
 x=sort(double(d.sol.true_per_ue_rate(:)));
 plot(x,(1:numel(x))'/numel(x),'LineWidth',1.5,'DisplayName',names{j});
end
xlabel('Static true-channel UE SE (bit/s/Hz)');ylabel('Empirical CDF');
title('17 groups, 210 UE: bounded, unconverged diagnostic');
legend('Location','southeast');grid on;xlim([0 2]);
exportgraphics(f,fullfile(beDir,'fullscale_bounded_cdf.png'),'Resolution',160);
close(f);
end
