function plot_original_cdf()
% Plot only completed full-size WMMSE/equal-power CDF data.
root=fileparts(fileparts(mfilename('fullpath')));
pyfile=fullfile(root,'reproduction','python','original_cdf_2000.mat');
mlfile=fullfile(root,'reproduction','matlab','matlab_cdf_2000.mat');
assert(isfile(pyfile) && isfile(mlfile), ...
    'Both formal 2000-setup input MAT files are required.');
p=load(pyfile); m=load(mlfile);
assert(all(isfinite(p.SE_MR_WMMSE_ADMM(:))));
assert(all(isfinite(m.curves.MR_sumSE(:))));
outdir=fullfile(root,'reproduction','comparison');
cases={ ...
 'MR sum-SE','SE_MR_WMMSE_ADMM','MR_sumSE'; ...
 'RZF sum-SE','SE_RZF_WMMSE_ADMM','RZF_sumSE'; ...
 'MR PF','SE_MR_WMMSE_PF_ADMM','MR_PF'; ...
 'RZF PF','SE_RZF_WMMSE_PF_ADMM','RZF_PF'};
for i=1:size(cases,1)
    pyvalues=p.(cases{i,2}); mlvalues=m.curves.(cases{i,3});
    py=sort(pyvalues(:)); ml=sort(mlvalues(:));
    y=linspace(0,1,numel(py))'; % main_DNN.py Yvals definition
    f=figure('Visible','off','Color','w');
    plot(py,y,'LineWidth',1.5); hold on;
    plot(ml,y,'LineWidth',1.5);
    grid on; xlabel('SE per UE (bit/s/Hz)'); ylabel('CDF');
    title([cases{i,1} ' — independent random draws']);
    legend('Original Python','Independent MATLAB','Location','southeast');
    stem=strrep(cases{i,1},' ','_');
    savefig(f,fullfile(outdir,['cdf_' stem '.fig']));
    exportgraphics(f,fullfile(outdir,['cdf_' stem '.png']));
    close(f);
    writematrix([py ml y],fullfile(outdir,['cdf_' stem '.csv']));
    pySum=sort(sum(pyvalues,1)).';
    mlSum=sort(sum(mlvalues,1)).';
    ysum=linspace(0,1,numel(pySum))'; % original Yvals_sum
    f=figure('Visible','off','Color','w');
    plot(pySum,ysum,'LineWidth',1.5);hold on;
    plot(mlSum,ysum,'LineWidth',1.5);
    grid on; xlabel('Network sum SE (bit/s/Hz)'); ylabel('CDF');
    title([cases{i,1} ' network sum — independent draws']);
    legend('Original Python','Independent MATLAB','Location','southeast');
    savefig(f,fullfile(outdir,['sum_cdf_' stem '.fig']));
    exportgraphics(f,fullfile(outdir,['sum_cdf_' stem '.png']));
    close(f);
    writematrix([pySum mlSum ysum], ...
        fullfile(outdir,['sum_cdf_' stem '.csv']));
end
end
