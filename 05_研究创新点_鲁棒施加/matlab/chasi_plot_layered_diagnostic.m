function summary=chasi_plot_layered_diagnostic(result,outDir)
% Diagnostic only; CDF x-axis is individual user's resource-weighted SE.
assert(isfield(result,'true_per_ue_rate') && isfield(result,'input_manifest_sha256'));
if ~isfolder(outDir),mkdir(outDir);end
ue=double(result.true_per_ue_rate(:)); assert(all(isfinite(ue)) && all(ue>=0));
sorted=sort(ue);cdf=(1:numel(sorted))'/numel(sorted);
fig=figure('Visible','off');plot(sorted,cdf,'LineWidth',1.5);grid on;
xlabel('Per-user resource-weighted true SE (bit/s/Hz)');ylabel('Empirical CDF');
title('Diagnostic; one frozen input');
exportgraphics(fig,fullfile(outDir,'per_user_se_cdf.png'));close(fig);
summary=struct('input_manifest_sha256',result.input_manifest_sha256, ...
 'formal_result',false,'true_weighted_sum_se',sum(ue), ...
 'mean_user_se',mean(ue),'users',numel(ue), ...
 'metric','per-user true SE over resource weights; not protocol rate');
fid=fopen(fullfile(outDir,'STATUS.json'),'wt');assert(fid>=0);
fprintf(fid,'%s',jsonencode(summary,PrettyPrint=true));fclose(fid);
end
