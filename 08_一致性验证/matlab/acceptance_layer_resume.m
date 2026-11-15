function report=acceptance_layer_resume(root,outRoot)
% Compare uninterrupted two-outer run to disk-checkpoint 1+1 continuation.
assert(~isfolder(outRoot));mkdir(outRoot);
inp=fullfile(root,'03_研究场景构造','inputs','structured_main_v1');
cfg=fullfile(root,'05_研究创新点_鲁棒施加','config','study.json');
splitDir=fullfile(outRoot,'split');continuousDir=fullfile(outRoot,'continuous');
first=chasi_run_layered_small(inp,cfg,'fixed',1,1,splitDir);
resumed=chasi_run_layered_small(inp,cfg,'fixed',1,2,splitDir);
continuous=chasi_run_layered_small(inp,cfg,'fixed',1,2,continuousDir);
v1=resumed.precoder{1};v2=continuous.precoder{1};
report=struct('mode','fixed','input_manifest_sha256',resumed.identity.input_manifest_sha256, ...
 'first_status',first.status,'resumed_status',resumed.status, ...
 'continuous_status',continuous.status, ...
 'resumed_iterations',resumed.iterations,'continuous_iterations',continuous.iterations, ...
 'max_precoder_abs_diff',max(abs(v1-v2),[],'all'), ...
 'true_weighted_se_abs_diff',abs(resumed.true_weighted_sum_se-continuous.true_weighted_sum_se), ...
 'history_length_match',numel(resumed.history)==numel(continuous.history));
fid=fopen(fullfile(outRoot,'comparison.json'),'wt');assert(fid>=0);
fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));fclose(fid);
end
