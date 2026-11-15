function acceptance_layer_masks(root,outDir)
% Exact protected service masks and artificial policy covariance export.
assert(~isfolder(outDir));mkdir(outDir);
s=chasi_load_frozen_binary(fullfile(root,'03_研究场景构造','inputs','structured_main_v1'));
cfg=jsondecode(fileread(fullfile(root,'05_研究创新点_鲁棒施加','config','study.json')));
layer=chasi_layer1_calibrate(s.resources,cfg);
for mode={'uniform','variable','fixed','mapped'}
 name=mode{1};[r,~]=chasi_layer_protection(s.resources(1),layer,cfg,name,cfg.fixed_level);
 writebin(fullfile(outDir,[name '.mask.bin']),uint8(r.servingMask),'uint8');
 writebin(fullfile(outDir,[name '.policy.real.bin']),real(r.policyCov),'double');
 writebin(fullfile(outDir,[name '.policy.imag.bin']),imag(r.policyCov),'double');
end
end
function writebin(path,x,precision)
fid=fopen(path,'wb');assert(fid>=0);fwrite(fid,x,precision);fclose(fid);
end
