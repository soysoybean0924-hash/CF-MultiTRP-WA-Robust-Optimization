function run_original_matlab_setup()
% Independent random draw, 1 setup x 100 realizations; diagnostic only.
root=fileparts(fileparts(mfilename('fullpath')));
tic; scene=generate_original_setup(4873256,100);
channelSeconds=toc;
stats=statistics_from_channels(scene.H,scene.Hhat,scene.p);
rng(4873256,'twister');
cases={'sumSE_MR','sumSE_RZF','PF_MR','PF_RZF'};
results=struct(); tic;
for i=1:numel(cases)
    name=cases{i};
    if contains(name,'RZF'), suffix='RZF'; else, suffix='MR'; end
    if startsWith(name,'PF'), mode='PF'; else, mode='sumSE'; end
    a=stats.(['a_' suffix]); B=stats.(['B_' suffix]);
    draws=rand(size(a,1),size(a,2),100);
    results.(name)=wmmse_admm_power(a,B,scene.Pmax,mode,draws,100,5000);
    fprintf('%s objective=%g outer=%d stop=%s maxAPpower=%g\n', ...
        name,results.(name).objective(end), ...
        numel(results.(name).objective)-1, ...
        results.(name).stopReason,results.(name).maxPower);
end
optimizerSeconds=toc;
save(fullfile(root,'reproduction','matlab','independent_original_setup.mat'), ...
    'scene','stats','results','channelSeconds','optimizerSeconds','-v7.3');
fprintf('channel_s=%g optimizer_s=%g\n',channelSeconds,optimizerSeconds);
end
