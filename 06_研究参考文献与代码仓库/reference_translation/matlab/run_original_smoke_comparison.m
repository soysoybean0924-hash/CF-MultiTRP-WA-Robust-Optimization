function run_original_smoke_comparison()
root=fileparts(fileparts(mfilename('fullpath')));
s=load(fullfile(root,'reproduction','inputs','original_smoke_nonideal_100.mat'));
cases={'sumSE_MR','sumSE_RZF','PF_MR','PF_RZF'};
for index=1:numel(cases)
    name=cases{index};
    if contains(name,'RZF'), suffix='RZF'; else, suffix='MR'; end
    if startsWith(name,'PF'), mode='PF'; else, mode='sumSE'; end
    draws=s.(['draws_' name]);
    a=s.(['a_' suffix]); B=s.(['B_' suffix]);
    result=wmmse_admm_power(a,B,s.Pmax,mode,draws,size(draws,3),5000);
    python=s.(['mu_' suffix '_WMMSE_' ...
        choose(startsWith(name,'PF'),'PF_ADMM','ADMM')]);
    absolute=max(abs(result.mu(:)-python(:)));
    relative=absolute/max(1,max(abs(python(:))));
    fprintf('%s mu_max_abs=%.15g mu_rel_scaled=%.15g ',name,absolute,relative);
    fprintf('outer=%d/%d inner_cap=%d power=%.15g stop=%s\n', ...
        numel(result.objective)-1,size(draws,3), ...
        any(result.innerCount==5000),result.maxPower,result.stopReason);
    save(fullfile(root,'reproduction','matlab',[name '_comparison.mat']),'result');
end
end

function output=choose(condition,a,b)
if condition, output=a; else, output=b; end
end
