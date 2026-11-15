function run_ideal_csi_comparison()
% Paired one-setup ideal-CSI diagnostic; original nonideal Hhat remains separate.
root=fileparts(fileparts(mfilename('fullpath')));
source=load(fullfile(root,'reproduction','inputs','original_smoke_nonideal_100.mat'),'H');
target=load(fullfile(root,'reproduction','inputs','original_smoke_ideal_100.mat'));
stats=statistics_from_channels(source.H,source.H,100);
fields={'a_MR','B_MR','a_RZF','B_RZF'};
for j=1:numel(fields)
    name=fields{j}; truth=target.(name); calculated=stats.(name);
    absolute=max(abs(truth(:)-calculated(:)));
    relative=absolute/max(1,max(abs(truth(:))));
    fprintf('ideal %s abs=%.15g rel_scaled=%.15g\n',name,absolute,relative);
    assert(relative<1e-9);
end
cases={'sumSE_MR','sumSE_RZF','PF_MR','PF_RZF'};
results=struct();
for j=1:numel(cases)
    label=cases{j};
    if contains(label,'RZF'), beam='RZF'; else, beam='MR'; end
    if startsWith(label,'PF'), mode='PF'; else, mode='sumSE'; end
    draws=target.(['draws_' label]);
    result=wmmse_admm_power(stats.(['a_' beam]),stats.(['B_' beam]), ...
        double(target.Pmax),mode,draws,size(draws,3),5000);
    python=target.(['mu_' label]);
    muAbsolute=max(abs(result.mu(:)-python(:)));
    muRelative=muAbsolute/max(1,max(abs(python(:))));
    se=evaluate_original_se(stats.(['a_' beam]),stats.(['B_' beam]), ...
        result.mu,double(target.Pmax),double(target.prelogFactor));
    pythonSE=target.(['SE_' label]);
    seAbsolute=max(abs(se(:)-pythonSE(:)));
    fprintf('ideal %s mu_abs=%.15g mu_rel_scaled=%.15g SE_abs=%.15g outer=%d/%d max_raw_AP_power=%.15g stop=%s\n', ...
        label,muAbsolute,muRelative,seAbsolute,numel(result.objective)-1, ...
        size(draws,3),result.maxPower,result.stopReason);
    assert(muRelative<1e-9 && seAbsolute<1e-9);
    results.(label)=result;
end
save(fullfile(root,'reproduction','matlab','ideal_csi_comparison.mat'), ...
    'results','-v7.3');
end
