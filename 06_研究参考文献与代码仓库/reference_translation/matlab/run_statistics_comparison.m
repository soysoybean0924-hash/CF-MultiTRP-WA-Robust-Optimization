function run_statistics_comparison()
root=fileparts(fileparts(mfilename('fullpath')));
s=load(fullfile(root,'reproduction','inputs','original_smoke_nonideal_100.mat'));
if isfield(s,'Np')
    recreated=pilot_mmse_from_channels(s.H,s.CorrR,s.Np,s.pilotIndex,100,10);
    fprintf('Hhat max_abs=%.15g\n',max(abs(recreated(:)-s.Hhat(:))));
end
t=statistics_from_channels(s.H,s.Hhat,100);
fields={'a_MR','B_MR','a_RZF','B_RZF'};
for i=1:numel(fields)
    name=fields{i};
    x=s.(name); y=t.(name);
    absolute=max(abs(x(:)-y(:)));
    relative=absolute/max(1,max(abs(x(:))));
    fprintf('%s max_abs=%.15g max_rel_scaled=%.15g\n',name,absolute,relative);
end
save(fullfile(root,'reproduction','matlab','statistics_comparison.mat'),'t');
end
