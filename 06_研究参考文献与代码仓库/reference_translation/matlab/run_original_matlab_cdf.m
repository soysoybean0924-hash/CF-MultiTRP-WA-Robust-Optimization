function run_original_matlab_cdf(target,batchSize)
% Independent MATLAB 2000-setup x 100-realization WMMSE CDF data runner.
% No Python calculation. Checkpointed. DNN curves not included.
if nargin==0, target=2000; end
if nargin<2, batchSize=100; end
assert(target>=1 && target<=2000 && target==floor(target));
assert(batchSize>=1 && batchSize==floor(batchSize));
root=fileparts(fileparts(mfilename('fullpath')));
outdir=fullfile(root,'reproduction','matlab');
checkpoint=fullfile(outdir,'matlab_cdf_checkpoint.mat');
names={'MR_equal','RZF_equal','MR_sumSE','RZF_sumSE','MR_PF','RZF_PF'};
if isfile(checkpoint)
    s=load(checkpoint);
    completed=s.completed; curves=s.curves; rng(s.randomState);
    apFixed=s.apFixed; stopCounts=s.stopCounts; innerCap=s.innerCap;
else
    completed=0; curves=struct();
    for i=1:numel(names),curves.(names{i})=nan(20,2000);end
    stopCounts=nan(4,2000);innerCap=false(4,2000);
end
tic;
endAt=min(target,completed+batchSize);
for setup=completed+1:endAt
    if setup==1
        scene=generate_original_setup(4873256,100); apFixed=scene.ap;
    else
        scene=generate_original_setup([],100,apFixed);
    end
    stats=statistics_from_channels(scene.H,scene.Hhat,scene.p);
    L=16; K=20;
    equal=sqrt(scene.Pmax/K)*ones(K,L);
    curves.MR_equal(:,setup)=evaluate_original_se( ...
        stats.a_MR,stats.B_MR,equal,scene.Pmax,scene.prelog);
    curves.RZF_equal(:,setup)=evaluate_original_se( ...
        stats.a_RZF,stats.B_RZF,equal,scene.Pmax,scene.prelog);
    branches={'MR_sumSE','RZF_sumSE','MR_PF','RZF_PF'};
    for j=1:4
        name=branches{j};
        if contains(name,'RZF'), suffix='RZF'; else, suffix='MR';end
        if contains(name,'PF'),mode='PF';else,mode='sumSE';end
        draws=rand(L,K,100);
        result=wmmse_admm_power(stats.(['a_' suffix]), ...
            stats.(['B_' suffix]),scene.Pmax,mode,draws,100,5000);
        stopCounts(j,setup)=numel(result.objective)-1;
        innerCap(j,setup)=any(result.innerCount==5000);
        curves.(name)(:,setup)=evaluate_original_se( ...
            stats.(['a_' suffix]),stats.(['B_' suffix]), ...
            result.mu,scene.Pmax,scene.prelog);
    end
    if mod(setup,25)==0 || setup==endAt
        completed=setup; randomState=rng; %#ok<NASGU>
        save(checkpoint,'completed','curves','randomState','apFixed', ...
            'stopCounts','innerCap');
        fprintf('completed=%d elapsed_s=%.2f\n',completed,toc);
    end
end
if target==2000 && completed==2000
    save(fullfile(outdir,'matlab_cdf_2000.mat'), ...
        'curves','stopCounts','innerCap');
    fprintf('FORMAL_MATLAB_CDF_DATA_COMPLETE=2000\n');
else
    fprintf('BATCH_COMPLETE=%d TARGET=%d\n',completed,target);
end
end
