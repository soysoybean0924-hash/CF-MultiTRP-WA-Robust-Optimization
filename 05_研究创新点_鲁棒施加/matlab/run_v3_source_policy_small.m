function state=run_v3_source_policy_small(root,baseMode,policy,outPath)
%RUN_V3_SOURCE_POLICY_SMALL One resource group, one WA update, paired input.
% Smoke only. No claim of convergence, protocol experience, or validation.
assert(ismember(baseMode,{'nominal','physical','source'}));
startup;
inputPath=fullfile(root,'03_研究场景构造','inputs', ...
    'srs_overlap_jakes_v3_main','input.mat');
assert(isfile(inputPath));
inputSha=chasi_sha256_file(inputPath);
assert(strcmp(inputSha,strtrim(fileread([inputPath '.sha256']))));
x=load(inputPath,'data'); data=x.data;
assert(strcmp(data.modelId,'srs_overlap_jakes_v3'));
if ~strcmp(baseMode,'source')
    assert(~any(policy.enabled(:)),'Only source mode permits policy switches');
end
[resources,audit]=chasi_v3_source_policy(data,1,policy);
physicalC=resources.C;
if strcmp(baseMode,'nominal')
    resources.C=zeros(size(physicalC),'like',physicalC);
    resources.policyCov=[];
elseif strcmp(baseMode,'physical')
    resources.policyCov=[];
end
% The native solver has a scalar cfg.budget, so use the most conservative
% available TRP cap. The policy audit retains the uncompressed vector.
budget=min(resources.powerWattPerTrp);
c=data.config;
cfg=struct('rank',c.rank,'budget',budget,'outer_iterations',1, ...
    'outer_rel_tol',c.outerObjectiveRelativeTolerance, ...
    'outer_precoder_tol',c.outerCovarianceRelativeTolerance, ...
    'admm_iterations',c.admmIterations,'cg_iterations',c.cgIterations, ...
    'eps_abs',c.admmAbsoluteTolerance,'eps_rel',c.admmRelativeTolerance, ...
    'pf_floor',c.pfFloor,'users',data.sceneConfig.users);
if strcmp(baseMode,'nominal'),cfg.branch='matrix_pf_nominal';
else,cfg.branch='matrix_pf_statistical';end
out=chasi_matrix_wmmse(resources,cfg,[]);
assert(all(isfinite(out.power(:))) && all(out.power(:)<=budget+1e-6));
state=struct('status','one_outer_unconverged','baseMode',baseMode, ...
    'group',1,'inputSha256',inputSha,'policyCodeSha256', ...
    chasi_sha256_file(which('chasi_v3_source_policy')), ...
    'solverCodeSha256',chasi_sha256_file(which('chasi_matrix_wmmse')), ...
    'policy',policy,'audit',audit,'physicalC',physicalC, ...
    'powerCapW',budget,'actualPowerW',out.power, ...
    'serviceMask',resources.servingMask,'precoder',{out.precoder}, ...
    'truePerUeRate',out.true_per_ue_rate,'solverReason',out.reason, ...
    'trajectory',out.trajectory);
if nargin>=4 && ~isempty(outPath)
    folder=fileparts(outPath);if ~isfolder(folder),mkdir(folder);end
    assert(~isfile(outPath),'Output already exists');
    save(outPath,'state','-v7.3');
end
end
