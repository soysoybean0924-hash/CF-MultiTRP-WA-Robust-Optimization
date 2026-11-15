function state=run_v3_policy_batch(root,inputPath,configPath,outDir,targetOuter,maxGroups)
% Resume-safe v3 policy runner. Default maxGroups=17; tests may use one.
if nargin<6 || isempty(maxGroups),maxGroups=17;end
startup;assert(targetOuter>=1 && maxGroups>=1 && maxGroups<=17);
assert(isfile(inputPath) && isfile(configPath));
inputHash=chasi_sha256_file(inputPath);
assert(strcmp(inputHash,strtrim(fileread([inputPath '.sha256']))));
x=load(inputPath,'data');data=x.data;
assert(ismember(data.modelId,{'srs_overlap_jakes_v3','srs_csi_joint_v4','srs_csi_joint_v5'}));
spec=jsondecode(fileread(configPath));mode=char(spec.mode);
assert(ismember(mode,{'nominal','physical','uniform','source'}));
policy=spec.policy;policy.enabled=logical(policy.enabled);
if strcmp(mode,'source')
    if strcmp(policy.strengthMode,'link_source')
        assert(isfield(policy,'calibration') && policy.calibration.gainErrorScale>0);
    end
else
    policy.enabled=false(4,3);
end
groups=1:maxGroups;[resources,audit]=chasi_v3_source_policy(data,groups,policy);
physicalC=cell(1,maxGroups);
for g=groups
    physicalC{g}=resources(g).C;
    switch mode
        case 'nominal'
            resources(g).C=zeros(size(resources(g).C),'like',resources(g).C);
            resources(g).policyCov=[];
        case 'physical'
            resources(g).policyCov=[];
        case 'uniform'
            assert(spec.uniformExtra>=0 && spec.uniformExtra<=1);
            resources(g).policyCov=spec.uniformExtra*double(resources(g).C);
    end
end
c=data.config;budget=min([resources.powerWattPerTrp]);
cfg=struct('rank',c.rank,'budget',budget,'outer_iterations',1, ...
    'outer_rel_tol',c.outerObjectiveRelativeTolerance, ...
    'outer_precoder_tol',c.outerCovarianceRelativeTolerance, ...
    'admm_iterations',c.admmIterations,'cg_iterations',c.cgIterations, ...
    'eps_abs',c.admmAbsoluteTolerance,'eps_rel',c.admmRelativeTolerance, ...
    'pf_floor',c.pfFloor,'users',data.sceneConfig.users);
if strcmp(mode,'nominal'),cfg.branch='matrix_pf_nominal';
else,cfg.branch='matrix_pf_statistical';end
identity=struct('inputSha256',inputHash,'configSha256', ...
    chasi_sha256_file(configPath),'policySha256', ...
    chasi_sha256_file(which('chasi_v3_source_policy')), ...
    'runnerSha256',chasi_sha256_file(which('run_v3_policy_batch')), ...
    'solverSha256',chasi_sha256_file(which('chasi_matrix_wmmse')), ...
    'linkPowerSha256',chasi_sha256_file(which('chasi_v3_link_power')), ...
    'trueRateSha256',chasi_sha256_file(which('chasi_v3_true_rates')), ...
    'maxGroups',maxGroups,'mode',mode,'modelId',data.modelId);
if ~isfolder(outDir),mkdir(outDir);end
checkpoint=fullfile(outDir,'latest.mat');
if isfile(checkpoint)
    y=load(checkpoint,'state');state=y.state;
    assert(isequaln(state.identity,identity),'Checkpoint identity mismatch');
    assert(state.iter==numel(state.history));
    if state.complete || state.iter>=targetOuter,return;end
else
    state=struct('identity',identity,'iter',0,'history',struct([]), ...
        'precoder',{cell(1,maxGroups)},'complete',false, ...
        'stopReason','initial','physicalC',{physicalC},'audit',audit, ...
        'sceneSeed',data.sceneSeed,'measurementSeed',data.measurementSeed, ...
        'initialSeed',data.initialSeed,'trafficSeed',data.trafficSeed);
    state.precoder=[];
    saveCheckpoint(checkpoint,state);
end
for step=state.iter+1:targetOuter
    sol=chasi_matrix_wmmse(resources,cfg,state.precoder);
    v=sol.precoder; power=sol.power;
    linkPower=[];
    if strcmp(mode,'source') && any(policy.enabled(:,3))
        risk=cell(1,maxGroups);
        for g=groups
            sh=audit(g).sourceShare;lam=audit(g).appliedStrength;
            risk{g}=zeros(size(sh,1),size(sh,2));
            for q=1:4
                risk{g}=risk{g}+lam(:,:,q,3).*sh(:,:,q);
            end
            risk{g}=min(1,max(0,risk{g}));
        end
        pc=chasi_v3_link_power(resources,v,budget*ones(size(sol.power)),risk, ...
            spec.powerControl);
        v=pc.precoder;power=pc.actualTrpPowerW;linkPower=pc.linkPowerW;
    end
    assert(all(power<=budget+1e-6));
    admm=sol.last_admm_history(end);
    row=struct('iter',step,'solverReason',sol.reason, ...
        'objective',sol.optimization_objective, ...
        'relativeObjectiveChange',sol.trajectory(end).relative_objective_change, ...
        'maxPowerViolationW',max(power-budget), ...
        'admmPrimal',admm.primal,'admmDual',admm.dual, ...
        'epsPrimal',admm.eps_primal,'epsDual',admm.eps_dual, ...
        'innerIterations',admm.inner, ...
        'admmConverged',sol.trajectory(end).previous_admm_converged);
    if isempty(state.history),state.history=row;
    else,state.history(end+1)=row;end
    state.iter=step;state.precoder=v;state.actualTrpPowerW=power;
    state.linkPowerW=linkPower;state.serviceMask={resources.servingMask};
    state.policyCov={resources.policyCov};
    state.truePerUeRate=chasi_v3_true_rates(resources,v,data.sceneConfig.users);
    state.trueRateAfterLinkPowerAvailable=true;
    state.complete=strcmp(sol.reason,'relative_objective_precoder_and_admm_tolerance');
    if state.complete,state.stopReason='solver_converged';
    elseif ismember(sol.reason,{'optimization_objective_decrease', ...
            'inexact_precoder_step_no_proxy_improvement'})
        state.complete=true;state.stopReason=['solver_failed_' sol.reason];
    elseif ~row.admmConverged
        state.stopReason='inner_inexact_budget';
    else,state.stopReason='budget_unconverged';end
    saveCheckpoint(checkpoint,state);
    if state.complete,break;end
end
end

function saveCheckpoint(path,state)
tmp=[path '.partial'];assert(~isfile(tmp));
save(tmp,'state','-v7.3');[ok,msg]=movefile(tmp,path,'f');assert(ok,msg);
end
