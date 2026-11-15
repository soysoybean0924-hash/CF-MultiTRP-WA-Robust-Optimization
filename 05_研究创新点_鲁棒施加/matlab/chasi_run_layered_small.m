function result=chasi_run_layered_small(inputRoot,configPath,mode,maxGroups,maxOuter,outDir,ablation,trigger,powerCap,level)
% Native MATLAB bounded campaign with immutable identity and per-step checkpoints.
if nargin<4 || isempty(maxGroups),maxGroups=2;end
if nargin<5 || isempty(maxOuter),maxOuter=2;end
if nargin<6,outDir='';end
if nargin<7 || isempty(ablation),ablation='none';end
if nargin<8 || isempty(trigger),trigger='proposed';end
assert(maxGroups>=1 && maxGroups<=17 && maxOuter>=1);
full=chasi_load_frozen_binary(inputRoot);cfg=jsondecode(fileread(configPath));
assert(maxOuter<=cfg.max_outer);
if nargin<10 || isempty(level),level=cfg.fixed_level;end
layer=chasi_layer1_calibrate(full.resources,cfg);
resources=full.resources(1:maxGroups);
if nargin>=9 && ~isempty(powerCap)
 assert(powerCap>0 && powerCap<=full.config.powerWattPerTrp);
 for g=1:numel(resources),resources(g).powerWattPerTrp=powerCap;end
end
[resources,mechanism]=chasi_layer_protection(resources,layer,cfg,mode, ...
 level,ablation,trigger);
c=full.config;
scfg=struct('rank',c.rank,'budget',mechanism.effective_cap_w, ...
 'outer_iterations',1,'outer_rel_tol',c.outerObjectiveRelativeTolerance, ...
 'outer_precoder_tol',c.outerCovarianceRelativeTolerance, ...
 'admm_iterations',c.admmIterations,'cg_iterations',c.cgIterations, ...
 'eps_abs',c.admmAbsoluteTolerance,'eps_rel',c.admmRelativeTolerance, ...
 'pf_floor',c.pfFloor,'users',21*c.usersPerCell,'branch','matrix_statistical');
identity=struct('input_manifest_sha256',full.input_manifest_sha256, ...
 'config_sha256',chasi_sha256_file(configPath), ...
 'solver_sha256',chasi_sha256_file(which('chasi_matrix_wmmse')), ...
 'runner_sha256',chasi_sha256_file([mfilename('fullpath') '.m']), ...
 'mode',mode,'level',level,'ablation',ablation,'trigger',trigger, ...
 'maxGroups',maxGroups,'powerCap',[]);
if nargin>=9,identity.powerCap=powerCap;end
v=[];history=struct([]);consecutive=0;startOuter=1;
if ~isempty(outDir)
 if ~isfolder(outDir),mkdir(outDir);end
 files=dir(fullfile(outDir,'checkpoint_*.mat'));
 if ~isempty(files)
  [~,idx]=sort({files.name});x=load(fullfile(outDir,files(idx(end)).name),'state');
  state=x.state;assert(isequaln(state.identity,identity),'checkpoint identity mismatch');
  v=state.precoder;history=state.history;consecutive=state.consecutive;
  startOuter=state.iterations+1;
  % Python layered runner treats convergence, failure, and inner inexactness
  % as terminal. Reopening that checkpoint must not append another update.
  if ~strcmp(state.status,'limit_unconverged')
   result=state;return
  end
 end
end
status='limit_unconverged';sol=[];
for outer=startOuter:maxOuter
 previous=v;sol=chasi_matrix_wmmse(resources,scfg,v);v=sol.precoder;
 inner=false;
 if numel(sol.trajectory)>=2 && isfield(sol.trajectory(1),'admm_converged')
  inner=logical(sol.trajectory(1).admm_converged);
 end
 rel=sol.trajectory(end).relative_objective_change;
 if isempty(previous),move=inf;else,move=chasi_covariance_change(previous,v);end
 meets=inner && isfinite(rel) && rel<=scfg.outer_rel_tol && ...
  move<=scfg.outer_precoder_tol;
 if meets,consecutive=consecutive+1;else,consecutive=0;end
 row=struct('outer',outer,'objective',sol.optimization_objective, ...
  'reason',sol.reason,'inner_converged',inner,'relative_objective_change',rel, ...
  'covariance_change',move,'true_weighted_sum_se',sum(sol.true_per_ue_rate), ...
  'power',sol.power);
 if isempty(history),history=row;else,history(end+1)=row;end %#ok<AGROW>
 if strcmp(sol.reason,'optimization_objective_decrease') || ...
    strcmp(sol.reason,'inexact_precoder_step_no_proxy_improvement')
  status='failed';
 elseif ~inner,status='inner_inexact';
 elseif consecutive>=cfg.consecutive_outer,status='converged';
 end
 state=struct('identity',identity,'mechanism',mechanism,'status',status, ...
  'iterations',outer,'consecutive',consecutive,'history',history, ...
  'precoder',{v},'true_per_ue_rate',sol.true_per_ue_rate, ...
  'true_weighted_sum_se',sum(sol.true_per_ue_rate), ...
  'mean_user_se',mean(sol.true_per_ue_rate), ...
  'input_manifest_sha256',full.input_manifest_sha256);
 if ~isempty(outDir)
  path=fullfile(outDir,sprintf('checkpoint_%04d.mat',outer));
  assert(~isfile(path));save(path,'state','-v7.3');
 end
 if ~strcmp(status,'limit_unconverged'),break;end
end
if isempty(sol)
 result=state;
else
 result=state;
end
end
