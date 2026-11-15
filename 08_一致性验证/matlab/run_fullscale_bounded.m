function report=run_fullscale_bounded(branches)
% Full 17-group/210-UE input; bounded solver diagnostics, never acceptance.
if nargin<1, branches={'A','B','C','D','E'};end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);startup();
inputRoot=fullfile(root,'03_研究场景构造','inputs','structured_main_v1');
s=chasi_load_frozen_binary(inputRoot);
r=s.resources; Kall=s.manifest.srsSenderCount;
assert(numel(r)==17 && Kall==210 && sum([r.rb_count])==273);
for g=1:17
 K=size(r(g).H,1);L=size(r(g).H,2);
 gain=reshape(sum(abs(r(g).H).^2,[3 4]),K,L);
 gain(~r(g).servingMask')=-inf;
 [~,r(g).home_trp]=max(gain,[],2);
 r(g).omitted_gain=zeros(K,L);
end
cfg=struct('rank',2,'budget',s.config.powerWattPerTrp,...
 'outer_iterations',1,'outer_rel_tol',s.config.outerObjectiveRelativeTolerance,...
 'outer_precoder_tol',s.config.outerCovarianceRelativeTolerance,...
 'admm_iterations',10,'cg_iterations',10,'eps_abs',s.config.admmAbsoluteTolerance,...
 'eps_rel',s.config.admmRelativeTolerance,'pf_floor',s.config.pfFloor,...
 'users',Kall,'rho',.05);
out=fullfile(root,'07_实验结果与进展','local_runs',...
 ['fullscale_bounded_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);
mkdir(out);
report=struct('classification','fullscale_bounded_unconverged_diagnostic',...
 'input_manifest_sha256',s.input_manifest_sha256,'matlab_version',version,...
 'scene_seed',s.config.sceneSeed,'measurement_seed',s.config.measurementSeed,...
 'outer_cap',cfg.outer_iterations,'admm_cap',cfg.admm_iterations,...
 'groups',17,'users',210,'output',out,'branches',struct());
for j=1:numel(branches)
 name=branches{j};rr=r;
 switch name
 case 'A',cfg.branch='matrix_pf_nominal';
 case 'B',cfg.branch='ideal_pf';
 case 'C',cfg.branch='matrix_pf_nominal';
 case 'D',cfg.branch='matrix_pf_statistical';
 case 'E'
  cfg.branch='matrix_pf_statistical';
  % Diagnostic source-weighted beam proxy; not a validated three-layer policy.
  for g=1:17
   p=double(rr(g).sourceCov);
   rr(g).policyCov=.01*(.1*p(:,:,:,:,1)+.5*p(:,:,:,:,2)+...
     .1*p(:,:,:,:,3)+.3*p(:,:,:,:,4));
  end
 otherwise,error('Unknown branch');
 end
 t=tic;
 if strcmp(name,'A'),sol=chasi_no_cooperation(rr,cfg);
 else,sol=chasi_matrix_wmmse(rr,cfg,{rr.initial});end
 row=struct('reason',sol.reason,'converged',sol.converged,...
  'runtime_s',toc(t),'mean_se',mean(sol.true_per_ue_rate),...
  'r5_se',prctile(sol.true_per_ue_rate,5),...
  'total_power_w',sum(sol.power),'max_power_violation_w',sol.max_power_violation);
 report.branches.(name)=row;
 save(fullfile(out,[name '.mat']),'row','sol','-v7.3');
 fid=fopen(fullfile(out,'summary.json'),'w');fwrite(fid,jsonencode(report,PrettyPrint=true));fclose(fid);
 fprintf('%s: R5=%.6g mean=%.6g reason=%s elapsed=%.1fs\n',...
  name,row.r5_se,row.mean_se,row.reason,row.runtime_s);
end
end
