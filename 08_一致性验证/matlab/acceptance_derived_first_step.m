function acceptance_derived_first_step(root,outPath)
% Probe exact first U/W/quadratic/ADMM candidate on derived frozen slice.
s=chasi_load_frozen_binary(fullfile(root,'03_研究场景构造','inputs','structured_main_v1'),1);
r=s.resources(1);u=1:2;l=find(any(r.servingMask(:,u),2));
r.H=r.H(u,l,:,:);r.H_hat=r.H_hat(u,l,:,:);
r.H_true=r.H_true(u,l,:,:);r.C=r.C(u,l,:,:);r.C_stat=r.C_stat(u,l,:,:);
r.servingMask=r.servingMask(l,u);r.ueIds=r.ueIds(u);
initial=r.initial(l,u,:,:);r.initial=initial;c=s.config;
cfg=struct('rank',c.rank,'budget',c.powerWattPerTrp,'outer_iterations',1, ...
 'outer_rel_tol',c.outerObjectiveRelativeTolerance, ...
 'outer_precoder_tol',c.outerCovarianceRelativeTolerance, ...
 'admm_iterations',c.admmIterations,'cg_iterations',c.cgIterations, ...
 'eps_abs',c.admmAbsoluteTolerance,'eps_rel',c.admmRelativeTolerance, ...
 'pf_floor',c.pfFloor,'users',21*c.usersPerCell, ...
 'branch','matrix_statistical','capture_first_step',true);
sol=chasi_matrix_wmmse(r,cfg,{initial});p=sol.firstProbe;
A=p.A{1};W=p.Omega{1};rate=p.rate{1};
M=p.parts{1}.M;B=p.parts{1}.B;Q=p.parts{1}.D;
V=p.admmCandidate{1};admmHistory=p.admmHistory;
save(outPath,'A','W','rate','M','B','Q','V','admmHistory','-v7');
end
