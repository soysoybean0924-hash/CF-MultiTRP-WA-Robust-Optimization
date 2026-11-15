function validate_first_step_arrays(inputRoot,outPath)
% Export the first WMMSE state and capped ADMM candidate for Python parity.
s=chasi_load_frozen_binary(inputRoot,1);r=s.resources(1);c=s.config;
cfg=struct('rank',c.rank,'budget',c.powerWattPerTrp,'outer_iterations',1, ...
 'admm_iterations',3,'cg_iterations',10, ...
 'eps_abs',c.admmAbsoluteTolerance,'eps_rel',c.admmRelativeTolerance, ...
 'branch','matrix_statistical','pf_floor',c.pfFloor, ...
 'users',21*c.usersPerCell,'capture_first_step',true);
sol=chasi_matrix_wmmse(r,cfg,{r.initial});p=sol.firstProbe;
A=p.A{1};W=p.Omega{1};rate=p.rate{1};
M=p.parts{1}.M;B=p.parts{1}.B;Q=p.parts{1}.D;
V=p.admmCandidate{1};
save(outPath,'A','W','rate','M','B','Q','V','-v7');
end
