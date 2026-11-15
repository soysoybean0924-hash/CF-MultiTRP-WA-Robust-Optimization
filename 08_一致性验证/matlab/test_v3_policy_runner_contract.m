function test_v3_policy_runner_contract(root)
% A one-group, one-outer call is intentionally a smoke test, not convergence.
if nargin<1,root=fileparts(fileparts(fileparts(mfilename('fullpath'))));end
startup;
out=fullfile(root,'08_一致性验证','v3_policy_smoke_20261113.mat');
if isfile(out),delete(out);end
p=struct('enabled',false(4,3),'strength',ones(4,3), ...
    'powerReserveFraction',0.1,'serviceMargin',0);
s=run_v3_source_policy_small(root,'nominal',p,out);
assert(strcmp(s.status,'one_outer_unconverged') && s.group==1);
assert(isfile(out)); x=load(out,'state');
assert(isequaln(x.state.inputSha256,s.inputSha256));
assert(all(s.actualPowerW<=s.powerCapW+1e-6));
physical=run_v3_source_policy_small(root,'physical',p,'');
sourceOff=run_v3_source_policy_small(root,'source',p,'');
assert(norm(physical.truePerUeRate-sourceOff.truePerUeRate)<1e-9);
assert(norm(physical.actualPowerW-sourceOff.actualPowerW)<1e-9);
fprintf('V3_POLICY_RUNNER status=%s maxPowerViolation=%.6g\n', ...
    s.status,max(s.actualPowerW-s.powerCapW));
end
