function test_v6_formal_unjudged_gate(root)
% Full-window-like synthetic decisions; no input or BLER curve is created.
addpath(fullfile(root,'05_研究创新点_鲁棒施加','matlab'));
addpath(fullfile(root,'08_一致性验证','matlab'));
addpath(fullfile(root,'03_研究场景构造','matlab'));
prb=chasi_exact_half_prb_schedule([16*ones(1,16) 17],1000);
assert(isequal(prb.usedRbByTti(1:4),[136 137 136 137]) && ...
    prb.averageFraction==.5 && ...
    all(sum(prb.allocatedRbByGroupTti,1)==prb.usedRbByTti));
base=struct('rateMbps',[1;2],'valid',true(2,1), ...
    'totalAirPowerW',10);
ideal=base;ideal.rateMbps=[3;4];
method=base;method.rateMbps=[2.8;3.8];
method.totalAirPowerW=6;
proxy=chasi_two_objective_decision(base,ideal,method,struct());
assert(strcmp(proxy.edgeStatus,'pass') && ...
    strcmp(proxy.powerStatus,'pass'));
[e,p,edge,power]=chasi_v6_formal_decision_gate(proxy, ...
    {'missing_validated_sinr_mcs_bler_curve'},true);
assert(strcmp(e,'unjudged_shared_dependency') && ...
    strcmp(p,e) && isempty(edge) && isempty(power));
[e,p]=chasi_v6_formal_decision_gate(proxy,{},false);
assert(strcmp(e,'unjudged_wa_not_converged') && strcmp(p,e));
base.valid(2)=false;
missing=chasi_two_objective_decision(base,ideal,method,struct());
assert(strcmp(missing.edgeStatus,'missing_samples') && ...
    strcmp(missing.powerStatus,'missing_samples'));
[e,p]=chasi_v6_formal_decision_gate(missing, ...
    {'no_complete_common_ue_second_samples'},true);
assert(~strcmp(e,'pass') && ~strcmp(p,'pass'));
ideal.rateMbps=base.rateMbps;base.valid(:)=true;
noDen=chasi_two_objective_decision(base,ideal,method,struct());
assert(strcmp(noDen.edgeStatus,'nonpositive_denominator'));
end
