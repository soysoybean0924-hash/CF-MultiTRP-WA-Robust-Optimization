function report=chasi_check_problem(c)
% Concrete configuration invariants; research targets are not pass criteria.
assert(c.rank>=1 && c.powerWattPerTrp>0 && c.usersPerCell>=1);
assert(c.resourceGroups==c.revisitTti/c.hopIntervalTti);
assert(c.pilotSequenceResources==2*c.uePilotLabels);
assert(c.admmAbsoluteTolerance>0 && c.admmRelativeTolerance>0);
report=struct('rank',c.rank,'power_w_per_trp',c.powerWattPerTrp, ...
 'resource_groups',c.resourceGroups,'srs_revisit_tti',c.revisitTti, ...
 'checked',true);
end
