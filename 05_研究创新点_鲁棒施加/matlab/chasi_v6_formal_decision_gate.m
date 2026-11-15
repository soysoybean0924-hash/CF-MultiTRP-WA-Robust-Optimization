function [edgeStatus,powerStatus,edge,power]= ...
    chasi_v6_formal_decision_gate(proxy,blockingReasons,waConverged)
% Apply evidence gates after the unchanged two-objective formula.
% A proxy or missing BLER/channel/common sample can never produce pass.
assert(isstruct(proxy) && iscell(blockingReasons) && ...
    islogical(waConverged) && isscalar(waConverged));
edge=[];power=[];
if ~waConverged
 edgeStatus='unjudged_wa_not_converged';
 powerStatus=edgeStatus;return
end
objectiveOnly={'nonpositive_ideal_r5_gain_denominator', ...
    'nonpositive_independent_mean_or_power_denominator'};
shared=blockingReasons(~ismember(blockingReasons,objectiveOnly));
if ~isempty(shared)
 edgeStatus='unjudged_shared_dependency';
 powerStatus=edgeStatus;return
end
edge=proxy;power=proxy;
edgeStatus=proxy.edgeStatus;
powerStatus=proxy.powerStatus;
end
