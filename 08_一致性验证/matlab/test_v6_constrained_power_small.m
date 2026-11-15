function test_v6_constrained_power_small(root)
% Source only until MATLAB licensing and an independent calibration exist.
addpath(fullfile(root,'05_研究创新点_鲁棒施加','matlab'));
r=struct('H',complex(ones(1,1,1,1)), ...
    'C',complex(zeros(1,1,1,1)), ...
    'ueIds',1,'resourceWeightNu',1);
base={complex(ones(1,1,1,1))};
source={sqrt(2)*complex(ones(1,1,1,1))};
opts=struct('maximumMeanLoss',.1,'powerBudgetW',2);
out=chasi_v6_constrained_air_power(r,source,base,1,1e6,opts,[]);
assert(out.feasible && out.finalTotalAirPowerW< ...
    out.initialTotalAirPowerW && out.finalProxyLoss<.1);
assert(~out.calibrationValidated && ~isempty(out.trajectory));
r.H_true=r.H;
failed=false;
try
 chasi_v6_estimated_experience_proxy(r,source,1,1e6);
catch
 failed=true;
end
assert(failed,'Causal power proxy accepted H_true');
end
