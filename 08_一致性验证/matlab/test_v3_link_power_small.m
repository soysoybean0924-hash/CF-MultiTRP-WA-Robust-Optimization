function out=test_v3_link_power_small()
r=struct('H',ones(2,2,1,2),'resourceWeightNu',1);
v={ones(2,2,2,1)/sqrt(2)};risk={zeros(2,2)};
out=chasi_v3_link_power(r,v,[1;1],risk, ...
    struct('powerPrice',2,'iterations',20,'step',.1));
assert(all(out.actualTrpPowerW<=1+1e-9));
assert(any(out.linkPowerW(:)<out.initialLinkPowerW(:)));
assert(size(out.linkPowerW,1)==2 && size(out.linkPowerW,2)==2 && ...
    size(out.linkPowerW,3)==1);
fprintf('LINK_POWER max=%.6g min=%.6g\n', ...
    max(out.actualTrpPowerW),min(out.linkPowerW(:)));
end
