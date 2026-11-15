function pl=chasi_v6_tr38901_outdoor_pathloss(scenario,isLos,fc,d2,d3,hbs)
% TR 38.901 v17.1.0 7.4.1, outdoor UMi street canyon / UMa.
assert(ismember(scenario,{'UMi','UMa'}) && fc>0 && ...
    d2>=10 && d3>=d2 && hbs>1.5);
bp=4*(hbs-1)*.5*fc*1e9/299792458;
if strcmp(scenario,'UMi')
 los1=32.4+21*log10(d3)+20*log10(fc);
 los2=32.4+40*log10(d3)+20*log10(fc)- ...
     9.5*log10(bp^2+(hbs-1.5)^2);
 if d2<=bp,los=los1;else,los=los2;end
 nlos=max(los,35.3*log10(d3)+22.4+21.3*log10(fc));
else
 los1=28+22*log10(d3)+20*log10(fc);
 los2=28+40*log10(d3)+20*log10(fc)- ...
     9*log10(bp^2+(hbs-1.5)^2);
 if d2<=bp,los=los1;else,los=los2;end
 nlos=max(los,13.54+39.08*log10(d3)+20*log10(fc));
end
if isLos,pl=los;else,pl=nlos;end
end
