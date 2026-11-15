function p=chasi_v6_tr38901_lsp_table(scenario,isLos,carrierGhz,d2dMeters)
% TR 38.901 v17.1.0 Tables 7.5-6 and 7.5-7/8; outdoor UMi/UMa only.
assert(ismember(scenario,{'UMi','UMa'}) && islogical(isLos) && ...
    isscalar(isLos) && carrierGhz>0 && d2dMeters>=10);
f=carrierGhz;if strcmp(scenario,'UMa'),f=max(f,6);else,f=max(f,2);end
z=log10(1+f);q=log10(f);
if strcmp(scenario,'UMi')
 if isLos
  mu=[-0.24*z-7.14,-0.05*z+1.21,-0.08*z+1.73,-0.1*z+0.73];
  sd=[0.38,0.41,0.014*z+0.28,-0.04*z+0.34];
  p.clusters=12;p.rTau=3;p.xprDb=[9 3];p.clusterSpread=[5 3 17 7];
  p.lspDistanceMeters=[10 15 7 8 8 12 12];
 else
  mu=[-0.24*z-6.83,-0.23*z+1.53,-0.08*z+1.81,-0.04*z+0.92];
  sd=[0.16*z+0.28,0.11*z+0.33,0.05*z+0.3,-0.07*z+0.41];
  p.clusters=19;p.rTau=2.1;p.xprDb=[8 3];p.clusterSpread=[11 10 22 7];
  p.lspDistanceMeters=[13 NaN 10 10 9 10 10];
 end
 dz=d2dMeters/1000;
 if isLos,p.zsdLogMean=max(-0.21,-14.8*dz+0.01*8.5+0.83);
 else,p.zsdLogMean=max(-0.5,-3.1*dz+0.2);end
 p.zsdLogStd=0.35;
else
 if isLos
  mu=[-6.955-0.0963*q,1.06+0.1114*q,1.81,0.95];
  sd=[0.66,0.28,0.20,0.16];
  p.clusters=12;p.rTau=2.5;p.xprDb=[8 4];
  p.clusterSpread=[max(.25,6.5622-3.4084*q),5 11 7];
  p.lspDistanceMeters=[37 12 30 18 15 15 15];
 else
  mu=[-6.28-0.204*q,1.5-0.1144*q,2.08-0.27*q, ...
      -0.3236*q+1.512];
  sd=[0.39,0.28,0.11,0.16];
  p.clusters=20;p.rTau=2.3;p.xprDb=[7 3];
  p.clusterSpread=[max(.25,6.5622-3.4084*q),2 15 7];
  p.lspDistanceMeters=[50 NaN 40 50 50 50 50];
 end
 % Table 7.5-7 ZSD; h_UT=1.5 m for this scene.
 dz=d2dMeters/1000;
 if isLos,p.zsdLogMean=max(-.5,-2.1*dz+.75);p.zsdLogStd=.40;
 else,p.zsdLogMean=max(-.5,-2.1*dz+.9);p.zsdLogStd=.49;end
end
p.muLog=[mu(1:3),p.zsdLogMean,mu(4)];
p.sigmaLog=[sd(1:3),p.zsdLogStd,sd(4)];
p.raysPerCluster=20;p.clusterShadowDb=3;
p.scenario=scenario;p.isLos=isLos;p.carrierGhz=carrierGhz;
end
