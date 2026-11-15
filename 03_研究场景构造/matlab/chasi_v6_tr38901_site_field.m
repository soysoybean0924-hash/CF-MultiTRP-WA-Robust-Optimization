function out=chasi_v6_tr38901_site_field(mode,arg1,arg2,arg3,arg4)
% TR 38.901 v17.1.0 7.6.3.1 and 7.6.3.4, outdoor UMi/UMa.
% One frozen field per site and propagation condition. Co-sited sectors
% evaluate the same modes; different sites use independent seed streams.
if strcmp(mode,'create')
 seed=arg1;scenario=arg2;site=arg3;fc=arg4;
 assert(seed>=0 && seed==fix(seed) && site>=1 && site==fix(site));
 out=struct();out.scenario=scenario;out.site=site;
 stream=RandStream('mt19937ar','Seed',seed+30000+site);
 out.los=chasi_v6_tr38901_spectral_field('create',stream,50,1,64);
 out.condition=cell(1,2);
 for condition=0:1
  p=chasi_v6_tr38901_lsp_table(scenario,logical(condition),fc,10);
  N=p.clusters;
  assert(N>=1 && N==fix(N) && numel(p.lspDistanceMeters)==7);
  distance=p.lspDistanceMeters;distance(~isfinite(distance))=1;
  stream=RandStream('mt19937ar','Seed',seed+50000+10*site+condition);
  template=chasi_v6_tr38901_lsp(p,stream,zeros(7,1));
  cross=template.correlation;
  lspJoint=chasi_v6_tr38901_lsp_field( ...
      'create',stream,distance,cross,256);
  clusterDistance=ternary(strcmp(scenario,'UMi'), ...
      ternary(condition==1,12,15),ternary(condition==1,40,50));
  delayDistance=299792458*10^(p.muLog(1)+2*p.sigmaLog(1));
  % 7.6.3.2 B: delay/AOD/ZOD share the specified decorrelation length;
  % AOA/ZOA use 50 m, and cluster shadowing uses the SF length.
  % Fixed ray coupling/signs are retained from the first drop at reads.
  out.condition{condition+1}=struct('lspJoint',lspJoint, ...
      'lspFactor',chol(cross,'lower'), ...
      'delay',chasi_v6_tr38901_spectral_field( ...
          'create',stream,delayDistance,N,32), ...
      'angleArrival',chasi_v6_tr38901_spectral_field( ...
          'create',stream,50,2*N,32), ...
      'angleDeparture',chasi_v6_tr38901_spectral_field( ...
          'create',stream,delayDistance,2*N,32), ...
      'shadow',chasi_v6_tr38901_spectral_field( ...
          'create',stream,p.lspDistanceMeters(1),N,32), ...
      'small',chasi_v6_tr38901_spectral_field( ...
          'create',stream,clusterDistance,169*p.clusters,32), ...
      'clusters',p.clusters,'distanceMeters',clusterDistance);
 end
elseif strcmp(mode,'evaluate')
 field=arg1;xy=double(arg2(:).');condition=logical(arg3);
 assert(numel(xy)==2 && all(isfinite(xy)) && isscalar(condition));
 c=field.condition{double(condition)+1};N=c.clusters;
 correlated=chasi_v6_tr38901_lsp_field( ...
     'evaluate',c.lspJoint,xy);
 z=c.lspFactor\correlated;
 delayNormal=chasi_v6_tr38901_spectral_field( ...
     'evaluate',c.delay,xy);
 shadowNormal=chasi_v6_tr38901_spectral_field( ...
     'evaluate',c.shadow,xy);
 arrivalNormal=chasi_v6_tr38901_spectral_field( ...
     'evaluate',c.angleArrival,xy);
 departureNormal=chasi_v6_tr38901_spectral_field( ...
     'evaluate',c.angleDeparture,xy);
 angleNormal=[arrivalNormal(1:N),departureNormal(1:N), ...
     arrivalNormal(N+1:2*N),departureNormal(N+1:2*N)];
 row=chasi_v6_tr38901_spectral_field('evaluate',c.small,xy);
 assert(numel(row)==169*N);
 cursor=0;
 centerSign=reshape(row(cursor+1:cursor+4*N),N,4);cursor=cursor+4*N;
 centerJitter=reshape(row(cursor+1:cursor+4*N),N,4);cursor=cursor+4*N;
 couplingScore=reshape(row(cursor+1:cursor+60*N),20,3,N);
 cursor=cursor+60*N;
 xprNormal=row(cursor+1:cursor+20*N);cursor=cursor+20*N;
 phaseNormal=reshape(row(cursor+1:cursor+80*N),20*N,4);
 cursor=cursor+80*N;
 clusterSign=row(cursor+1:cursor+N);
 uniform=.5*(1+erf(delayNormal/sqrt(2)));
 out=struct('lspNormal',z, ...
     'ssp',struct('delayUniform',min(1-eps,max(eps,uniform)), ...
         'shadowNormal',shadowNormal, ...
         'angleUniform',min(1-eps,max(eps, ...
             .5*(1+erf(angleNormal/sqrt(2)))))), ...
     'ray',struct('centerSign',centerSign, ...
         'centerJitter',centerJitter,'couplingScore',couplingScore, ...
         'xprNormal',xprNormal,'phaseNormal',phaseNormal, ...
         'clusterSign',clusterSign));
elseif strcmp(mode,'los')
 field=arg1;xy=double(arg2(:).');trp=double(arg3(:).');
 assert(numel(xy)==2 && numel(trp)==2);
 d=max(10,norm(xy-trp));
 if strcmp(field.scenario,'UMi'),scale=36;else,scale=63;end
 probability=min(1,18/d+(1-18/d)*exp(-d/scale));
 normal=chasi_v6_tr38901_spectral_field('evaluate',field.los,xy);
 out=.5*(1+erf(normal/sqrt(2)))<probability;
else
 error('Unknown site-field mode');
end
end
function v=ternary(flag,a,b)
if flag,v=a;else,v=b;end
end
