function out=chasi_v6_tr38901_clusters(p,lsp,stream,spatialDraw)
% TR 38.901 7.5 Steps 5-6 or 7.6.3.2 Procedure B Steps 5-7.
N=p.clusters;M=p.raysPerCluster;
assert(N>=2 && M==20 && lsp.dsSeconds>0);
if nargin>=4 && isfield(spatialDraw,'angleUniform')
 out=procedureB(p,lsp,spatialDraw);return
end
if nargin<4
 uniform=rand(stream,N,1);shadowNormal=randn(stream,N,1);
else
 uniform=spatialDraw.delayUniform(:);
 shadowNormal=spatialDraw.shadowNormal(:);
 assert(isequal(size(uniform),[N,1]) && ...
     isequal(size(shadowNormal),[N,1]) && ...
     all(uniform>0 & uniform<1));
end
raw=-p.rTau*lsp.dsSeconds*log(max(realmin,uniform));
delay=sort(raw-min(raw));
powerDelay=delay;Ctau=1;
if p.isLos
 K=10^(lsp.kDb/10);
 Ctau=.7705-.0433*lsp.kDb+.0002*lsp.kDb^2+ ...
     .000017*lsp.kDb^3;
 assert(Ctau>0);delay=delay/Ctau;
else,K=0;end
shadowDb=p.clusterShadowDb*shadowNormal;
unnormalized=exp(-powerDelay*(p.rTau-1)/(p.rTau*lsp.dsSeconds)).* ...
    10.^(-shadowDb/10);
power=unnormalized/sum(unnormalized);
anglePower=power/(K+1);
if p.isLos,anglePower(1)=anglePower(1)+K/(K+1);end
keep=anglePower>=max(anglePower)*10^(-25/10);
if p.isLos,keep(1)=true;end
power=power(keep);anglePower=anglePower(keep);
delay=delay(keep);shadowDb=shadowDb(keep);
rayDelay=repelem(delay,M,1);
[~,order]=sort(power,'descend');
if nargin>=4 && isfield(spatialDraw,'fixedStrongestOriginalIndices')
 fixed=spatialDraw.fixedStrongestOriginalIndices(:);
 original=find(keep);
 [present,local]=ismember(fixed,original);
 order=local(present);
end
sub=[zeros(8,1);ones(4,1);2*ones(4,1);ones(2,1);zeros(2,1)];
for j=1:min(2,numel(order))
 idx=(order(j)-1)*M+(1:M);
 rayDelay(idx)=rayDelay(idx)+1e-9*p.clusterSpread(1)*1.28*sub;
end
out=struct('delaySeconds',delay,'power',power, ...
    'delayPowerFactor',Ctau,'keptOriginalIndices',find(keep), ...
    'anglePower',anglePower,'clusterShadowDb',shadowDb, ...
    'rayDelaySeconds',rayDelay,'rayPower',repelem(power/M,M,1), ...
    'strongestClusterIndices',order(1:min(2,numel(order))));
end

function out=procedureB(p,lsp,draw)
% 7.6-15..17f: position-dependent delays, angle offsets and joint powers.
N=p.clusters;M=p.raysPerCluster;
assert(isequal(size(draw.angleUniform),[N,4]) && ...
    isequal(size(draw.delayUniform),[N,1]));
delayMaximum=10^(p.muLog(1)+2*p.sigmaLog(1));
delay=delayMaximum*draw.delayUniform(:);
delay=delay-min(delay);
if p.isLos,delay(1)=0;end
spread=[lsp.asaDeg,lsp.asdDeg,lsp.zsaDeg,lsp.zsdDeg];
mu=[p.muLog(3),p.muLog(2),p.muLog(5),p.muLog(4)];
sigma=[p.sigmaLog(3),p.sigmaLog(2), ...
    p.sigmaLog(5),p.sigmaLog(4)];
angle=(2*draw.angleUniform-1).*repmat(10.^(mu+2*sigma),N,1);
if p.isLos,angle(1,:)=0;end
K=0;if p.isLos,K=10^(lsp.kDb/10);end
effectiveDS=lsp.dsSeconds*(1+2*K);
effectiveSpread=spread*(1+K);
shadowDb=p.clusterShadowDb*draw.shadowNormal(:);
logPower=-delay/effectiveDS- ...
    sqrt(2)*sum(abs(angle)./effectiveSpread,2)- ...
    log(10)*shadowDb/10;
power=exp(logPower-max(logPower));power=power/sum(power);
anglePower=power/(1+K);
if p.isLos,anglePower(1)=anglePower(1)+K/(1+K);end
rayDelay=repelem(delay,M,1);
[~,order]=sort(power,'descend');
if isfield(draw,'fixedStrongestOriginalIndices')
 order=draw.fixedStrongestOriginalIndices(:);
end
sub=[zeros(8,1);ones(4,1);2*ones(4,1);ones(2,1);zeros(2,1)];
for j=1:min(2,numel(order))
 idx=(order(j)-1)*M+(1:M);
 rayDelay(idx)=rayDelay(idx)+1e-9*p.clusterSpread(1)*1.28*sub;
end
out=struct('delaySeconds',delay,'power',power, ...
    'delayPowerFactor',1,'keptOriginalIndices',(1:N)', ...
    'anglePower',anglePower,'clusterShadowDb',shadowDb, ...
    'rayDelaySeconds',rayDelay,'rayPower',repelem(power/M,M,1), ...
    'strongestClusterIndices',order(1:min(2,numel(order))), ...
    'centerOffsetDeg',angle,'procedureB',true);
end
