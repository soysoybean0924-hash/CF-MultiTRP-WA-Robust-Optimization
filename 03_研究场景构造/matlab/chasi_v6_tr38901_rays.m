function [ray,rx,tx]=chasi_v6_tr38901_rays(p,lsp,clusters,stream, ...
    trpXY,ueXY,velocityMps,bearingRad,carrierHz,spatialDraw)
% TR 38.901 v17.1.0 7.5 Steps 7-10, Tables 7.5-2/3/4.
% Outdoor 1.5 m UE and 10 m UMi / 25 m UMa BS, no O2I.
N=numel(clusters.power);M=20;lambda=299792458/carrierHz;
if strcmp(p.scenario,'UMi'),hbs=10;else,hbs=25;end
delta=[ueXY-trpXY,1.5-hbs];d3=norm(delta);d2=norm(delta(1:2));
azTx=atan2d(delta(2),delta(1));azRx=mod(azTx+180,360);
zeTx=acosd(delta(3)/d3);zeRx=180-zeTx;
keys=[4 5 8 10 11 12 14 15 16 19 20 25];
values=[.779 .860 1.018 1.090 1.123 1.146 1.190 1.211 ...
    1.226 1.273 1.289 1.358];
phi=interp1(keys,values,p.clusters,'linear');
theta=interp1([8 10 11 12 15 19 20 25], ...
    [.889 .957 1.031 1.104 1.1088 1.184 1.178 1.282], ...
    p.clusters,'linear');
if p.isLos
 k=lsp.kDb;
 phi=phi*(1.1035-.028*k-.002*k^2+.0001*k^3);
 theta=theta*(1.3086+.0339*k-.0077*k^2+.0002*k^3);
end
if ~isfield(clusters,'procedureB')
 assert(phi>0 && theta>0);
end
phi=max(phi,eps);theta=max(theta,eps);
ratio=max(realmin,clusters.anglePower/max(clusters.anglePower));
azScale=2*sqrt(-log(ratio))/phi/1.4;
zeScale=-log(ratio)/theta;
assert(nargin>=10,'Spatial ray fields are required for generated trace');
pick=clusters.keptOriginalIndices(:);
indexMatrix=(pick-1)*M+(1:M);
rayPick=reshape(indexMatrix.',[],1);
aoa=center(azScale*lsp.asaDeg,lsp.asaDeg,azRx, ...
    spatialDraw.centerSign(pick,1),spatialDraw.centerJitter(pick,1),p.isLos);
aod=center(azScale*lsp.asdDeg,lsp.asdDeg,azTx, ...
    spatialDraw.centerSign(pick,2),spatialDraw.centerJitter(pick,2),p.isLos);
zoa=center(zeScale*lsp.zsaDeg,lsp.zsaDeg,zeRx, ...
    spatialDraw.centerSign(pick,3),spatialDraw.centerJitter(pick,3),p.isLos);
zodOffset=0;
if ~p.isLos
 if strcmp(p.scenario,'UMi')
  zodOffset=-10^(-1.5*log10(max(10,d2))+3.3);
 else
  fc=max(6,carrierHz/1e9);a=.208*log10(fc)-.782;
  e=7.66*log10(fc)-5.96;
  zodOffset=e-10^(a*log10(max(25,d2)));
 end
end
zod=center(zeScale*lsp.zsdDeg,lsp.zsdDeg,zeTx+zodOffset, ...
    spatialDraw.centerSign(pick,4),spatialDraw.centerJitter(pick,4),p.isLos);
if isfield(clusters,'procedureB') && clusters.procedureB
 offset=clusters.centerOffsetDeg;
 assert(isequal(size(offset),[N,4]));
 aoa=azRx+offset(:,1);aod=azTx+offset(:,2);
 zoa=zeRx+offset(:,3);zod=zeTx+zodOffset+offset(:,4);
end
alpha=[.0447 .1413 .2492 .3715 .5129 .6797 .8844 ...
    1.1481 1.5195 2.1551];
alpha=reshape([alpha;-alpha],1,[]);
aoaRay=zeros(N,M);aodRay=aoaRay;zoaRay=aoaRay;zodRay=aoaRay;
for n=1:N
 aoaRay(n,:)=aoa(n)+p.clusterSpread(3)*alpha;
 strong=ismember(n,clusters.strongestClusterIndices);
 aodRay(n,:)=aod(n)+p.clusterSpread(2)* ...
     alpha(coupling(spatialDraw.couplingScore(:,1,pick(n)),strong));
 zoaRay(n,:)=zoa(n)+p.clusterSpread(4)* ...
     alpha(coupling(spatialDraw.couplingScore(:,2,pick(n)),strong));
 zodRay(n,:)=zod(n)+.375*10^p.zsdLogMean* ...
     alpha(coupling(spatialDraw.couplingScore(:,3,pick(n)),strong));
end
aoaRay=mod(aoaRay,360);aodRay=mod(aodRay,360);
zoaRay=abs(mod(zoaRay+180,360)-180);
zodRay=abs(mod(zodRay+180,360)-180);
ar=reshape(aoaRay.',[],1);at=reshape(aodRay.',[],1);
zr=reshape(zoaRay.',[],1);zt=reshape(zodRay.',[],1);
arrival=unitvec(ar,zr);departure=unitvec(at,zt);
power=clusters.rayPower;
if p.isLos,power=power/(1+10^(lsp.kDb/10));end
delay=clusters.rayDelaySeconds;
xpr=10.^((p.xprDb(1)+p.xprDb(2)* ...
    spatialDraw.xprNormal(rayPick))/10);
phase=-pi+2*pi*.5*(1+erf( ...
    spatialDraw.phaseNormal(rayPick,:)/sqrt(2)));
if p.isLos
 K=10^(lsp.kDb/10);
 power(end+1,1)=K/(K+1);
 delay(end+1,1)=0; % excess delay; geometric carrier phase below
 xpr(end+1,1)=Inf;
 carrierPhase=-2*pi*d3/lambda;
 phase(end+1,:)=[carrierPhase 0 0 carrierPhase+pi];
 arrival(end+1,:)=unitvec(azRx,zeRx);
 departure(end+1,:)=unitvec(azTx,zeTx);
 ar(end+1)=azRx;at(end+1)=azTx;zr(end+1)=zeRx;zt(end+1)=zeTx;
end
Q=numel(power);tx=struct('positionMeters',zeros(64,3), ...
    'fieldTheta',zeros(64,Q),'fieldPhi',zeros(64,Q));
rx=struct('positionMeters',zeros(4,3), ...
    'fieldTheta',zeros(4,Q),'fieldPhi',zeros(4,Q));
for j=1:Q
 panel=chasi_tr38901_array_panel(bearingRad,deg2rad(at(j)), ...
     deg2rad(zt(j)));
 tx.positionMeters=panel.positionWavelength.'*lambda;
 tx.fieldTheta(:,j)=panel.fieldTheta(:);
 tx.fieldPhi(:,j)=panel.fieldPhi(:);
 rx.positionMeters=([-.25 0 0;-.25 0 0;.25 0 0;.25 0 0])*lambda;
 rx.fieldTheta(:,j)=[1;1;1;1]/sqrt(2);
 rx.fieldPhi(:,j)=[-1;1;-1;1]/sqrt(2);
end
ray=struct('power',power,'delaySeconds',delay, ...
    'arrivalUnit',arrival,'departureUnit',departure, ...
    'xprLinear',xpr,'phasesRad',phase,'velocityMps',velocityMps, ...
    'spatialPhasesRad',phase, ...
    'clusterIndex',[repelem((1:N)',M,1);zeros(double(p.isLos),1)], ...
    'clusterSign',2*(spatialDraw.clusterSign(pick)>0)-1, ...
    'clusterArrivalUnit',unitvec(aoa,zoa), ...
    'clusterDepartureUnit',unitvec(aod,zod), ...
    'isLos',p.isLos, ...
    'clusterShadowDb',clusters.clusterShadowDb, ...
    'shadowFieldIndex',clusters.keptOriginalIndices, ...
    'delaySpreadSeconds',lsp.dsSeconds,'delayScaling',p.rTau, ...
    'delayPowerFactor',clusters.delayPowerFactor);
end
function a=center(offset,spread,geometric,signDraw,jitter,isLos)
signs=2*(signDraw(:)>0)-1;
a=signs.*offset+spread/7*jitter(:);
if isLos,a=a-a(1);end
a=a+geometric;
end
function v=unitvec(azDeg,zeDeg)
v=[sind(zeDeg(:)).*cosd(azDeg(:)), ...
   sind(zeDeg(:)).*sind(azDeg(:)),cosd(zeDeg(:))];
end
function index=coupling(score,strong)
if ~strong,[~,index]=sort(score);index=index.';return;end
index=1:20;
groups={[1:8 19:20],[9:12 17:18],13:16};
for j=1:3
 g=groups{j};[~,order]=sort(score(g));index(g)=g(order);
end
end
