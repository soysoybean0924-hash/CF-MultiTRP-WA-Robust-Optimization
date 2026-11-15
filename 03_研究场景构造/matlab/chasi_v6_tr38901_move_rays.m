function moved=chasi_v6_tr38901_move_rays(ray,timeSeconds,referenceDistanceMeters,carrierHz)
% TR 38.901 v17.1.0 7.6.3.2 Procedure A, static BS, straight UE path.
% Cluster-center discretization of 7.6-9 and 7.6-11..14. Fixed ray
% offsets rotate with each cluster. Time-varying spatial fields remain a gate.
assert(timeSeconds>=0 && isfinite(timeSeconds) && ...
    referenceDistanceMeters>0 && carrierHz>0);
c0=299792458;lambda=c0/carrierHz;
v=double(ray.velocityMps(:));moved=ray;
if timeSeconds==0,return;end
steps=max(1,ceil(norm(v)*timeSeconds/.5));dt=timeSeconds/steps;
N=numel(ray.power);clusterCount=numel(ray.clusterShadowDb);
absoluteDelay=double(ray.delaySeconds(:))+referenceDistanceMeters/c0;
if ray.isLos
 losInitial=double(ray.arrivalUnit(N,:).');
 losPosition=referenceDistanceMeters*losInitial;
end
for step=1:steps
 for j=1:clusterCount
  members=find(ray.clusterIndex==j);
  a=double(moved.clusterArrivalUnit(j,:).');
  d=double(moved.clusterDepartureUnit(j,:).');
  [phiA,thetaA]=angles(a);[phiD,thetaD]=angles(d);
  R=rotZ(phiD+pi)*rotY(pi/2-thetaD)* ...
      diag([1 ray.clusterSign(j) 1])*rotY(pi/2-thetaA)*rotZ(-phiA);
  vr=R*v;vt=-v;
  radius=max(c0*mean(absoluteDelay(members)),referenceDistanceMeters/10);
  phiD=phiD+dot(vr,basisPhi(phiD))/ ...
      max(radius*sin(thetaD),1e-9)*dt;
  thetaD=thetaD+dot(vr,basisTheta(thetaD,phiD))/radius*dt;
  phiA=phiA+dot(vt,basisPhi(phiA))/ ...
      max(radius*sin(thetaA),1e-9)*dt;
  thetaA=thetaA+dot(vt,basisTheta(thetaA,phiA))/radius*dt;
  newA=unit(phiA,thetaA);newD=unit(phiD,thetaD);
  Ra=rotationBetween(a,newA);Rd=rotationBetween(d,newD);
  absoluteDelay(members)=absoluteDelay(members)-dot(a,v)*dt/c0;
  for n=members(:).'
   oldA=double(moved.arrivalUnit(n,:).');
   moved.phasesRad(n,:)=moved.phasesRad(n,:)+ ...
       2*pi*dot(oldA,v)*dt/lambda;
   moved.arrivalUnit(n,:)=(Ra*oldA).';
   moved.departureUnit(n,:)=(Rd* ...
       double(moved.departureUnit(n,:).')).';
  end
  moved.clusterArrivalUnit(j,:)=newA.';
  moved.clusterDepartureUnit(j,:)=newD.';
 end
 if ray.isLos
  oldDistance=norm(losPosition);
  losPosition=losPosition-v*dt;
  newDistance=norm(losPosition);
  assert(newDistance>0,'UE trajectory crosses the TRP location');
  newA=losPosition/newDistance;
  moved.arrivalUnit(N,:)=newA.';
  moved.departureUnit(N,:)=-newA.';
  absoluteDelay(N)=newDistance/c0;
  moved.phasesRad(N,:)=moved.phasesRad(N,:)+ ...
      2*pi*(oldDistance-newDistance)/lambda;
 end
end
moved.delaySeconds=max(0,absoluteDelay-min(absoluteDelay));
if isfield(ray,'clusterIndex') && isfield(ray,'clusterShadowDb')
 weights=zeros(clusterCount,1);
 shadowDb=ray.clusterShadowDb;
 if isfield(ray,'shadowField')
  initial=chasi_v6_tr38901_spectral_field( ...
      'evaluate',ray.shadowField,ray.initialXY);
  displaced=ray.initialXY+timeSeconds*v(1:2).';
  current=chasi_v6_tr38901_spectral_field( ...
      'evaluate',ray.shadowField,displaced);
  idx=ray.shadowFieldIndex;
  shadowDb=shadowDb+ray.clusterShadowSigmaDb* ...
      (current(idx)-initial(idx));
 end
 diffuse=ray.clusterIndex>0;
 for j=1:clusterCount
  members=ray.clusterIndex==j;
  delay=min(moved.delaySeconds(members))*ray.delayPowerFactor;
  weights(j)=exp(-delay*(ray.delayScaling-1)/ ...
      (ray.delayScaling*ray.delaySpreadSeconds))* ...
      10^(-shadowDb(j)/10);
 end
 weights=weights/sum(weights);
 diffuseTotal=sum(ray.power(diffuse));
 for j=1:clusterCount
  members=ray.clusterIndex==j;
  moved.power(members)=diffuseTotal*weights(j)/nnz(members);
 end
end
end
function [phi,theta]=angles(v)
phi=atan2(v(2),v(1));theta=acos(max(-1,min(1,v(3))));
end
function v=unit(phi,theta)
v=[sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];
end
function v=basisPhi(phi)
v=[-sin(phi);cos(phi);0];
end
function v=basisTheta(theta,phi)
v=[cos(theta)*cos(phi);cos(theta)*sin(phi);-sin(theta)];
end
function R=rotZ(a)
R=[cos(a) -sin(a) 0;sin(a) cos(a) 0;0 0 1];
end
function R=rotY(a)
R=[cos(a) 0 sin(a);0 1 0;-sin(a) 0 cos(a)];
end
function R=rotationBetween(a,b)
axis=cross(a,b);s=norm(axis);c=max(-1,min(1,dot(a,b)));
if s<1e-12
 if c>0,R=eye(3);
 else
  [~,i]=min(abs(a));basis=zeros(3,1);basis(i)=1;
  axis=cross(a,basis);axis=axis/norm(axis);
  R=2*(axis*axis.')-eye(3);
 end
 return
end
K=[0 -axis(3) axis(2);axis(3) 0 -axis(1); ...
   -axis(2) axis(1) 0];
R=eye(3)+K+K*K*((1-c)/s^2);
end
