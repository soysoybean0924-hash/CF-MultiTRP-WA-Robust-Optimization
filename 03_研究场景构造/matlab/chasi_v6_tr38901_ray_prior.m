function [Q,meanH,components]=chasi_v6_tr38901_ray_prior( ...
    ray,rx,tx,carrierHz,frequencyHz)
% LOS conditional mean and covariance of independently phased diffuse
% polarization components. Each of the four phases is one generating term.
if nargin<5,frequencyHz=0;end
N=numel(ray.power);Nr=size(rx.positionMeters,1);
Nt=size(tx.positionMeters,1);
assert(Nr==4 && Nt==64 && size(rx.fieldTheta,2)==N);
lambda=299792458/carrierHz;
Q=complex(zeros(Nr,Nr));meanH=complex(zeros(Nr,Nt));
components=complex(zeros(Nr,Nt,N,4));
polarization=[1 1;1 2;2 1;2 2];
for n=1:N
 receiveArray=exp(1i*2*pi/lambda* ...
     (double(rx.positionMeters)*double(ray.arrivalUnit(n,:)).'));
 transmitArray=exp(1i*2*pi/lambda* ...
     (double(tx.positionMeters)*double(ray.departureUnit(n,:)).'));
 receiveField=[rx.fieldTheta(:,n),rx.fieldPhi(:,n)];
 transmitField=[tx.fieldTheta(:,n),tx.fieldPhi(:,n)];
 shared=sqrt(double(ray.power(n)))* ...
     exp(-1i*2*pi*frequencyHz*double(ray.delaySeconds(n)))* ...
     (receiveArray*transmitArray.');
 for q=1:4
  a=polarization(q,1);b=polarization(q,2);
  coupling=1;
  if a~=b,coupling=1/sqrt(double(ray.xprLinear(n)));end
  h=shared.*(receiveField(:,a)*transmitField(:,b).')* ...
      coupling*exp(1i*double(ray.phasesRad(n,q)));
  components(:,:,n,q)=h;
  if ray.clusterIndex(n)==0
   meanH=meanH+h;
  else
   Q=Q+(h*h')/Nt;
  end
 end
end
Q=(Q+Q')/2;
assert(all(isfinite(Q(:))) && min(real(eig(Q)))>0, ...
    'Generated diffuse receive covariance must be positive definite');
end
