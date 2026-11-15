function report=validate_srs_source_components_small(outPath)
%VALIDATE_SRS_SOURCE_COMPONENTS_SMALL Independent generating-draw check.
% Controlled synthetic five-sender pilot, no fitted-model error generator.
stream=RandStream('mt19937ar','Seed',2261);
rb=17; n=(0:6*rb-1)'; shifts=(0:7);
seq=exp(2i*pi*n*shifts/8); gram=blkdiag(seq'*seq/(6*rb),seq'*seq/(6*rb));
gram=(gram+gram')/2; T=16; Nr=4; K=5; draws=12000;
pilot=[1 9;1 9;2 10;8 16;3 11];
energy=[1 1;1.4 .8;.5 1.2;1.1 .7;.7 1.5]; rows=[1 2];
gain=[.6 1.3 .9 .4 1.1]; base=.35.^abs((1:Nr)'-(1:Nr));
R=complex(zeros(Nr,Nr,K)); M=complex(zeros(T,Nr,K));
for j=1:K
    R(:,:,j)=gain(j)*base;
    for p=1:2
        M(:,rows(p),j)=M(:,rows(p),j)+ ...
            sqrt(energy(j,p))*gram(:,pilot(j,p));
    end
end
noiseCov=.02*gram; measure=complex(zeros(T));
for j=1:K, measure=measure+M(:,:,j)*R(:,:,j)*M(:,:,j)'; end
k=1; J=M(:,:,k); A=(R(:,:,k)*J')/(measure+noiseCov); a=.78;
own=eye(Nr)-A*J;
analytic=cat(3,a^2*own*R(:,:,k)*own', ...
    a^2*A*(measure-J*R(:,:,k)*J')*A', ...
    a^2*A*noiseCov*A', (1-a^2)*R(:,:,k));
X=complex(zeros(Nr,draws,K));
for j=1:K
    X(:,:,j)=chol(R(:,:,j),'lower')*cn(stream,Nr,draws);
end
other=complex(zeros(T,draws));
for j=1:K
    if j~=k, other=other+M(:,:,j)*X(:,:,j); end
end
noise=chol(noiseCov,'lower')*cn(stream,T,draws);
innovation=chol(R(:,:,k),'lower')*cn(stream,Nr,draws);
term=cat(3,a*own*X(:,:,k),-a*A*other, ...
    -a*A*noise,sqrt(1-a^2)*innovation);
relative=zeros(1,4);
for q=1:4
    sample=term(:,:,q)*term(:,:,q)'/draws;
    relative(q)=norm(sample-analytic(:,:,q),'fro')/ ...
        max(norm(analytic(:,:,q),'fro'),eps);
end
total=sum(term,3); totalEmp=total*total'/draws;
totalAnalytic=sum(analytic,3);
totalRelative=norm(totalEmp-totalAnalytic,'fro')/ ...
    norm(totalAnalytic,'fro');
cross=totalEmp;
for q=1:4, zq=term(:,:,q); cross=cross-zq*zq'/draws; end
crossRelative=norm(cross,'fro')/norm(totalAnalytic,'fro');
assert(max(relative)<.12 && totalRelative<.08, ...
    'Independent source draws disagree with analytic covariance');
assert(crossRelative<.08,'Independent generating terms have large cross covariance');
report=struct('seed',2261,'draws',draws,'rb',rb, ...
    'pilotId',pilot,'energy',energy,'relativeErrorBySource',relative, ...
    'totalRelativeError',totalRelative,'crossRelativeError',crossRelative, ...
    'scope','controlled Gaussian pilot model; not formal field calibration');
save(outPath,'report','-v7');
fprintf('SOURCE_MC sourceRel=%s totalRel=%.5g\n', ...
    mat2str(relative,4),totalRelative);
end

function x=cn(stream,varargin)
x=(randn(stream,varargin{:})+1i*randn(stream,varargin{:}))/sqrt(2);
end
