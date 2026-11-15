function report=test_v6_joint_srs_posterior_monte_carlo(root)
% Independent Gaussian draws test the joint SRS posterior statistically.
% An algebraic five-source reconstruction alone cannot pass this test.
addpath(fullfile(root,'03_研究场景构造','matlab'));
K=2;L=1;Nr=2;Nt=2;N=1500;noise=.2;
Q=kron([1 .8;.8 1],eye(Nr));factor=chol(Q,'lower');
R=repmat(eye(Nr),1,1,K,L);pilot=[1 9;1 9];
energy=ones(K,2);priorMean=complex(zeros(K,L,Nr,Nt));
stream=RandStream('mt19937ar','Seed',20261115);
observed=complex(zeros(Nr,Nr,K));predicted=observed;
withoutCross=observed;crossSize=0;
for sample=1:N
 H=complex(zeros(K,L,Nr,Nt));
 for port=1:Nt
  draw=factor*complexNormal(stream,K*Nr,1);
  columns=reshape(draw,Nr,K);
  for k=1:K
   H(k,1,:,port)=reshape(columns(:,k),1,1,Nr,1);
  end
 end
 z=complexNormal(stream,Nt,L,16);
 snap=srs_overlap_snapshot(H,R,pilot,energy,noise,z, ...
     0,1,H,[1 2],16,3,3.5,priorMean,Q);
 r=struct('H',snap.H_hat,'C',permute(snap.C,[3 4 1 2]), ...
     'ueIds',(1:K)','srs_rows',[1 2], ...
     'srs_diagnostic',struct('pilotGram',snap.gram,'a',1), ...
     'age_ms',0);
 data=struct('modelId','srs_overlap_jakes_v3', ...
     'resources',r,'pilotId',pilot,'pilotEnergy',energy, ...
     'priorRx',R,'jointPrior',Q,'pilotNoiseVariance',noise);
 sources=srs_overlap_source_statistics(data,1);
 for k=1:K
  e=reshape(H(k,1,:,:)-snap.H_hat(k,1,:,:),Nr,Nt);
  observed(:,:,k)=observed(:,:,k)+e*e'/Nt/N;
  predicted(:,:,k)=predicted(:,:,k)+snap.C(:,:,k,1)/N;
  pieces=reshape(sum(sources.covariance(k,1,:,:,:),5),Nr,Nr);
  withoutCross(:,:,k)=withoutCross(:,:,k)+pieces/N;
 end
 crossSize=crossSize+norm(sources.crossSourceCovariance(:))/N;
end
relative=zeros(K,1);ablated=zeros(K,1);
for k=1:K
 target=predicted(:,:,k);actual=observed(:,:,k);
 relative(k)=norm(actual-target,'fro')/max(norm(target,'fro'),eps);
 ablated(k)=norm(actual-withoutCross(:,:,k),'fro')/ ...
     max(norm(target,'fro'),eps);
end
report=struct('independentDraws',N,'seed',20261115, ...
    'priorCrossFro',norm(Q(1:Nr,Nr+1:end),'fro'), ...
    'sourceCrossFro',crossSize, ...
    'empiricalErrorCovariance',observed, ...
    'predictedPosteriorCovariance',predicted, ...
    'empiricalRelativeMismatch',relative, ...
    'crossOmittedRelativeMismatch',ablated, ...
    'tolerance',.15,'passed',all(relative<.15) && ...
        crossSize>1e-5 && mean(ablated)>mean(relative));
assert(report.passed, ...
    'Joint SRS posterior did not statistically match independent draws');
end

function z=complexNormal(stream,varargin)
z=(randn(stream,varargin{:})+1i*randn(stream,varargin{:}))/sqrt(2);
end
