function report=validate_srs_overlap_posterior(outPath)
%VALIDATE_SRS_OVERLAP_POSTERIOR Independent conditional-error Monte Carlo.
seed=1702; stream=RandStream('mt19937ar','Seed',seed);
K=3; L=1; Nr=4; Nt=2; T=16; trials=600;
R=complex(zeros(Nr,Nr,K,L));
for k=1:K, R(:,:,k,1)=eye(Nr)*(0.5+k/4); end
errorSum=complex(zeros(Nr,Nr,K)); predicted=zeros(Nr,Nr,K);
for t=1:trials
 H=complex(zeros(K,L,Nr,Nt)); innovation=H;
 for k=1:K
  h=sqrt(real(R(1,1,k,1))/2)*(randn(stream,Nr,Nt)+1i*randn(stream,Nr,Nt));
  e=sqrt(real(R(1,1,k,1))/2)*(randn(stream,Nr,Nt)+1i*randn(stream,Nr,Nt));
  H(k,1,:,:)=reshape(h,1,1,Nr,Nt);
  innovation(k,1,:,:)=reshape(e,1,1,Nr,Nt);
 end
 z=(randn(stream,Nt,L,T)+1i*randn(stream,Nt,L,T))/sqrt(2);
 s=srs_overlap_snapshot(H,R,[1 9;1 9;2 10],ones(K,2), ...
     0.2,z,20,120,innovation,[1 2],17,3,3.5);
 for k=1:K
  e=reshape(s.H_true(k,1,:,:)-s.H_hat(k,1,:,:),Nr,Nt);
  errorSum(:,:,k)=errorSum(:,:,k)+e*e'/Nt;
  predicted(:,:,k)=real(s.C(:,:,k,1));
 end
end
empirical=errorSum/trials;
relative=norm(empirical(:)-predicted(:))/norm(predicted(:));
report=struct('seed',seed,'trials',trials,'relativeFrobenius',relative, ...
 'empirical',empirical,'predicted',predicted, ...
 'scope','small Gaussian Jakes-lag posterior; not formal scene validation');
save(outPath,'report','-v7');
end
