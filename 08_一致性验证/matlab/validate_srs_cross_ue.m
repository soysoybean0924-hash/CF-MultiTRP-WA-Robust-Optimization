function report=validate_srs_cross_ue(outPath)
%VALIDATE_SRS_CROSS_UE Empirical cross-user posterior block, small problem.
stream=RandStream('mt19937ar','Seed',1723);
K=3; L=1; Nr=4; Nt=2; T=16; trials=12000;
pilotId=[1 9;1 9;2 10]; energy=ones(K,2);
R=complex(zeros(Nr,Nr,K,L));
for k=1:K, R(:,:,k,1)=eye(Nr)*(0.5+k/4); end
empirical=complex(zeros(Nr,Nr)); predicted=[];
for t=1:trials
 H=complex(zeros(K,L,Nr,Nt)); innovation=H;
 for k=1:K
  scale=sqrt(real(R(1,1,k,1))/2);
  H(k,1,:,:)=reshape(scale*(randn(stream,Nr,Nt)+1i*randn(stream,Nr,Nt)),1,1,Nr,Nt);
  innovation(k,1,:,:)=reshape(scale*(randn(stream,Nr,Nt)+1i*randn(stream,Nr,Nt)),1,1,Nr,Nt);
 end
 z=(randn(stream,Nt,L,T)+1i*randn(stream,Nt,L,T))/sqrt(2);
 s=srs_overlap_snapshot(H,R,pilotId,energy,0.2,z,20,120, ...
     innovation,[1 2],17,3,3.5);
 e1=reshape(s.H_true(1,1,:,:)-s.H_hat(1,1,:,:),Nr,Nt);
 e2=reshape(s.H_true(2,1,:,:)-s.H_hat(2,1,:,:),Nr,Nt);
 empirical=empirical+e1*e2'/Nt;
 if t==1
  data=struct('modelId','srs_overlap_jakes_v3','pilotId',pilotId, ...
   'pilotEnergy',energy,'priorRx',R, ...
   'config',struct('ulPortToDlRowsOneBased',[1 2]));
  data.resources=struct('srs_diagnostic',struct('pilotGram',s.gram, ...
      'noiseCov',s.noiseCov,'a',s.a));
  predicted=srs_overlap_cross_ue(data,1,1,1,2);
 end
end
empirical=empirical/trials;
report=struct('seed',1723,'trials',trials, ...
 'predicted',predicted,'empirical',empirical, ...
 'relativeFrobenius',norm(predicted-empirical,'fro')/max(norm(predicted,'fro'),eps), ...
 'scope','small Gaussian posterior cross-UE block, not formal scene calibration');
save(outPath,'report','-v7');
end
