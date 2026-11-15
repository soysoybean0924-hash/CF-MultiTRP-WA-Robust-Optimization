function test_v6_joint_srs_posterior(root)
% Source-level regression to run when MATLAB is licensed.
addpath(fullfile(root,'03_研究场景构造','matlab'));
Nr=4;Nt=64;K=2;L=1;
H=complex(zeros(K,L,Nr,Nt));H(1,1,1,1)=1;
H(2,1,1,1)=2;
R=repmat(eye(Nr),1,1,K,L);
pilot=[1 9;1 9];energy=ones(K,2);z=complex(zeros(Nt,L,16));
joint=kron([1 .35;.35 1],eye(Nr));
args={H,R,pilot,energy,.1,z,0,1,H,[1 2],16,3,3.5, ...
    complex(zeros(size(H)))};
independent=srs_overlap_snapshot(args{:});
correlated=srs_overlap_snapshot(args{:},joint);
assert(norm(correlated.H_hat(:)-independent.H_hat(:))>1e-6);
assert(norm(correlated.C(:)-independent.C(:))>1e-6);
assert(correlated.priorCrossFro(1,2,1)>0);
resources=struct('H',correlated.H_hat,'C', ...
    permute(correlated.C,[3 4 1 2]),'ueIds',(1:K)', ...
    'srs_rows',[1 2], ...
    'srs_diagnostic',struct('pilotGram',correlated.gram,'a',1), ...
    'age_ms',0);
data=struct('modelId','srs_overlap_jakes_v3', ...
    'resources',resources,'pilotId',pilot,'pilotEnergy',energy, ...
    'priorRx',R,'jointPrior',joint,'pilotNoiseVariance',.1);
sources=srs_overlap_source_statistics(data,1);
assert(sources.maxRelativeTotalMismatch<1e-8);
assert(norm(sources.crossSourceCovariance(:))>0);
data.jointPrior=[];
data.resources.C=permute(independent.C,[3 4 1 2]);
independentSources=srs_overlap_source_statistics(data,1);
assert(norm(independentSources.crossSourceCovariance(:))==0 && ...
    independentSources.maxRelativeTotalMismatch<1e-8, ...
    'Independent UE phases must not create an artificial cross term');
for k=1:K
 physical=reshape(correlated.C(:,:,k,1),Nr,Nr);
 pieces=reshape(sum(sources.covariance(k,1,:,:,:),5),Nr,Nr)+ ...
     reshape(sources.crossSourceCovariance(k,1,:,:),Nr,Nr);
 assert(norm(physical-pieces,'fro')<1e-8);
end
end
