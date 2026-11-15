function report=validate_v5_csi_sources_small(root)
startup;p=fullfile(root,'03_研究场景构造','inputs', ...
    'srs_csi_joint_v5_prb50_attribution_smoke','input.mat');
assert(strcmp(chasi_sha256_file(p),strtrim(fileread([p '.sha256']))));
x=load(p,'data');d=x.data;assert(strcmp(d.modelId,'srs_csi_joint_v5'));
assert(d.prbEffective==.5 && all(d.prbSchedule.usedRbByTti(1:2:end)==136) ...
    && all(d.prbSchedule.usedRbByTti(2:2:end)==137));
r=d.resources(2);cdg=r.csi_diagnostic;
assert(numel(cdg.eventTti)>=1 && all(cdg.eventSourceMismatch<1e-6));
assert(all(cdg.eventMinPhysicalEigen>0));
[K,L,Nr,~,~]=size(r.sourceCov5);mismatch=0;white=zeros(K*L,1);n=0;
for k=1:K,for l=1:L
 C=reshape(double(r.C(k,l,:,:)),Nr,Nr);
 total=reshape(sum(double(r.sourceCov5(k,l,:,:,:)),5),Nr,Nr);
 mismatch=max(mismatch,norm(C-total,'fro')/max(norm(C,'fro'),eps));
 [F,flag]=chol((C+C')/2,'lower');assert(flag==0);
 e=reshape(double(r.H_true(k,l,:,:))-double(r.H(k,l,:,:)),Nr,[]);
 n=n+1;z=F\e;white(n)=mean(abs(z(:)).^2);
end,end
assert(mismatch<1e-6);
% Independent innovation/observation Monte Carlo checks the sequential
% Gaussian update and the finite-sample cross term, not external accuracy.
stream=RandStream('mt19937ar','Seed',721);draws=12000;
P=diag([.3 .5 .7 1]);N=.4*eye(4);J=P/(P+N);F=eye(4)-J;
e=chol(P,'lower')*cn(stream,4,draws);
noise=chol(N,'lower')*cn(stream,4,draws);
own=F*e;obs=-J*noise;total=own+obs;
pred=F*P*F'+J*N*J';emp=total*total'/draws;
cross=own*obs'/draws+obs*own'/draws;
empRel=norm(emp-pred,'fro')/norm(pred,'fro');
crossRel=norm(cross,'fro')/norm(pred,'fro');
assert(empRel<.05 && crossRel<.05);
report=struct('inputSha256',chasi_sha256_file(p),'prbEffective', ...
    d.prbEffective,'eventMismatchMax',max(cdg.eventSourceMismatch), ...
    'finalMismatch',mismatch,'whiteMean',mean(white(1:n)), ...
    'monteCarloRelative',empRel,'crossRelative',crossRel, ...
    'csiNoiseSeparated',true);
fprintf('V5_CSI prb=%.6f event=%.6g final=%.6g white=%.6g mc=%.6g cross=%.6g\n', ...
    report.prbEffective,report.eventMismatchMax,mismatch, ...
    report.whiteMean,empRel,crossRel);
end

function z=cn(stream,varargin)
z=(randn(stream,varargin{:})+1i*randn(stream,varargin{:}))/sqrt(2);
end
