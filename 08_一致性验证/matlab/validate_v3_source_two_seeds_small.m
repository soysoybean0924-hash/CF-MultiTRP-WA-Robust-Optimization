function report=validate_v3_source_two_seeds_small(root)
% One group from each frozen input; internal consistency, not field accuracy.
startup;
names={'srs_overlap_jakes_v3_main','srs_overlap_jakes_v3_test_seed41'};
report=cell(1,2);
for i=1:2
    p=fullfile(root,'03_研究场景构造','inputs',names{i},'input.mat');
    assert(strcmp(chasi_sha256_file(p),strtrim(fileread([p '.sha256']))));
    x=load(p,'data');d=x.data;s=srs_overlap_source_statistics(d,1);
    r=d.resources(1);Nr=size(r.C,3);K=size(r.C,1);L=size(r.C,2);
    minEig=inf;maxHerm=0;white=zeros(K*L,1);n=0;
    for k=1:K,for l=1:L
        C=reshape(double(r.C(k,l,:,:)),Nr,Nr);
        maxHerm=max(maxHerm,norm(C-C','fro'));
        C=(C+C')/2;minEig=min(minEig,min(real(eig(C))));
        e=reshape(double(r.H_true(k,l,:,:))-double(r.H(k,l,:,:)),Nr,[]);
        [F,flag]=chol(C,'lower');assert(flag==0);
        z=F\e;n=n+1;white(n)=mean(abs(z(:)).^2);
        for q=1:4
            Q=reshape(s.covariance(k,l,:,:,q),Nr,Nr);
            assert(norm(Q-Q','fro')<1e-8*max(1,norm(Q,'fro')));
            assert(min(real(eig((Q+Q')/2)))>-1e-8);
        end
    end,end
    assert(s.maxRelativeTotalMismatch<1e-5 && minEig>0);
    report{i}=struct('input',names{i},'inputSha256',chasi_sha256_file(p), ...
        'links',n,'maxRelativeReconstruction',s.maxRelativeTotalMismatch, ...
        'minPhysicalEigenvalue',minEig,'maxHermitianResidual',maxHerm, ...
        'whitenedMean',mean(white(1:n)));
end
report=[report{:}];
disp(report);
end
