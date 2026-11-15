function report=validate_srs_overlap_input(root,outPath)
%VALIDATE_SRS_OVERLAP_INPUT Verify formal input identity and available stats.
path=fullfile(root,'03_研究场景构造','inputs', ...
 'srs_overlap_jakes_v3_main','input.mat');
assert(strcmp(strtrim(fileread([path '.sha256'])),chasi_sha256_file(path)));
s=load(path,'data'); d=s.data; r=d.resources;
assert(strcmp(d.modelId,'srs_overlap_jakes_v3'));
assert(numel(r)==17 && d.sceneConfig.users==210 && d.sceneConfig.trps==21);
assert(isequal(d.srsSenderIds(:),(1:210)'));
assert(d.sceneSeed==19 && d.measurementSeed==1701 && d.initialSeed==901);
age=zeros(1,17); coefficient=zeros(1,17); minEig=inf;
width=zeros(1,17);
for g=1:17
 q=r(g); assert(strcmp(q.srs_diagnostic.ageModel, ...
     'jakes_isotropic_single_lag_extension'));
 assert(isequal(q.srs_diagnostic.srsSenderIds(:),(1:210)'));
 assert(size(q.H,2)==21 && isequal(size(q.H),size(q.H_true)));
 age(g)=q.srs_diagnostic.ageMs;
 width(g)=q.srs_rb_count;
 coefficient(g)=q.srs_diagnostic.a;
 fd=(d.sceneConfig.speedKmh/3.6)*(d.config.frequencyGhz*1e9)/299792458;
 assert(abs(coefficient(g)-besselj(0,2*pi*fd*age(g)/1000))<1e-12);
 c=squeeze(double(q.C(1,1,:,:)));
 minEig=min(minEig,min(real(eig((c+c')/2))));
end
assert(minEig>=-1e-7);
assert(isequal(width,[repmat(16,1,16),17]) && sum(width)==273);
report=struct('modelId',d.modelId,'sha256',chasi_sha256_file(path), ...
 'groups',numel(r),'users',d.sceneConfig.users,'trps',d.sceneConfig.trps, ...
 'ageMs',age,'srsRbCount',width,'jakesCoefficient',coefficient, ...
 'minimumSampleCovarianceEigenvalue',minEig, ...
 'sceneSeed',d.sceneSeed,'measurementSeed',d.measurementSeed, ...
 'initialSeed',d.initialSeed,'srsSendersPerGroup',210);
save(outPath,'report','-v7');
end
