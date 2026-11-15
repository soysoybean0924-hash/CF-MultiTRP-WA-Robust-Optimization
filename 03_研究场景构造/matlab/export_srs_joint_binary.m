function export_srs_joint_binary(inputPath,outDir)
% Portable column-major array export; exact same frozen MATLAB input.
assert(isfile(inputPath) && ~isfolder(outDir)); mkdir(outDir);
x=load(inputPath,'data'); d=x.data;
meta=struct('modelId',d.modelId,'inputSha256',chasi_sha256_file(inputPath), ...
    'configSha256',d.configSha256,'config',d.config,'groups',numel(d.resources), ...
    'pilotNoiseVariance',double(d.pilotNoiseVariance), ...
    'srsSenderCount',numel(d.srsSenderIds));
put(outDir,'pilotId',double(d.pilotId));
put(outDir,'pilotEnergy',double(d.pilotEnergy));
put(outDir,'priorRx',double(d.priorRx));
for g=1:numel(d.resources)
    r=d.resources(g); prefix=sprintf('g%02d_',g);
    put(outDir,[prefix 'H_true'],double(r.H_true));
    put(outDir,[prefix 'H_hat'],double(r.H));
    put(outDir,[prefix 'C'],double(r.C));
    if isfield(r,'sourceCov') && ~isempty(r.sourceCov)
        assert(max(abs(sum(double(r.sourceCov),5)-double(r.C)),[],'all')<1e-5, ...
            'Source covariance must reconstruct total posterior C');
        put(outDir,[prefix 'sourceCov'],double(r.sourceCov));
    end
    put(outDir,[prefix 'serve'],double(r.servingMask));
    put(outDir,[prefix 'ueIds'],double(r.ueIds(:)));
    put(outDir,[prefix 'initial'],double(d.initialPrecoder{g}));
    errorPower=sum(abs(double(r.H_true)-double(r.H)).^2,'all');
    truePower=sum(abs(double(r.H_true)).^2,'all');
    theoryErrorPower=size(r.H,4)*real(sum( ...
        double(r.C(:,:,1,1))+double(r.C(:,:,2,2))+ ...
        double(r.C(:,:,3,3))+double(r.C(:,:,4,4)),'all'));
    meta.resources(g)=struct('users',numel(r.ueIds), ...
        'trps',size(r.H,2),'rx',size(r.H,3),'tx',size(r.H,4), ...
        'rank',d.config.rank,'resourceWeightNu',double(r.resourceWeightNu), ...
        'rbCount',double(r.rb_count),'ageMs',double(r.age_ms), ...
        'measurementTti',double(r.srs_diagnostic.measurementTti), ...
        'decisionTti',double(r.srs_diagnostic.decisionTti), ...
        'measurementTimeMs',double(r.srs_diagnostic.measurementTimeMs), ...
        'decisionTimeMs',double(r.srs_diagnostic.decisionTimeMs), ...
        'correlationA',double(r.srs_diagnostic.a), ...
        'sampleNmse',errorPower/max(truePower,realmin), ...
        'theoryNmse',theoryErrorPower/max(truePower,realmin), ...
        'theoryErrorPower',theoryErrorPower); %#ok<AGROW>
end
files=dir(fullfile(outDir,'*.bin'));
for j=1:numel(files)
    meta.binaryFiles(j)=struct('name',files(j).name, ...
        'sha256',chasi_sha256_file(fullfile(outDir,files(j).name))); %#ok<AGROW>
end
fid=fopen(fullfile(outDir,'manifest.json'),'wt'); assert(fid>=0);
fprintf(fid,'%s',jsonencode(meta)); fclose(fid);
end

function put(dir,name,a)
fid=fopen(fullfile(dir,[name '.real.bin']),'wb'); assert(fid>=0);
fwrite(fid,real(a(:)),'double'); fclose(fid);
if ~isreal(a)
    fid=fopen(fullfile(dir,[name '.imag.bin']),'wb'); assert(fid>=0);
    fwrite(fid,imag(a(:)),'double'); fclose(fid);
end
end

