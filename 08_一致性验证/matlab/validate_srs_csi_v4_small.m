function report=validate_srs_csi_v4_small(root)
startup;
p=fullfile(root,'03_研究场景构造','inputs', ...
    'srs_csi_joint_v4_prb50_smoke','input.mat');
assert(strcmp(chasi_sha256_file(p),strtrim(fileread([p '.sha256']))));
x=load(p,'data');d=x.data;
assert(strcmp(d.modelId,'srs_csi_joint_v4') && d.schemaVersion==4);
assert(d.prbRequested==.5 && abs(d.prbEffective-.5)<.04);
assert(numel(d.resources)==2 && ~d.fullScaleInput);
mis=zeros(1,2);events=zeros(1,2);
for g=1:2
 r=d.resources(g);s=srs_overlap_source_statistics(d,g);
 assert(size(r.C,3)==4 && size(r.H,4)==64);
 assert(s.maxRelativeTotalMismatch<1e-5);
 e=r.csi_diagnostic.eventTti;
 assert(all(e>r.srs_diagnostic.measurementTti) && ...
     all(e<=r.srs_diagnostic.decisionTti));
 assert(r.csi_diagnostic.lastCsiTti<=d.config.decisionTti);
 mis(g)=s.maxRelativeTotalMismatch;events(g)=numel(e);
end
assert(events(2)>0,'Second group must exercise CSI-RS update');
report=struct('inputSha256',chasi_sha256_file(p),'modelId',d.modelId, ...
    'prbRequested',d.prbRequested,'prbEffective',d.prbEffective, ...
    'csiEventCount',events,'sourceMismatch',mis);
fprintf('V4_SMALL prb=%.6g events=%s mismatch=%s\n', ...
    report.prbEffective,mat2str(events),mat2str(mis,6));
end
