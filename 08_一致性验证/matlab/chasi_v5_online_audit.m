function audit=chasi_v5_online_audit(inputPath,checkpoint)
% Audit executed online allocation/events and final five-to-four C view.
x=load(inputPath,'data');d=x.data;y=load(checkpoint,'state');s=y.state;
assert(strcmp(s.identity.inputSha256,chasi_sha256_file(inputPath)));
T=s.iter;G=numel(d.resources);c=d.config;
actual=s.actualRbByGroupTti(:,1:T);
planned=d.prbSchedule.allocatedRbByGroupTti(1:G,1:T);
allocationError=max(abs(actual-planned),[],'all');
assert(allocationError<1e-8);
e=s.eventLog;assert(numel(e)==T);
for t=1:T
    absolute=c.decisionTti+t;
    assert(e(t).tti==absolute);
    expectedSrs=0;
    if mod(absolute,c.hopIntervalTti)==0
        expectedSrs=mod(absolute/c.hopIntervalTti,c.resourceGroups)+1;
    end
    assert(e(t).srsGroup==expectedSrs);
    assert(e(t).csiRs==(mod(absolute,c.csiRsPeriodTti)==0));
    assert(e(t).waTrigger==(e(t).csiRs || ...
        (e(t).srsGroup>0 && e(t).srsGroup<=G)));
    if e(t).waTrigger
        assert(e(t).waEffectiveTti==absolute);
    else
        assert(e(t).waEffectiveTti==0);
    end
end
sourceError=0;policyError=0;
for g=1:G
 r=s.resources(g);[K,L,Nr,~,~]=size(r.sourceCov5);
 for k=1:K,for l=1:L
  physical=reshape(double(r.C(k,l,:,:)),Nr,Nr);
  five=reshape(sum(double(r.sourceCov5(k,l,:,:,:)),5),Nr,Nr);
  four=reshape(sum(double(r.sourceCov(k,l,:,:,:)),5),Nr,Nr);
  sourceError=max(sourceError,norm(physical-five,'fro')/max(norm(physical,'fro'),eps));
  policyError=max(policyError,norm(physical-four,'fro')/max(norm(physical,'fro'),eps));
 end,end
 assert(r.srs_diagnostic.measurementTti<=c.decisionTti+T);
 assert(r.csi_diagnostic.lastCsiTti<=c.decisionTti+T);
end
assert(sourceError<1e-5 && policyError<1e-5);
logs=s.waLog;assert(all([logs.maxPowerViolationW]<1e-6));
waMask=~strcmp({logs.branch},'independent');
audit=struct('inputSha256',s.identity.inputSha256,'ttiCount',T, ...
    'allocationErrorRb',allocationError, ...
    'actualRbMean',mean(sum(actual,1)), ...
    'actualPrbFraction',mean(sum(actual,1))/273, ...
    'srsEvents',nnz([e.srsGroup]>0), ...
    'csiEvents',nnz([e.csiRs]), ...
    'waUpdates',numel(logs)/6,'fiveSourceError',sourceError, ...
    'fourSourceError',policyError, ...
    'nonconvergedUpdates',nnz(~[logs(waMask).converged]), ...
    'maxPowerViolationW',max([logs.maxPowerViolationW]));
fprintf('V5_AUDIT tti=%d rb=%g prb=%g srs=%d csi=%d wa=%d source=%g policy=%g nonconv=%d power=%g\n', ...
 T,audit.actualRbMean,audit.actualPrbFraction,audit.srsEvents, ...
 audit.csiEvents,audit.waUpdates,sourceError,policyError, ...
 audit.nonconvergedUpdates,audit.maxPowerViolationW);
end
