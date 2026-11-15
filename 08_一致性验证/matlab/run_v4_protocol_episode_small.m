function out=run_v4_protocol_episode_small(root,inputPath,waCheckpoint,ttiCount)
% Short post-decision online channel/queue smoke before next CSI/SRS event.
% Explicitly not a formal 1-second protocol experiment.
startup;assert(ttiCount>=2 && ttiCount<=19);
x=load(inputPath,'data');d=x.data;assert(strcmp(d.modelId,'srs_csi_joint_v4'));
y=load(waCheckpoint,'state');s=y.state;
assert(strcmp(s.identity.inputSha256,chasi_sha256_file(inputPath)));
resources=d.resources(1:numel(s.precoder));U=d.sceneConfig.users;
stream=RandStream('mt19937ar','Seed',d.measurementSeed+803);
arrival=zeros(U,ttiCount);arrival(:,1)=12000;capacity=zeros(U,ttiCount);
rate=zeros(U,ttiCount);a=besselj(0,2*pi*(d.config.speedKmh/3.6)* ...
    (d.config.frequencyGhz*1e9)/299792458*(d.config.ttiMs/1000));
for t=1:ttiCount
 for g=1:numel(resources)
  r=resources(g);[K,L,Nr,Nt]=size(r.H_true);
  for k=1:K,for l=1:L
   prior=d.priorRx(:,:,double(r.ueIds(k)),l);
   z=(randn(stream,Nr,Nt)+1i*randn(stream,Nr,Nt))/sqrt(2);
   old=reshape(double(r.H_true(k,l,:,:)),Nr,Nt);
   r.H_true(k,l,:,:)=single(reshape(a*old+sqrt(1-a^2)* ...
       chol(prior,'lower')*z,1,1,Nr,Nt));
  end,end
  resources(g)=r;
 end
 rate(:,t)=chasi_v3_true_rates(resources,s.precoder,U);
 capacity(:,t)=floor(max(0,rate(:,t))*d.sceneConfig.bandwidthHz* ...
     d.config.ttiMs/1000);
end
bler=.05*ones(U,ttiCount);events=chasi_v4_rlc_event_sim( ...
    arrival,capacity,bler,d.trafficSeed);
counter=chasi_acked_rlc_seconds(events,d.config.ttiMs);
out=struct('inputSha256',chasi_sha256_file(inputPath), ...
    'waCheckpoint',waCheckpoint,'ttiCount',ttiCount, ...
    'trueRateByTti',rate,'capacityBits',capacity, ...
    'events',events,'counter',counter, ...
    'formalProtocolValidated',false, ...
    'scope','short online truth evolution, fixed beam until next measurement; modeled BLER');
fprintf('V4_EPISODE tti=%d attempts=%d acks=%d valid=%d\n', ...
    ttiCount,nnz(events.attemptedBits),nnz(events.ackedBits),nnz(counter.valid));
end
