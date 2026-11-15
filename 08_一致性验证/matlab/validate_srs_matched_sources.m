function report=validate_srs_matched_sources(outPath)
%VALIDATE_SRS_MATCHED_SOURCES Controlled conditional RMS source check.
% This validates a flat-channel matched-filter reference, not the LTE
% frequency-selective estimator or the formal 2T4R SRS input.
seed=1701; stream=RandStream('mt19937ar','Seed',seed);
trials=2000;
noiseLevels=[0.001 0.01 0.1]; powers=[0 0.25 1 4];
rows=zeros(4*numel(noiseLevels)*numel(powers),9); ix=0;
for rb=[16 17]
N=6*rb; n=(0:N-1)'; target=ones(N,1);
for reuse=[false true]
 if reuse, interferer=target;
 else, interferer=exp(2i*pi*n/8); end
 rho=mean(conj(target).*interferer);
 for nv=noiseLevels
 for power=powers
  ix=ix+1; errors=complex(zeros(trials,1));
  for t=1:trials
   h=(randn(stream)+1i*randn(stream))/sqrt(2);
   w=sqrt(nv/2)*(randn(stream,N,1)+1i*randn(stream,N,1));
   y=target+sqrt(power)*interferer*h+w;
   errors(t)=mean(conj(target).*y)-1;
  end
  predicted=sqrt(power*abs(rho)^2+nv/N);
  empirical=sqrt(mean(abs(errors).^2));
  rows(ix,:)=[rb,double(reuse),nv,power,abs(rho),predicted,empirical, ...
    abs(predicted-empirical)/predicted,mean(abs(errors))];
 end
end
end
end
report=struct('model','matched_filter_flat_channel', ...
 'seed',seed,'observationsPerCondition',trials, ...
 'columns',{{'resourceBlocks','samePilot','noiseVariance','interfererPower','pilotOverlap', ...
 'predictedComplexRms','empiricalComplexRms','relativeRmsError', ...
 'meanSingleDrawMagnitude'}},'rows',rows, ...
 'maxRelativeRmsError',max(rows(:,8)));
save(outPath,'report','-v7');
end
