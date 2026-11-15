function report=validate_srs_jakes_age(outPath)
%VALIDATE_SRS_JAKES_AGE Independent Gaussian single-lag correlation check.
seed=1722; stream=RandStream('mt19937ar','Seed',seed);
speedKmh=3; frequencyGhz=3.5; c=299792458;
fd=(speedKmh/3.6)*(frequencyGhz*1e9)/c;
ageMs=[0 20 100 320]; trials=10000;
h0=(randn(stream,trials,1)+1i*randn(stream,trials,1))/sqrt(2);
empirical=zeros(size(ageMs)); predicted=empirical;
for i=1:numel(ageMs)
 a=besselj(0,2*pi*fd*ageMs(i)/1000);
 innovation=(randn(stream,trials,1)+1i*randn(stream,trials,1))/sqrt(2);
 h1=a*h0+sqrt(1-a^2)*innovation;
 empirical(i)=real(mean(h1.*conj(h0)))/mean(abs(h0).^2);
 predicted(i)=a;
end
report=struct('seed',seed,'trials',trials,'speedKmh',speedKmh, ...
 'frequencyGhz',frequencyGhz,'dopplerHz',fd,'ageMs',ageMs, ...
 'predictedCorrelation',predicted,'empiricalCorrelation',empirical, ...
 'maxAbsoluteError',max(abs(predicted-empirical)), ...
 'scope','single-lag Jakes extension, not a full time-series generator');
save(outPath,'report','-v7');
end
