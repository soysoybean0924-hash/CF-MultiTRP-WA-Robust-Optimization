function result=calibrate_arvidson_srs(cfg,outputPath)
%CALIBRATE_ARVIDSON_SRS Fit Eq.23 against an independent pilot estimator.
% This is a project-specific reference estimator, not Arvidson's LTE simulator.
% Each row uses a fresh channel and receiver-noise realization. No Eq.23
% perturbation is used to generate the training or validation target.
arguments
 cfg struct
 outputPath (1,:) char
end
required={'seed','users','trps','rxPorts','txPorts','resourceBlocks', ...
 'pilotLabels','noiseVariance','pilotEnergy','gain','trainSamples','testSamples'};
for i=1:numel(required), assert(isfield(cfg,required{i}),required{i}); end
assert(ismember(cfg.resourceBlocks,[16 17]) && cfg.users>cfg.pilotLabels);
assert(isequal(size(cfg.gain),[cfg.users cfg.trps]));
assert(cfg.noiseVariance>0 && cfg.pilotEnergy>0);
stream=RandStream('mt19937ar','Seed',cfg.seed);
% N frequency observations are projected onto a finite delay window.
N=12*cfg.resourceBlocks; D=min(12,N);
n=(0:N-1)'; labels=mod((0:cfg.users-1)',cfg.pilotLabels);
seq=exp(2i*pi*n*labels'/cfg.pilotLabels);
% Relative cyclic shifts are deliberately nonorthogonal after truncation.
delayBasis=exp(-2i*pi*n*(0:D-1)/N)/sqrt(N);
features=zeros(cfg.trainSamples+cfg.testSamples,4);
errors=zeros(size(features,1),1);
for t=1:size(features,1)
 k=randi(stream,cfg.users); l=randi(stream,cfg.trps);
 h=(randn(stream,cfg.users,1)+1i*randn(stream,cfg.users,1))/sqrt(2);
 h=h.*sqrt(cfg.gain(:,l));
 phase=exp(-2i*pi*n*randi(stream,D,1,cfg.users)/N);
 pilots=seq.*phase;
 y=sqrt(cfg.pilotEnergy)*(pilots*h) + ...
   sqrt(cfg.noiseVariance/2)*(randn(stream,N,1)+1i*randn(stream,N,1));
 matched=conj(pilots(:,k)).*y/sqrt(cfg.pilotEnergy);
 % Finite-resolution channel estimate: retained delay taps then reconstruct.
 impulse=delayBasis'*matched;
 reconstructed=delayBasis*impulse;
 estimate=reconstructed(1);
 truth=h(k);
 same=(labels==labels(k)); same(k)=false;
 other=~same; other(k)=false;
 features(t,:)=[abs(truth),sum(abs(h(same))),sum(abs(h(other))),1];
 errors(t)=abs(estimate-truth);
end
train=1:cfg.trainSamples; test=cfg.trainSamples+(1:cfg.testSamples);
% Nonnegative least squares: alpha, same-label beta, other-label delta, gamma.
coef=lsqnonneg(features(train,:),errors(train));
pred=features(test,:)*coef;
result=struct('model','Arvidson_Eq23_single_draw_probe_rejected', ...
 'validated',false, ...
 'fitTarget','single_draw_absolute_error_not_conditional_RMSE', ...
 'coefficientNames',{{'alpha_truncation','beta_same_label', ...
 'delta_other_label','gamma_background'}},'coefficients',coef, ...
 'trainMae',mean(abs(features(train,:)*coef-errors(train))), ...
 'testMae',mean(abs(pred-errors(test))), ...
 'testConstantMae',mean(abs(mean(errors(train))-errors(test))), ...
 'testNormalizedMae',mean(abs(pred-errors(test)))/mean(errors(test)), ...
 'testBias',mean(pred-errors(test)), ...
 'testCorrelation',corr(pred,errors(test)), ...
 'designRank',rank(features(train,:)), ...
 'config',cfg,'testPredicted',pred,'testObserved',errors(test), ...
 'testFeatures',features(test,:));
assert(result.designRank==4,'Calibration features are not identifiable');
save(outputPath,'result','-v7.3');
end
