function out=chasi_measure_target_group_static(resource,config)
%CHASI_MEASURE_TARGET_GROUP_STATIC Explicit project-hypothesis pilot chain.
% Uses only H_true, prior_gain and omitted_gain from a target scene group.
% The target's precomputed H and C are deliberately not reused. This is a
% static, calibrated-port, Gaussian model; it is not a physical SRS claim.
required={'seed','pilotSequences','pilotMatchedEnergy','pilotSnrDb', ...
    'csiRsSnrDb','srsRows','rxCorrelationRho'};
assert(all(isfield(config,required)));
assert(all(isfield(resource,{'H_true','prior_gain','omitted_gain', ...
    'home_trp','ueIds'})));
Htrue=double(resource.H_true);
[K,L,Nr,Nt]=size(Htrue);
assert(isequal(size(resource.prior_gain),[K,L]));
assert(isequal(size(resource.omitted_gain),[K,L]));
assert(numel(resource.home_trp)==K && numel(resource.ueIds)==K);
assert(config.pilotSequences>=1 && ...
    config.pilotSequences==floor(config.pilotSequences));
assert(K>=1);
assert(config.pilotMatchedEnergy>0 && isfinite(config.pilotMatchedEnergy));
assert(isfinite(config.pilotSnrDb) && isfinite(config.csiRsSnrDb));
assert(numel(config.srsRows)==2 && all(config.srsRows>=1) && ...
    all(config.srsRows<=Nr) && numel(unique(config.srsRows))==2);
assert(config.rxCorrelationRho>=0 && config.rxCorrelationRho<1);
assert(config.seed>=0 && config.seed==floor(config.seed));
gain=double(resource.prior_gain);
modeled=resource.omitted_gain==0;
gain(~modeled)=0;
assert(all(gain(:)>=0) && all(isfinite(gain(:))));
for k=1:K
    assert(all(Htrue(k,~modeled(k,:),:,:)==0,'all'));
end
Rbase=config.rxCorrelationRho.^abs((1:Nr)'-(1:Nr));
Rrx=complex(zeros(Nr,Nr,K,L));
for l=1:L
    for k=1:K
        Rrx(:,:,k,l)=gain(k,l)*Rbase;
    end
end
homeIndex=sub2ind([K,L],(1:K)',double(resource.home_trp(:)));
referenceGain=median(gain(homeIndex));
assert(referenceGain>0);
% Unit matched-filter energy means each of the two port pilot resources
% carries E=1. This is a configurable convention, not a pilot length.
pilotEnergy=repmat(double(config.pilotMatchedEnergy),K,2);
% A group with fewer users than available pilots uses one distinct pilot
% per user; the unoccupied resources do not affect its observation.
usedSequences=min(config.pilotSequences,K);
userPilot=chasi_assign_user_pilots(gain,usedSequences,config.seed);
pilotId=chasi_expand_port_pilots(userPilot,usedSequences,2);
pilotNoiseVariance=referenceGain/10^(config.pilotSnrDb/10);
% Conditional 4-row CSI-RS noise follows the target's gain-relative SNR
% convention on modeled links. Zero-gain omitted links have zero prior and
% remain zero; their arbitrary positive noise is immaterial.
csiNoiseVariance=max(gain,referenceGain*eps)/10^(config.csiRsSnrDb/10);
stream=RandStream('mt19937ar','Seed',config.seed+1);
Nc=(randn(stream,K,L,Nr,Nt)+1i*randn(stream,K,L,Nr,Nt))/sqrt(2);
Np=(randn(stream,Nt,L,2*usedSequences)+ ...
    1i*randn(stream,Nt,L,2*usedSequences))/sqrt(2);
[Hhat,Ce,crossCe,detail]=chasi_static_measurement_chain( ...
    Htrue,Rrx,csiNoiseVariance,Nc,pilotId,pilotEnergy, ...
    pilotNoiseVariance,config.srsRows,Np);
errorPower=abs(Htrue-Hhat).^2;
truePower=abs(Htrue).^2;
linkError=reshape(sum(errorPower,[3,4]),K,L);
linkTrue=reshape(sum(truePower,[3,4]),K,L);
linkNmse=nan(K,L);
valid=linkTrue>0;
linkNmse(valid)=linkError(valid)./linkTrue(valid);
out=struct('modelId','target_gaussian_static_pilot_lmmse_hypothesis_v1', ...
    'Hhat',Hhat,'Ce',Ce,'crossCe',crossCe,'pilotId',pilotId, ...
    'pilotEnergy',pilotEnergy,'pilotNoiseVariance',pilotNoiseVariance, ...
    'csiNoiseVariance',csiNoiseVariance,'referenceGain',referenceGain, ...
    'Rrx',Rrx,'detail',detail,'linkNmse',linkNmse, ...
    'globalNmse',sum(linkError(:))/sum(linkTrue(:)), ...
    'availableUserPilots',config.pilotSequences, ...
    'usedUserPilots',usedSequences, ...
    'config',config,'ueIds',resource.ueIds(:));
end
