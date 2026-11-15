function [sigma,parts]=arvidson_srs_sigma(targetAmplitude,interfererAmplitude,coeff)
%ARVIDSON_SRS_SIGMA Thesis Eq. 23, per RB received-amplitude model.
% Interferers include every SRS transmitter sharing this comb, scheduled or not.
% Coefficients must be calibrated for the channel/sequence/bandwidth setting.
% Amplitudes and sigma have identical channel-amplitude units. This function
% does not add a second noise draw to an existing joint SRS observation.
arguments
 targetAmplitude (:,1) double {mustBeNonnegative}
 interfererAmplitude (:,:) double {mustBeNonnegative}
 coeff struct
end
assert(isfield(coeff,'alpha') && isfield(coeff,'beta') && isfield(coeff,'gamma'));
assert(isscalar(coeff.alpha) && isscalar(coeff.gamma));
assert(numel(coeff.beta)==size(interfererAmplitude,2));
assert(all(isfinite([coeff.alpha coeff.beta(:)' coeff.gamma])));
assert(all([coeff.alpha coeff.beta(:)' coeff.gamma]>=0));
parts=struct('truncation',coeff.alpha*targetAmplitude,...
 'pilot',interfererAmplitude*coeff.beta(:),...
 'background',coeff.gamma*ones(size(targetAmplitude)));
sigma=parts.truncation+parts.pilot+parts.background;
% Eq. 23 adds empirical standard-deviation contributions, not variances.
end
