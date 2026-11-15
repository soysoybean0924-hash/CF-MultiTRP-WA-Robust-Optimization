function validate_small(root)
% Server only. Separate output; never touches formal checkpoints.
if nargin<1, root=fileparts(mfilename('fullpath')); end
addpath(fullfile(root,'pilot'),fullfile(root,'matlab','chasi'));
folder=fullfile(root,'validation');
if ~isfolder(folder), mkdir(folder); end
diary(fullfile(folder,'small_validation.log'));
guard=onCleanup(@()diary('off')); %#ok<NASGU>
test_pilot_assignment;
test_csi_rs_gaussian_prior;
test_static_measurement_chain;
test_static_pilot_lmmse;
% test_target_aged_small requires a generated formal input and is run
% separately after run_campaign('run') has prepared it.
% These tests cover covariance PSD, empirical/theoretical error, shared
% pilot cross covariance, zero age and no contamination within their models.
fprintf('SERVER_SMALL_VALIDATION_PASS\n');
end
