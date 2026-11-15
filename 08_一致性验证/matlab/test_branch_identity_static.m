function test_branch_identity_static()
% Regression source-level guard for the formal solve entry, pending MATLAB run.
source=fileread(fullfile(fileparts(mfilename('fullpath')), ...
    'chasi_run_target_aged_batch.m'));
assert(contains(source,"case 'ideal_pf'"));
block=extractBetween(source,"case 'ideal_pf'","case 'pilot_pf'");
assert(numel(block)==1 && contains(block,"resources(g).H=resources(g).H_true"));
assert(contains(block,"matrix_pf_nominal"));
end
