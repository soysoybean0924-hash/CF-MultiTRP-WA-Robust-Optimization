function test_srs_joint_branch_mapping(root)
% Deferred source guard for the actual formal runner branch mapping.
if nargin<1, root=fileparts(fileparts(mfilename('fullpath'))); end
source=fileread(fullfile(root,'pilot','run_srs_joint_batch.m'));
names={'ideal_sumse','nominal_sumse','robust_sumse', ...
    'ideal_pf','nominal_pf','robust_pf'};
for i=1:numel(names), assert(contains(source,['case ''' names{i} ''''])); end
block=extractBetween(source,"case 'ideal_pf'","case 'nominal_pf'");
assert(numel(block)==1 && contains(block,'resources(g).H=resources(g).H_true'));
assert(contains(block,"nativeCfg.branch='matrix_pf_nominal'"));
block=extractBetween(source,"case 'nominal_pf'","case 'robust_pf'");
assert(numel(block)==1 && ~contains(block,'H_true'));
end
