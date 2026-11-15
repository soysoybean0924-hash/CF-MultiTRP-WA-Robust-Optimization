function report=acceptance_dependency_paths(root,outPath)
% Assert active functions resolve to this D target, never to source copies.
names={'chasi_config','chasi_generate_scene','srs_joint_snapshot', ...
 'chasi_load_frozen_binary','chasi_matrix_wmmse','run_srs_joint_batch', ...
 'chasi_layer1_calibrate','chasi_layer_protection','chasi_run_layered_small'};
rows=struct('name',{},'resolved_path',{},'in_target',{});
for j=1:numel(names)
 p=which(names{j});assert(~isempty(p),['Missing ' names{j}]);
 rows(j)=struct('name',names{j},'resolved_path',p, ...
  'in_target',startsWith(p,[root filesep]));
end
report=struct('all_in_target',all([rows.in_target]), ...
 'history_09_on_matlab_path',contains(path,'09_历史版本与待确认项'), ...
 'functions',rows);
assert(report.all_in_target && ~report.history_09_on_matlab_path);
fid=fopen(outPath,'wt');assert(fid>=0);
fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));fclose(fid);
end
