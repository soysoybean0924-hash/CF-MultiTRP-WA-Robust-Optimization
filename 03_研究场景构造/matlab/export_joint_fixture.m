function export_joint_fixture(fixturePath,resultPath,h0,r,pilot_id,energy,noise,z,age_ms,tc_ms,innovation)
% Export exact random arrays and result for Python; caller supplies arrays.
assert(~isfile(fixturePath) && ~isfile(resultPath));
save(fixturePath,'h0','r','pilot_id','energy','noise','z', ...
    'age_ms','tc_ms','innovation','-v7');
out=srs_joint_snapshot(h0,r,pilot_id,energy,noise,z,age_ms,tc_ms,innovation,[1,2]);
H_true=out.H_true; H_hat=out.H_hat; C=out.C;
crossUE=out.crossUE; Y=out.Y;
save(resultPath,'H_true','H_hat','C','crossUE','Y','-v7');
end
