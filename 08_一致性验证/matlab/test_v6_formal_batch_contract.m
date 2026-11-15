function test_v6_formal_batch_contract(root)
% Preflight contract; formal results still require the full online run.
startup;
assert(exist('run_v6_formal_batch','file')==2);
c=chasi_config(struct('usersPerCell',10,'prbUtilization',.5));
assert(c.sites==7 && c.sectorsPerSite==3 && c.trps==21);
assert(c.users==210 && c.siteSpacingM==300);
assert(c.txPorts==64 && c.ueTxPorts==2 && c.ueRxPorts==4);
assert(c.bandwidthHz==100e6 && c.scsHz==30e3 && c.rbTotal==273);
assert(c.srsPeriodTti==340 && c.srsHopIntervalTti==20 && ...
    isequal(c.srsRbCounts,[repmat(16,1,16),17]));
assert(isfile(fullfile(root,'03_研究场景构造','config', ...
    'srs_csi_joint_v6_formal.json')));
end
