function report=validate_v3_source_policy_small(root)
% Bounded, single-group policy contract check; not a performance result.
if nargin<1,root=fileparts(fileparts(fileparts(mfilename('fullpath'))));end
startup;
input=fullfile(root,'03_研究场景构造','inputs','srs_overlap_jakes_v3_main','input.mat');
assert(strcmp(chasi_sha256_file(input),strtrim(fileread([input '.sha256']))));
x=load(input,'data'); data=x.data;
p=struct('enabled',false(4,3),'strength',ones(4,3), ...
    'powerReserveFraction',0.1,'serviceMargin',0);
[off,a0]=chasi_v3_source_policy(data,1,p);
base=data.resources(1);
assert(isequaln(off.C,base.C) && isequaln(off.servingMask,base.servingMask));
assert(norm(off.policyCov(:))==0 && off.powerWattPerTrp==data.config.powerWattPerTrp);
p.enabled(3,2)=true;p.strength(3,2)=0.2;
[beam,a1]=chasi_v3_source_policy(data,1,p);
assert(isequaln(beam.C,base.C) && norm(beam.policyCov(:))>0);
assert(isequaln(beam.servingMask,base.servingMask));
p.enabled(:)=false;p.enabled(2,1)=true;p.serviceMargin=0;
[service,~]=chasi_v3_source_policy(data,1,p);
assert(isequal(sum(service.servingMask,1),sum(base.servingMask,1)));
assert(isequaln(service.C,base.C) && norm(service.policyCov(:))==0);
p.enabled(:)=false;p.enabled(4,3)=true;
[power,a2]=chasi_v3_source_policy(data,1,p);
assert(all(power.powerWattPerTrp<=data.config.powerWattPerTrp) && ...
    all(power.powerWattPerTrp>0));
report=struct('modelId',data.modelId,'group',1,'allOffIdentical',true, ...
    'beamOnlyPhysicalCIdentical',true,'serviceCardinalityPreserved',true, ...
    'powerFeasible',true, ...
    'sourceMismatch',[a0.sourceMismatch,a1.sourceMismatch,a2.sourceMismatch]);
disp(report);
end
