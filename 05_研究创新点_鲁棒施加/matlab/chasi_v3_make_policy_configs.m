function paths=chasi_v3_make_policy_configs(root,outDir,input,developmentGroups)
% Lock development-only normalization and write reproducible ablations.
startup;
if nargin<3 || isempty(input)
    input=fullfile(root,'03_研究场景构造','inputs', ...
        'srs_overlap_jakes_v3_main','input.mat');
end
if nargin<4 || isempty(developmentGroups),developmentGroups=1:8;end
x=load(input,'data');d=x.data;
cal=chasi_v3_calibrate_strength(d,developmentGroups);
if ~isfolder(outDir),mkdir(outDir);end
allowed=logical([0 1 0;1 0 0;0 1 0;1 1 1]);
names={'nominal','physical','uniform','source_all', ...
    'source_service','source_beam','source_power', ...
    'source_residual','source_pollution','source_noise','source_age'};
paths=cell(size(names));
for i=1:numel(names)
 p=struct('enabled',false(4,3),'strength',.1*ones(4,3), ...
     'strengthMode','link_source','calibration',cal, ...
     'powerReserveFraction',.1,'serviceMargin',.01);
 s=struct('mode','source','uniformExtra',.1,'policy',p, ...
     'powerControl',struct('powerPrice',.01,'iterations',30,'step',.1));
 switch names{i}
  case 'nominal',s.mode='nominal';
  case 'physical',s.mode='physical';
  case 'uniform',s.mode='uniform';
  case 'source_all',p.enabled=allowed;
  case 'source_service',p.enabled(:,1)=allowed(:,1);
  case 'source_beam',p.enabled(:,2)=allowed(:,2);
  case 'source_power',p.enabled(:,3)=allowed(:,3);
  case 'source_residual',p.enabled(1,:)=allowed(1,:);
  case 'source_pollution',p.enabled(2,:)=allowed(2,:);
  case 'source_noise',p.enabled(3,:)=allowed(3,:);
  case 'source_age',p.enabled(4,:)=allowed(4,:);
 end
 s.policy=p;
 path=fullfile(outDir,[names{i} '.json']);
 assert(~isfile(path),'Refusing to overwrite locked policy config');
 fid=fopen(path,'w');assert(fid>=0);clean=onCleanup(@()fclose(fid));
 fwrite(fid,jsonencode(s,PrettyPrint=true),'char');clear clean
 paths{i}=path;
end
end
