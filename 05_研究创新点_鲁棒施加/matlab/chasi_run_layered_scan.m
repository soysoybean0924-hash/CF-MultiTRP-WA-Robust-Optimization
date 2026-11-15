function results=chasi_run_layered_scan(inputRoot,configPath,kind,outRoot,maxGroups,maxOuter)
% Explicit opt-in scan; each case has a distinct identity-checked directory.
assert(~isfolder(outRoot),'Choose a new scan output root');mkdir(outRoot);
cfg=jsondecode(fileread(configPath));
switch kind
 case 'strength'
  tags={'x0','x05','x1','x15','x2'};
  levels=cfg.fixed_level*[0 .5 1 1.5 2];
  modes=repmat({'mapped'},1,5);abl=repmat({'none'},1,5);
  trig=repmat({'proposed'},1,5);caps=nan(1,5);
 case 'ablation'
  tags={'none','beam','service','power'};
  levels=cfg.fixed_level*ones(1,4);modes=repmat({'mapped'},1,4);
  abl=tags;trig=repmat({'proposed'},1,4);caps=nan(1,4);
 case 'trigger'
  tags={'off','fixed','proposed'};
  levels=cfg.fixed_level*ones(1,3);modes=repmat({'mapped'},1,3);
  abl=repmat({'none'},1,3);trig=tags;caps=nan(1,3);
 case 'power_cap'
  values=cfg.power_caps_w(:)';tags=arrayfun(@(x)sprintf('%gW',x),values,'UniformOutput',false);
  caps=repelem(values,4);tags=cellfun(@(a,b)[a '_' b],repelem(tags,4), ...
   repmat({'uniform','variable','fixed','mapped'},1,numel(values)),'UniformOutput',false);
  modes=repmat({'uniform','variable','fixed','mapped'},1,numel(values));
  levels=cfg.fixed_level*ones(1,numel(caps));
  abl=repmat({'none'},1,numel(caps));trig=repmat({'proposed'},1,numel(caps));
 otherwise,error('Unknown scan kind');
end
results=cell(1,numel(tags));
for i=1:numel(tags)
 cap=[];if ~isnan(caps(i)),cap=caps(i);end
 results{i}=chasi_run_layered_small(inputRoot,configPath,modes{i}, ...
  maxGroups,maxOuter,fullfile(outRoot,tags{i}),abl{i},trig{i},cap,levels(i));
end
end
