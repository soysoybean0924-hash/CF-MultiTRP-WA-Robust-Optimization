function layer=chasi_layer1_calibrate(resources,cfg)
% Group holdout empirical squared radius. Source shares are descriptive.
dev=double(cfg.split.development_groups(:)'); test=double(cfg.split.test_groups(:)');
assert(isempty(intersect(dev,test)) && isequal(sort([dev test]),1:numel(resources)));
records=struct('group',{},'ue',{},'trp',{},'source_share',{},'normalized_error',{},'split',{});
train=[]; held=[];
for g=1:numel(resources)
 r=resources(g); K=size(r.H_true,1); L=size(r.H_true,2);
 assert(isequal(size(r.sourceCov),[K,L,4,4,4]) || size(r.sourceCov,5)==4);
 for k=1:K, for l=1:L
  traces=zeros(1,4);
  for q=1:4, traces(q)=real(trace(squeeze(r.sourceCov(k,l,:,:,q)))); end
  assert(all(traces>=-1e-7)); traces=max(traces,0);
  share=traces/max(sum(traces),1e-30);
  h=squeeze(r.H_true(k,l,:,:)-r.H_hat(k,l,:,:));
  observed=mean(abs(h(:)).^2);
  predicted=real(trace(squeeze(r.C(k,l,:,:))))/size(r.C,3);
  value=observed/max(predicted,1e-30);
  row=struct('group',g,'ue',double(r.ueIds(k)),'trp',l,'source_share',share, ...
   'normalized_error',value,'split','test');
  if ismember(g,dev),row.split='development';train(end+1)=value; %#ok<AGROW>
  else,held(end+1)=value; end %#ok<AGROW>
  records(end+1)=row; %#ok<AGROW>
 end, end
end
assert(~isempty(train) && ~isempty(held));
train=sort(train); idx=max(1,ceil(cfg.coverage_target*numel(train)));
layer=struct('radius_squared',train(idx),'records',records, ...
 'development_coverage',mean(train<=train(idx)), ...
 'test_coverage',mean(held<=train(idx)), ...
 'calibration_scope','one scene, held-out resource groups only');
end
