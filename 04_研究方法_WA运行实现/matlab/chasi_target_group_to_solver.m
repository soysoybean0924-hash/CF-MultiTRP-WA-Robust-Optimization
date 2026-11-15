function out=chasi_target_group_to_solver(resource,HdesignDense)
%CHASI_TARGET_GROUP_TO_SOLVER Structural adapter for ONE target RB group.
% Inputs have [groupUE,TRP,DL-Rx,TRP-Tx]. No unit conversion, scheduling,
% rank change, or noise scaling is performed. Caller must verify those.
required={'H_true','servingMask','ueIds','rb_count','resourceWeightNu'};
assert(all(isfield(resource,required)));
truth=resource.H_true;
K=size(truth,1); L=size(truth,2); Nr=size(truth,3); Nt=size(truth,4);
assert(isequal([size(HdesignDense,1),size(HdesignDense,2), ...
    size(HdesignDense,3),size(HdesignDense,4)],[K,L,Nr,Nt]));
assert(isequal(size(resource.servingMask),[L,K]));
assert(numel(resource.ueIds)==K && all(isfinite(truth(:))) && ...
    all(isfinite(HdesignDense(:))));
counts=sum(resource.servingMask,1);
if isempty(counts) || any(counts<1) || any(counts~=counts(1))
    error('chasi:VariableServiceCount', ...
        'Source solver requires one fixed serving-TRP count per group.');
end
S=counts(1);
serviceTRP=zeros(K,S);
Htrue=cell(K,L,1); Hdesign=cell(K,L,1);
for k=1:K
    serviceTRP(k,:)=find(resource.servingMask(:,k)).';
    for l=1:L
        Htrue{k,l,1}=reshape(double(truth(k,l,:,:)),Nr,Nt);
        Hdesign{k,l,1}=reshape(double(HdesignDense(k,l,:,:)),Nr,Nt);
    end
end
out=struct('Htrue',{Htrue},'Hdesign',{Hdesign}, ...
    'serviceTRP',serviceTRP,'ueIds',resource.ueIds(:), ...
    'rbCount',resource.rb_count,'resourceWeightNu',resource.resourceWeightNu, ...
    'normalization','unchanged; caller must align noise, power and RB meaning');
if isfield(resource,'tx_rb_count')
    out.txRbCount=resource.tx_rb_count;
end
end
