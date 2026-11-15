function result=arvidson_eq23_fit_trials(trials,split)
% Fit Eq.23 to conditional complex-channel RMSE, never single-draw |error|.
% trials: one row per independent sounding realization. Required columns:
% conditionId, truth, estimate, targetRx, interfererRx [N-by-J], channelType,
% rb, comb, cyclicShift. Complex truth/estimate; Rx entries are amplitudes.
% split: disjoint condition IDs train and test; no test-conditioned fitting.
required={'conditionId','truth','estimate','targetRx','interfererRx', ...
    'channelType','rb','comb','cyclicShift'};
for j=1:numel(required),assert(isfield(trials,required{j}),required{j});end
n=numel(trials.conditionId);assert(n>=4 && numel(trials.truth)==n);
assert(numel(trials.estimate)==n && numel(trials.targetRx)==n);
assert(size(trials.interfererRx,1)==n && all(trials.targetRx(:)>=0));
assert(all(trials.interfererRx(:)>=0));
ids=unique(trials.conditionId(:));train=split.train(:);test=split.test(:);
assert(~isempty(train) && ~isempty(test) && ...
    isempty(intersect(train,test)) && all(ismember([train;test],ids)));
channelTypes=string(trials.channelType(:));
assert(numel(channelTypes)==n && numel(unique(channelTypes))==1, ...
    'Fit one channel type at a time; original coefficients are channel-specific');
assert(numel(trials.rb)==n && numel(trials.comb)==n && ...
    numel(trials.cyclicShift)==n && ...
    all(isfinite(double(trials.rb(:)))));
mode='project_extension';
if isfield(split,'configurationMode'),mode=char(split.configurationMode);end
assert(ismember(mode,{'project_extension','original_lte_4rb_inputs'}));
if strcmp(mode,'original_lte_4rb_inputs')
    assert(all(trials.rb(:)==4), ...
        'Original LTE mode requires 4 RB independent estimator trials');
end
J=size(trials.interfererRx,2);X=zeros(numel(ids),J+2);y=zeros(numel(ids),1);
counts=zeros(numel(ids),1);
for i=1:numel(ids)
    rows=trials.conditionId(:)==ids(i);counts(i)=nnz(rows);
    assert(counts(i)>=2,'At least two independent trials per condition required');
    y(i)=sqrt(mean(abs(trials.estimate(rows)-trials.truth(rows)).^2));
    X(i,:)=[mean(trials.targetRx(rows)), ...
        mean(trials.interfererRx(rows,:),1),1];
end
tr=ismember(ids,train);te=ismember(ids,test);
assert(rank(X(tr,:))==size(X,2),'Training design is rank deficient');
[coef,solver]=nonnegativeFit(X(tr,:),y(tr));pred=X*coef;
result=struct('model','Arvidson_Eq23_conditional_RMSE_fit', ...
    'configurationMode',mode,'channelType',char(channelTypes(1)), ...
    'rbValues',unique(trials.rb(:))', ...
    'coefficient',coef,'fitSolver',solver, ...
    'conditionId',ids,'sampleCount',counts, ...
    'observedSigma',y,'predictedSigma',pred, ...
    'trainConditionIds',train,'testConditionIds',test, ...
    'trainMae',mean(abs(pred(tr)-y(tr))), ...
    'testMae',mean(abs(pred(te)-y(te))), ...
    'testRelativeMae',mean(abs(pred(te)-y(te)))/max(mean(y(te)),eps), ...
    'unit','received complex channel amplitude standard deviation', ...
    'originalConfigurationReproduced',false, ...
    'limitation','Requires independent LTE link-simulator trials with exact EPA/ETU/AWGN and original estimator; no such generator bundled');
end

function [c,info]=nonnegativeFit(A,y)
% Base-MATLAB projected NNLS, keeping fit independent of Optimization Toolbox.
Q=A'*A;b=A'*y;step=1/max(norm(Q,2),eps);
c=zeros(size(A,2),1);maxIterations=100000;tolerance=1e-9;
for it=1:maxIterations
    next=max(0,c-step*(Q*c-b));
    if norm(next-c)<=tolerance*(1+norm(c))
        c=next;break
    end
    c=next;
end
gradient=Q*c-b;
kkt=norm(c-max(0,c-gradient));
assert(it<maxIterations && kkt<=1e-6*(1+norm(b)), ...
    'Nonnegative Eq.23 fit did not meet KKT tolerance');
info=struct('algorithm','base_matlab_projected_nnls', ...
    'iterations',it,'kktResidual',kkt, ...
    'tolerance',tolerance,'maxIterations',maxIterations);
end
