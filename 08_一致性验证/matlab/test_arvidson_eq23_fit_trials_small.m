function result=test_arvidson_eq23_fit_trials_small()
% Synthetic contract test for conditional RMSE fitting, not LTE reproduction.
rng(337,'twister');conditions=12;repeats=400;n=conditions*repeats;
id=repelem((1:conditions)',repeats);t=zeros(n,1);j=zeros(n,1);
for c=1:conditions
    rows=id==c;t(rows)=.2+.1*c;j(rows)=.1+mod(c*5,11)*.13;
end
sigma=.04+.12*t+.20*j;
err=sigma.*(randn(n,1)+1i*randn(n,1))/sqrt(2);
trials=struct('conditionId',id,'truth',complex(ones(n,1)), ...
    'estimate',complex(ones(n,1))+err,'targetRx',t, ...
    'interfererRx',j,'channelType',repmat("AWGN",n,1), ...
    'rb',4*ones(n,1),'comb',zeros(n,1),'cyclicShift',zeros(n,1));
result=arvidson_eq23_fit_trials(trials,struct('train',(1:8)','test',(9:12)'));
assert(result.testRelativeMae<.1 && all(result.coefficient>=0));
fprintf('EQ23_CONDITIONAL synthetic_test_rel_mae=%.6g\n',result.testRelativeMae);
end
