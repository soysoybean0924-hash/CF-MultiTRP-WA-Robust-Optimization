function report=run_arvidson_public_4rb_approx(seed,trialsPerCondition)
% Public-conditions approximation: NOT the thesis proprietary LTE estimator.
% 4 RB, two UEs, 23-length ZC base repeated to 24 comb tones; AWGN and
% illustrative EPA/ETU delay-power profiles. Complex RMSE is per condition.
if nargin<1,seed=441;end
if nargin<2,trialsPerCondition=200;end
stream=RandStream('mt19937ar','Seed',seed);
N=24;n=(0:N-1)';m=mod(n,23);
base=exp(-1i*pi*m.*(m+1)/23);
types={'AWGN','EPA_approx','ETU_approx'};
profiles={[1],[0 -1 -2 -3 -8],[0 -1 -1 -3 -5 -7 -9 -10]};
conditions=12;rows=conditions*trialsPerCondition;
conditionId=repelem((1:conditions)',trialsPerCondition);
truth=complex(zeros(rows,1));estimate=truth;targetRx=zeros(rows,1);
interfererRx=zeros(rows,1);channelType=strings(rows,1);
shift=zeros(rows,1);
for c=1:conditions
 typ=1+mod(c-1,3);power1=10^((mod(c-1,4)-2)/10);
 power2=10^((mod(3*c,7)-3)/10);
 s=mod(c,8);pilot1=base;pilot2=base.*exp(2i*pi*s*n/8);
 for rep=1:trialsPerCondition
  row=(c-1)*trialsPerCondition+rep;
  h1=drawChannel(stream,profiles{typ},N);
  h2=drawChannel(stream,profiles{typ},N);
  noise=.1*(randn(stream,N,1)+1i*randn(stream,N,1))/sqrt(2);
  y=sqrt(power1)*pilot1.*h1+sqrt(power2)*pilot2.*h2+noise;
  raw=conj(pilot1).*y/sqrt(power1);
  taps=ifft(raw);taps(5:end)=0;hEst=fft(taps);
  truth(row)=h1(1);estimate(row)=hEst(1);
  targetRx(row)=abs(sqrt(power1)*h1(1));
  interfererRx(row)=abs(sqrt(power2)*h2(1));
  channelType(row)=types{typ};shift(row)=s;
 end
end
trials=struct('conditionId',conditionId,'truth',truth, ...
    'estimate',estimate,'targetRx',targetRx, ...
    'interfererRx',interfererRx,'channelType',channelType, ...
    'rb',4*ones(rows,1),'comb',zeros(rows,1),'cyclicShift',shift);
fit=arvidson_eq23_fit_trials(trials,struct('train',(1:8)','test',(9:12)'));
report=struct('fit',fit,'seed',seed,'trialsPerCondition',trialsPerCondition, ...
    'originalReproduction',false,'classification','public_conditions_approximation', ...
    'missing','Thesis exact link simulator, channel estimator and trial traces');
fprintf('ARVIDSON_4RB_APPROX testMae=%.6g relative=%.6g\n', ...
    fit.testMae,fit.testRelativeMae);
end

function H=drawChannel(stream,db,N)
p=10.^(db(:)/10);p=p/sum(p);h=complex(zeros(N,1));
for t=1:numel(p)
 h(t)=sqrt(p(t)/2)*(randn(stream)+1i*randn(stream));
end
H=fft(h);
end
