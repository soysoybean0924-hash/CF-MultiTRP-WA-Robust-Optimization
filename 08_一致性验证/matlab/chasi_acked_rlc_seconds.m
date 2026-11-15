function out=chasi_acked_rlc_seconds(ev,ttiMs)
% Completed burst assigned to the queue-empty ACK second. TS 28.554
% excludes both volume and time of the empty-buffer tail TTI.
% Open bursts at the observation end are censored, never zero-filled.
needed={'ackedBits','bufferBeforeBits','bufferAfterBits', ...
    'arrivalBits','attemptedBits','completedSduBits'};
for j=1:numel(needed),assert(isfield(ev,needed{j}),needed{j});end
[K,T]=size(ev.ackedBits);assert(ttiMs>0 && mod(1000,ttiMs)==0);
if isfield(ev,'sduAckBits')
 assert(isequal(ev.sduAckBits,ev.completedSduBits));
end
for j=1:numel(needed)
 a=ev.(needed{j});assert(isequal(size(a),[K,T]) && ...
     all(isfinite(a(:))) && all(a(:)>=0));
end
assert(all(ev.ackedBits(:)<=ev.bufferBeforeBits(:)+ev.arrivalBits(:)+1e-6));
slotsPerSecond=1000/ttiMs;S=ceil(T/slotsPerSecond);
bits=zeros(K,S);rawAck=zeros(K,S);excludedTail=zeros(K,S);
duration=zeros(K,S);
bursts=struct('ue',{},'firstTti',{},'completionTti',{}, ...
    'completionSecond',{},'ackedBits',{},'deliveredAckBits',{}, ...
    'excludedTailAckBits',{},'effectiveSeconds',{}, ...
    'crossSecond',{},'zeroDuration',{});
open=zeros(K,1);neverStarted=zeros(K,1);crossCount=0;zeroCount=0;
for k=1:K
 active=false;first=0;
 for t=1:T
  if ~active && ev.attemptedBits(k,t)>0
   active=true;first=t;
  end
  if active && ev.bufferAfterBits(k,t)==0
   sec=ceil(t/slotsPerSecond);
   payload=sum(ev.completedSduBits(k,first:t-1));
   delivered=sum(ev.completedSduBits(k,first:t));
   elapsed=(t-first)*ttiMs/1000;
   bits(k,sec)=bits(k,sec)+payload;
   rawAck(k,sec)=rawAck(k,sec)+delivered;
   excludedTail(k,sec)=excludedTail(k,sec)+ev.completedSduBits(k,t);
   duration(k,sec)=duration(k,sec)+elapsed;
   crosses=ceil(first/slotsPerSecond)~=sec;
   crossCount=crossCount+crosses;zeroCount=zeroCount+(elapsed==0);
   bursts(end+1)=struct('ue',k,'firstTti',first, ...
       'completionTti',t,'completionSecond',sec, ...
       'ackedBits',payload,'deliveredAckBits',delivered, ...
       'excludedTailAckBits',ev.completedSduBits(k,t), ...
       'effectiveSeconds',elapsed, ...
       'crossSecond',crosses,'zeroDuration',elapsed==0); %#ok<AGROW>
   active=false;
  end
 end
 if active || ev.bufferAfterBits(k,T)>0
  open(k)=1;
  neverStarted(k)=~active;
 end
end
valid=duration>0;rate=nan(K,S);
rate(valid)=bits(valid)./duration(valid)/1e6;
out=struct('ackedBits',bits,'effectiveSeconds',duration, ...
    'deliveredAckBits',rawAck,'excludedTailAckBits',excludedTail, ...
    'fragmentAckBits',ev.ackedBits, ...
    'rateMbps',rate,'valid',valid,'missingCount',nnz(~valid), ...
    'windowSeconds',S,'crossSecondBursts',crossCount, ...
    'zeroDurationBursts',zeroCount,'openBurstCount',sum(open), ...
    'neverStartedCensoredCount',sum(neverStarted), ...
    'bursts',bursts, ...
    'rule','completed burst to ACK/queue-empty second; volume and duration exclude last TTI; open bursts censored');
end
