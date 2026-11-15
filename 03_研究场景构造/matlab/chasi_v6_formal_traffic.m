function traffic=chasi_v6_formal_traffic(users,ttiCount,seed,busyTarget, ...
    payloadBits,tailPayloadFraction)
% Reuse chasi_traffic's burst intervals, not its static-rate proxy.
% At each burst start enqueue one SDU per planned busy TTI. The packet
% queue/ACK model decides the actual end and may extend past the plan.
assert(users>=1 && ttiCount>=1 && busyTarget>0 && busyTarget<1);
if nargin<5,payloadBits=12000;end
if nargin<6,tailPayloadFraction=.3;end
assert(payloadBits>0 && payloadBits==fix(payloadBits) && ...
    tailPayloadFraction>0 && tailPayloadFraction<1);
[busy,tail,terminal]=chasi_traffic(users,ttiCount,seed, ...
    'target_busy_fraction',busyTarget);
maxSdus=30;arrival=zeros(users,ttiCount,maxSdus,'uint32');
rows=zeros(0,8);totalBits=0;tailBits=0;burstCount=0;
for u=1:users
 t=1;
 while t<=ttiCount
  if ~busy(u,t),t=t+1;continue;end
  start=t;finish=find(terminal(u,t:end),1);
  assert(~isempty(finish),'Traffic burst has no terminal marker');
  finish=t+finish-1;n=finish-start+1;
  assert(n<=maxSdus && n>=1);
  tailCount=nnz(tail(u,start:finish));headCount=n-tailCount;
  assert(headCount>=1 && tailCount>=1);
  total=payloadBits*n;
  tailTotal=round(tailPayloadFraction*total);
  headTotal=total-tailTotal;
  sizes=[repmat(floor(headTotal/headCount),1,headCount), ...
      repmat(floor(tailTotal/tailCount),1,tailCount)];
  sizes(headCount)=sizes(headCount)+headTotal-sum(sizes(1:headCount));
  sizes(end)=sizes(end)+tailTotal-sum(sizes(headCount+1:end));
  arrival(u,start,1:n)=uint32(sizes);
  burstCount=burstCount+1;totalBits=totalBits+total;
  tailBits=tailBits+sum(sizes(headCount+1:end));
  rows(end+1,:)=[u,start,finish,n,headCount,tailCount,total, ...
      sum(sizes(headCount+1:end))]; %#ok<AGROW>
  t=finish+1;
 end
end
assert(nnz(arrival)>0 && totalBits==sum(double(arrival),'all'));
traffic=struct('arrivalSduBits',arrival,'plannedBusy',busy, ...
    'plannedTail',tail,'plannedTerminal',terminal, ...
    'burstRows',rows,'model','chasi_traffic_target_busy_fraction_batched_sdu_v6', ...
    'seed',seed,'targetBusyFraction',busyTarget, ...
    'actualPlannedBusyFraction',mean(busy,'all'), ...
    'actualTailTtiFraction',nnz(tail)/nnz(busy), ...
    'actualTailPayloadFraction',tailBits/totalBits, ...
    'burstCount',burstCount,'sduCount',nnz(arrival), ...
    'offeredPayloadBits',totalBits,'payloadBitsPerPlannedBusyTti',payloadBits, ...
    'scope','offered UE busy fraction is a proxy for unspecified 30% network load');
end
