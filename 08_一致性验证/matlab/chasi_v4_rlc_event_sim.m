function ev=chasi_v4_rlc_event_sim(arrivalSduBits,capacityBits,bler,seed,commonDraw)
% One or more SDU arrivals per UE/TTI ([UE,TTI,arrival slot]); one
% transport attempt per TTI. SDUs may span TTIs; a failed fragment is
% retransmitted before sending new bytes from the same SDU.
% ACK is emitted only on a successful attempted transmission. No PHY proxy
% is inferred here: capacityBits and BLER must come from an online channel.
[K,T,N]=size(arrivalSduBits);
assert(isequal(size(capacityBits),[K,T],size(bler)));
assert(all(arrivalSduBits(:)>=0) && all(capacityBits(:)>=0));
assert(all(bler(:)>=0 & bler(:)<=1));
stream=RandStream('mt19937ar','Seed',seed);
if nargin<5 || isempty(commonDraw),commonDraw=rand(stream,K,T);end
assert(isequal(size(commonDraw),[K,T]) && ...
    all(commonDraw(:)>=0 & commonDraw(:)<1));
ev=struct('ackedBits',zeros(K,T),'bufferBeforeBits',zeros(K,T), ...
    'bufferAfterBits',zeros(K,T),'arrivalBits',sum(arrivalSduBits,3), ...
    'arrivalSduCount',sum(arrivalSduBits>0,3), ...
    'queueSduCountBefore',zeros(K,T), ...
    'queueSduCountAfter',zeros(K,T), ...
    'attemptedBits',zeros(K,T),'retransmission',false(K,T), ...
    'completedSduBits',zeros(K,T),'completedSduCount',zeros(K,T), ...
    'success',false(K,T),'failure',false(K,T), ...
    'attemptDraw',nan(K,T), ...
    'attemptBler',nan(K,T),'attemptNumber',zeros(K,T), ...
    'sinrDb',nan(K,T),'mcsIndex',nan(K,T), ...
    'attemptPowerW',nan(K,T), ...
    'ackModel','fixed_bler_proxy_unvalidated');
for k=1:K
    queue=zeros(1,0);original=zeros(1,0);failedBits=0;attempts=0;
    for t=1:T
        ev.bufferBeforeBits(k,t)=sum(queue);
        ev.queueSduCountBefore(k,t)=numel(queue);
        for n=1:N
            if arrivalSduBits(k,t,n)>0
                queue(end+1)=arrivalSduBits(k,t,n); %#ok<AGROW>
                original(end+1)=arrivalSduBits(k,t,n); %#ok<AGROW>
            end
        end
        available=floor(capacityBits(k,t));
        if ~isempty(queue) && available>0 && ...
                (failedBits==0 || available>=failedBits)
            if failedBits>0,sendBits=failedBits;
            else,sendBits=min(queue(1),available);end
            ev.attemptedBits(k,t)=sendBits;
            ev.retransmission(k,t)=failedBits>0;
            attempts=attempts+1;
            ev.attemptNumber(k,t)=attempts;
            ev.attemptDraw(k,t)=commonDraw(k,t);
            ev.attemptBler(k,t)=bler(k,t);
            if commonDraw(k,t)>=bler(k,t)
                ev.ackedBits(k,t)=sendBits;ev.success(k,t)=true;
                queue(1)=queue(1)-sendBits;failedBits=0;attempts=0;
                if queue(1)==0
                    ev.completedSduBits(k,t)=original(1);
                    ev.completedSduCount(k,t)=1;
                    queue(1)=[];original(1)=[];
                end
            else
                failedBits=sendBits;
                ev.failure(k,t)=true;
            end
        end
        ev.bufferAfterBits(k,t)=sum(queue);
        ev.queueSduCountAfter(k,t)=numel(queue);
    end
end
ev.fragmentAckBits=ev.ackedBits;
ev.sduAckBits=ev.completedSduBits;
end
