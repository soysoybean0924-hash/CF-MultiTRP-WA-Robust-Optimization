function out=chasi_v6_srs_prior_blocks(jointPrior,targetIds,rxPorts)
% Compact audit rows: every scheduled UE against every SRS sender/TRP.
[rows,cols,L]=size(jointPrior);assert(rows==cols && ...
    mod(rows,rxPorts)==0 && numel(unique(targetIds))==numel(targetIds));
U=rows/rxPorts;targetIds=double(targetIds(:));
assert(all(targetIds>=1 & targetIds<=U & targetIds==fix(targetIds)));
blocks=complex(zeros(rxPorts,rxPorts,numel(targetIds),U,L,'single'));
for l=1:L,for k=1:numel(targetIds)
 ix=(targetIds(k)-1)*rxPorts+(1:rxPorts);
 for u=1:U
  jx=(u-1)*rxPorts+(1:rxPorts);
  blocks(:,:,k,u,l)=single(jointPrior(ix,jx,l));
 end
end,end,end
out=struct('targetUeIds',targetIds,'senderUeIds',(1:U)', ...
    'blocks',blocks,'channelUnit','noise_normalized_channel');
end
