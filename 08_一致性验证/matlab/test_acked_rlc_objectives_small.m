function test_acked_rlc_objectives_small()
e=struct('ackedBits',zeros(2,5),'bufferBeforeBits',zeros(2,5), ...
    'bufferAfterBits',zeros(2,5),'arrivalBits',zeros(2,5), ...
    'attemptedBits',zeros(2,5),'completedSduBits',zeros(2,5));
e.bufferBeforeBits(1,2:4)=[300 200 100];
e.bufferAfterBits(1,2:4)=[200 100 0];
e.ackedBits(1,2:4)=100;
e.bufferBeforeBits(2,2:3)=[200 100];
e.bufferAfterBits(2,2:3)=[100 0];e.ackedBits(2,2:3)=100;
e.attemptedBits=e.ackedBits;e.completedSduBits=e.ackedBits;
a=chasi_acked_rlc_seconds(e,1);
assert(a.ackedBits(1)==200 && a.effectiveSeconds(1)==.002);
assert(a.ackedBits(2)==100 && a.effectiveSeconds(2)==.001);
b=struct('rateMbps',[1;2],'valid',true(2,1),'totalAirPowerW',10);
i=b;i.rateMbps=[3;4];
m=b;m.rateMbps=[2.7;3.8];m.totalAirPowerW=6;
d=chasi_two_objective_decision(b,i,m,struct());
assert(strcmp(d.edgeStatus,'pass') && strcmp(d.powerStatus,'pass'));
i.rateMbps=b.rateMbps;z=chasi_two_objective_decision(b,i,m,struct());
assert(strcmp(z.edgeStatus,'nonpositive_denominator'));
fprintf('ACK_RLC bits=%s duration=%s objectives=%s/%s\n', ...
    mat2str(a.ackedBits(:)'),mat2str(a.effectiveSeconds(:)'), ...
    d.edgeStatus,d.powerStatus);
end
