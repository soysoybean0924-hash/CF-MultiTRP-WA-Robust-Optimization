function pilotId=chasi_expand_port_pilots(userPilot,tau,numPorts)
%CHASI_EXPAND_PORT_PILOTS Orthogonal SRS ports reusing user pilots.
% Pilot t on port q maps to observation index t+(q-1)*tau. Hence each
% user needs numPorts mutually orthogonal port resources and the observation
% has tau*numPorts dimensions. This is a project resource assumption.
userPilot=userPilot(:);
assert(tau>=1 && tau==floor(tau) && numPorts>=1 && ...
    numPorts==floor(numPorts));
assert(all(userPilot>=1) && all(userPilot<=tau) && ...
    all(userPilot==floor(userPilot)));
pilotId=userPilot+tau*(0:numPorts-1);
end
