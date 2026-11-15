function validate_srs_joint_small()
% Deferred deterministic identities; no fixture is generated during delivery.
K=2; L=1; Nr=4; Nt=2;
R=complex(zeros(Nr,Nr,K,L));
for k=1:K, R(:,:,k,1)=eye(Nr); end
H0=complex(zeros(K,L,Nr,Nt)); innovation=H0;
T=4; noiseVariance=1;
z=complex(zeros(Nt,L,T)); E=ones(K,2);
orth=[1 3;2 4]; shared=[1 3;1 3];
o=srs_joint_snapshot(H0,R,orth,E,noiseVariance,z,0,120,innovation,[1 2]);
s=srs_joint_snapshot(H0,R,shared,E,noiseVariance,z,0,120,innovation,[1 2]);
assert(norm(o.crossUE(:,:,1,2,1),'fro')<1e-12);
assert(norm(s.crossUE(:,:,1,2,1),'fro')>0);
for k=1:K
    assert(min(eig((s.C(:,:,k,1)+s.C(:,:,k,1)')/2))>=-1e-12);
    assert(real(trace(o.C(:,:,k,1)))<=real(trace(s.C(:,:,k,1)))+1e-12);
end
aged=srs_joint_snapshot(H0,R,shared,E,noiseVariance,z,20,120,innovation,[1 2]);
assert(real(trace(aged.C(:,:,1,1)))>=real(trace(s.C(:,:,1,1)))-1e-12);
% K=1 and L=1: MATLAB drops R's trailing singleton axes from size(R).
assert(size(z,1)==Nt && size(z,2)==L && size(z,3)==T);
assert(isequal(size(orth(1,:)),[1,2]) && ...
    isequal(size(E(1,:)),[1,2]) && all(orth(1,:)<=T));
assert(size(H0(1,:,:,:),4)==Nt && ...
    size(innovation(1,:,:,:),4)==Nt);
one=srs_joint_snapshot(H0(1,:,:,:),R(:,:,1,:),orth(1,:),E(1,:), ...
    noiseVariance,z,0,120,innovation(1,:,:,:),[1 2]);
assert(size(one.C,1)==Nr && size(one.C,2)==Nr && ...
    size(one.C,3)==1 && size(one.C,4)==1);
% A single port is outside this explicit 2T4R mapping contract.
try
    srs_joint_snapshot(H0(1,:,:,:),R(:,:,1,:),orth(1,1),E(1,1), ...
        noiseVariance,z,0,120,innovation(1,:,:,:),1);
    error('validate_srs_joint_small:MissingPortRejection', ...
        'Single-port call unexpectedly accepted.');
catch ME
    assert(strcmp(ME.identifier,'srs_joint_snapshot:UnsupportedPortCount'));
end
end
