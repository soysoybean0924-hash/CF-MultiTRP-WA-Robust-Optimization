function test_static_measurement_chain()
% Four-row CSI-RS first, one observed SRS port second, all noises explicit.
H=complex(zeros(1,1,4,1)); H(1,1,1,1)=2; H(1,1,2,1)=4;
R=reshape(eye(4),4,4,1,1);
Nc=complex(zeros(size(H))); Np=complex(zeros(1,1,1));
[Hhat,Ce,crossCe,detail]=chasi_static_measurement_chain( ...
    H,R,1,Nc,1,1,1,1,Np);
assert(abs(Hhat(1,1,1,1)-4/3)<1e-12);
assert(abs(Hhat(1,1,2,1)-2)<1e-12);
assert(abs(Ce(1,1,1,1)-1/3)<1e-12);
assert(abs(Ce(2,2,1,1)-1/2)<1e-12);
assert(abs(crossCe(1,1,1,1,1)-Ce(1,1,1,1))<1e-12);
assert(abs(detail.csiPriorMean(1,1,1,1)-1)<1e-12);
fprintf('test_static_measurement_chain PASS\n');
end
