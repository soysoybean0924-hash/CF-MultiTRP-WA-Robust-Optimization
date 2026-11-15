function test_csi_rs_gaussian_prior()
% Explicit full-row observation fixture; this is a CHASI adaptation.
H=complex(zeros(1,1,4,2));
H(1,1,:,1)=reshape([2;4;6;8],1,1,4,1);
H(1,1,:,2)=reshape([1;3;5;7],1,1,4,1);
R=reshape(eye(4),4,4,1,1);
noise=1; Nc=complex(zeros(size(H)));
[mean,C,observation]=chasi_csi_rs_gaussian_prior(H,R,noise,Nc);
assert(norm(mean(:)-H(:)/2)<1e-12);
assert(norm(reshape(C,4,4)-eye(4)/2,'fro')<1e-12);
assert(norm(observation(:)-H(:))<1e-12);
% One random non-diagonal covariance must match the target precision form.
R0=0.35.^abs((1:4)'-(1:4)); R=reshape(R0,4,4,1,1);
[mean,C]=chasi_csi_rs_gaussian_prior(H,R,noise,Nc);
expectedC=(R0\eye(4)+eye(4)/noise)\eye(4);
expectedMean=expectedC*reshape(H,4,2)/noise;
assert(norm(reshape(C,4,4)-expectedC,'fro')<1e-12);
assert(norm(reshape(mean,4,2)-expectedMean,'fro')<1e-12);
fprintf('test_csi_rs_gaussian_prior PASS\n');
end
