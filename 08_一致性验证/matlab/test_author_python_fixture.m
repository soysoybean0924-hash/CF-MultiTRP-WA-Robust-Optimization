function test_author_python_fixture()
% The input and Python expected output are generated from main_DNN.py:270-298.
data=load(fullfile('reproduction','chasi_pilot_mmse','inputs', ...
    'author_python_isotropic_fixture.mat'));
[M,~,K,L]=size(data.H); T=double(data.tau_p);
H=complex(zeros(K,L,1,M)); R=complex(zeros(1,1,K,L));
Np=complex(zeros(M,L,T));
for l=1:L
    for k=1:K
        H(k,l,1,:)=reshape(data.H(:,1,k,l),1,1,1,M);
        R(1,1,k,l)=data.CorrR(1,1,k,l);
    end
    for t=1:T
        Np(:,l,t)=data.Np(:,1,l,t);
    end
end
[estimate,Ce]=chasi_pilot_lmmse_static(H,R, ...
    double(data.pilotIndex(:))+1,repmat(double(data.p)*T,K,1),1,1,Np);
expected=complex(zeros(K,L,1,M));
for l=1:L
    for k=1:K
        expected(k,l,1,:)=reshape(data.HhatPython(:,1,k,l),1,1,1,M);
    end
end
relative=norm(estimate(:)-expected(:))/max(1,norm(expected(:)));
assert(relative<1e-12);
assert(all(real(Ce(:))>=-1e-12));
fprintf('author Python / MATLAB relative error %.4g\n',relative);
fprintf('test_author_python_fixture PASS\n');
end
