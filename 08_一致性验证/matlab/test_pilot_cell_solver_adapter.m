function test_pilot_cell_solver_adapter()
% One RB group, identical solver input contract for ideal and estimated CSI.
K=2; L=2; Nr=4; Nt=2;
Htrue=cell(K,L,1); R=complex(zeros(Nr,Nr,K,L));
for l=1:L
    for k=1:K
        h=complex(zeros(Nr,Nt));
        h(1,:)=[1+0.1*k,0.2*l];
        h(2,:)=[0.2*k,1+0.1*l];
        Htrue{k,l,1}=h;
        R(:,:,k,l)=eye(Nr);
    end
end
pilots=[1 2;1 2]; energy=ones(K,2); noise=1;
Np=complex(zeros(Nt,L,2));
[Hhat,Ce,crossCe]=chasi_pilot_estimate_cell_group( ...
    Htrue,R,pilots,energy,noise,[1 2],Np);
assert(isequal(size(Hhat),size(Htrue)));
assert(isequal(size(Hhat{1,1,1}),[Nr,Nt]));
assert(norm(Hhat{1,1,1}-Htrue{1,1,1},'fro')>1e-3);
assert(isequal(size(Ce),[Nr,Nr,K,L]));
assert(isequal(size(crossCe),[Nr,Nr,K,K,L]));
[knownMean,~]=chasi_pilot_estimate_cell_group( ...
    Htrue,R,pilots,energy,noise,[1 2],Np,Htrue);
for l=1:L
    for k=1:K
        assert(norm(knownMean{k,l,1}-Htrue{k,l,1},'fro')<1e-12);
    end
end
service=[1;2]; Pmax=1; streams=1;
[X0,~]=joint_matrix_common_initial(service,Nt,streams,1,Pmax,901);
opts={'maxOuter',2,'maxInner',40,'rho',0.1,'initialX',X0, ...
    'objectiveKind','sumse','stopPolicy','paired'};
ideal=joint_matrix_wmmse_admm(Htrue,service,1,Pmax,streams,opts{:});
nominal=joint_matrix_wmmse_admm(Hhat,service,1,Pmax,streams,opts{:});
assert(isequal(size(ideal.X),size(nominal.X)));
assert(all(isfinite(ideal.userSE)) && all(isfinite(nominal.userSE)));
[trueIdeal,~]=joint_matrix_rate(Htrue,ideal.X,service,1);
[trueNominal,~]=joint_matrix_rate(Htrue,nominal.X,service,1);
assert(isfinite(trueIdeal) && isfinite(trueNominal));
fprintf('test_pilot_cell_solver_adapter PASS: Htrue SE %.6g / Hhat-design SE %.6g\n', ...
    trueIdeal,trueNominal);
end
