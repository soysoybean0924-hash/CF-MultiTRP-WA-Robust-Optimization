function export_solver_parity(inputPath,outDir)
% Deferred MATLAB oracle from one immutable input and its shared initial w.
assert(isfile(inputPath) && ~isfolder(outDir)); mkdir(outDir);
loaded=load(inputPath,'data'); data=loaded.data;
r0=data.resources(1); ids=1:min(3,numel(r0.ueIds));
links=find(any(r0.servingMask(:,ids),2));
r=struct('H',double(r0.H(ids,links,:,:)), ...
    'H_true',double(r0.H_true(ids,links,:,:)), ...
    'C',double(r0.C(ids,links,:,:)), ...
    'C_stat',double(r0.C_stat(ids,links,:,:)), ...
    'servingMask',r0.servingMask(links,ids), ...
    'ueIds',double(r0.ueIds(ids)), ...
    'resourceWeightNu',double(r0.resourceWeightNu), ...
    'rb_count',double(r0.rb_count), ...
    'tx_rb_count',double(r0.tx_rb_count));
w=double(data.initialPrecoder{1}(links,ids,:,:));
modelId=data.modelId; inputSha256=chasi_sha256_file(inputPath);
configSha256=data.configSha256; pfFloor=data.config.pfFloor;
budget=data.config.powerWattPerTrp; users=data.sceneConfig.users;
epsAbs=data.config.admmAbsoluteTolerance;
epsRel=data.config.admmRelativeTolerance;
admmIterations=data.config.admmIterations;
cgIterations=data.config.cgIterations;
save(fullfile(outDir,'fixture.mat'),'r','w','modelId','inputSha256', ...
    'configSha256','pfFloor','budget','users','epsAbs','epsRel', ...
    'admmIterations','cgIterations','-v7');
names={'ideal_sumse','nominal_sumse','robust_sumse', ...
    'ideal_pf','nominal_pf','robust_pf'};
native={'ideal','matrix_nominal','matrix_statistical', ...
    'matrix_pf_nominal','matrix_pf_nominal','matrix_pf_statistical'};
for b=1:numel(names)
    rb=r;
    if strcmp(names{b},'ideal_pf'), rb.H=rb.H_true; end
    cfg=struct('rank',size(w,4),'budget',budget, ...
        'outer_iterations',1,'admm_iterations',admmIterations, ...
        'cg_iterations',cgIterations,'eps_abs',epsAbs,'eps_rel',epsRel, ...
        'outer_rel_tol',data.config.outerObjectiveRelativeTolerance, ...
        'outer_precoder_tol',data.config.outerCovarianceRelativeTolerance, ...
        'pf_floor',pfFloor,'users',users,'branch',native{b}, ...
        'capture_first_step',true);
    out=chasi_matrix_wmmse(rb,cfg,{w}); probe=out.firstProbe;
    H=probe.H{1}; C=probe.C{1}; A=probe.A{1};
    Omega=probe.Omega{1}; E=probe.E{1};
    F=complex(zeros(numel(ids),numel(ids),size(H,3),size(w,4)));
    Ryy=complex(zeros(numel(ids),size(H,3),size(H,3)));
    for k=1:numel(ids)
        ry=eye(size(H,3));
        for i=1:numel(ids)
            f=complex(zeros(size(H,3),size(w,4)));
            for l=1:numel(links)
                f=f+reshape(H(k,l,:,:),size(H,3),size(H,4))* ...
                    reshape(w(l,i,:,:),size(w,3),size(w,4));
            end
            F(k,i,:,:)=reshape(f,1,1,size(H,3),size(w,4));
            ry=ry+f*f';
        end
        for l=1:numel(links)
            p=sum(abs(w(l,:,:,:)).^2,'all');
            ry=ry+p*reshape(C(k,l,:,:),size(H,3),size(H,3));
        end
        Ryy(k,:,:)=reshape(ry,1,size(H,3),size(H,3));
    end
    part=probe.parts{1}; M=part.M; B=part.B; D=part.D;
    Qw=complex(zeros(size(w)));
    for k=1:numel(ids)
        mk=reshape(M(k,:,:),size(H,3),size(H,3));
        for i=1:numel(ids)
            block=reshape(F(k,i,:,:),size(H,3),size(w,4));
            for l=1:numel(links)
                updated=reshape(Qw(l,i,:,:),size(w,3),size(w,4))+ ...
                    reshape(H(k,l,:,:),size(H,3),size(H,4))'*mk*block;
                Qw(l,i,:,:)=reshape(updated,1,1,size(w,3),size(w,4));
            end
        end
    end
    for l=1:numel(links), Qw(l,:,:,:)=Qw(l,:,:,:)+D(l)*w(l,:,:,:); end
    Qw=Qw.*reshape(r.servingMask,numel(links),numel(ids),1,1);
    userWeights=probe.userWeights{1}; perUE=probe.perUE;
    initialObjective=probe.objective; designObjective=out.optimization_objective;
    trueSE=out.true_per_ue_rate; trueSumSE=out.true_sum_rate;
    truePF=sum(log1p(trueSE(double(r.ueIds))/pfFloor));
    power=out.power; admmCandidate=probe.admmCandidate{1};
    admmHistory=probe.admmHistory; updatedW=out.precoder{1};
    file=fullfile(outDir,[names{b} '.mat']);
    save(file,'H','C','F','Ryy','A','E','Omega','M','B','D','Qw', ...
        'userWeights','perUE','initialObjective','designObjective', ...
        'trueSE','trueSumSE','truePF','power','admmCandidate', ...
        'admmHistory','updatedW','-v7');
end
fid=fopen(fullfile(outDir,'identity.json'),'wt');
projectRoot=fileparts(fileparts(mfilename('fullpath')));
matlabSolverSha256=chasi_sha256_file(fullfile(projectRoot, ...
    'matlab','chasi','chasi_matrix_wmmse.m'));
fprintf(fid,'%s',jsonencode(struct('modelId',modelId, ...
    'inputSha256',inputSha256,'configSha256',configSha256, ...
    'matlabSolverSha256',matlabSolverSha256)));
fclose(fid);
end
