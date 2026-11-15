function result=evaluate_srs_joint(inputPath,checkpointPath)
% Recompute true-channel rates/power without the optimization solver.
d=load(inputPath,'data'); d=d.data;
x=load(checkpointPath,'state'); s=x.state;
assert(strcmp(s.inputSha256,chasi_sha256_file(inputPath)) && ...
    strcmp(s.configSha256,d.configSha256) && strcmp(s.modelId,d.modelId));
R=d.sceneConfig.users; ue=zeros(R,1); L=d.sceneConfig.trps;
power=zeros(L,1);
for g=1:numel(d.resources)
    r=d.resources(g); w=s.precoder{g}; H=double(r.H_true);
    K=numel(r.ueIds); Nr=size(H,3); Nt=size(H,4); D=size(w,4);
    nu=double(r.resourceWeightNu);
    power=power+nu*reshape(sum(abs(w).^2,[2,3,4]),L,1);
    for k=1:K
        allCov=eye(Nr); interf=eye(Nr);
        for i=1:K
            F=complex(zeros(Nr,D));
            for l=1:L
                F=F+reshape(H(k,l,:,:),Nr,Nt)*reshape(w(l,i,:,:),Nt,D);
            end
            allCov=allCov+F*F';
            if i~=k, interf=interf+F*F'; end
        end
        rate=max(0,(log(det(allCov))-log(det(interf)))/log(2));
        ue(double(r.ueIds(k)))=ue(double(r.ueIds(k)))+nu*real(rate);
    end
end
scheduled=unique(vertcat(d.resources.ueIds));
result=struct('branchId',s.branchId,'status',s.stopReason, ...
    'perUserTrueSE',ue,'trueSumSE',sum(ue), ...
    'truePF',sum(log1p(ue(scheduled)/d.config.pfFloor)), ...
    'trueRate5',p5(ue),'powerByTrp',power, ...
    'maxPowerViolation',max(0,max(power)-d.config.powerWattPerTrp));
end

function value=p5(x)
x=sort(x(:)); t=(numel(x)-1)*.05; lo=floor(t)+1; hi=ceil(t)+1;
value=x(lo)+(t-floor(t))*(x(hi)-x(lo));
end
