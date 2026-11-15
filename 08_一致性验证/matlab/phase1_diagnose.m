function phase1_diagnose()
% Read-only diagnosis of the frozen phase-1 MATLAB run.
p='C:\Users\sypxx\Documents\Codex\2026-11-08\d-f-3gb-1-qq-2\outputs';
a=load(fullfile(p,'input.mat'));
names={'ideal_pf','matrix_pf_nominal','matrix_pf_statistical'};
for j=1:numel(names)
    d=load(fullfile(p,[names{j} '.mat']));s=d.sol;m=d.met;
    fprintf('\n%s reason=%s objective=%.9g R5=%.9g power=%.9g\n',names{j},s.reason,s.optimization_objective,m.experience_bottom5_mbps_valid_ues,sum(s.power));
    fprintf('UE expMbps trueSE estSE PFweight powerW groups service\n');
    for u=1:4
        g=ceil(u/2);k=mod(u-1,2)+1;
        rr=a.r(g);v=s.precoder{g};nr=size(rr.H_true,3);nt=size(rr.H_true,4);rk=size(v,4);
        signal=0;interference=0;
        for i=1:size(v,2)
            eff=complex(zeros(nr,rk));
            for l=1:size(v,1)
                eff=eff+reshape(rr.H_true(k,l,:,:),nr,nt)*reshape(v(l,i,:,:),nt,rk);
            end
            if i==k,signal=norm(eff,'fro')^2;else,interference=interference+norm(eff,'fro')^2;end
        end
        fprintf('%d %.9g %.9g %.9g %.9g %.9g %d %d signal=%.6g interference=%.6g noise=%.6g\n',u,m.per_ue_experience_mbps(u),s.true_per_ue_rate(u),s.estimated_per_ue_rate(u),s.userWeight{g}(k),m.per_ue_weighted_tx_power_watt(u),m.per_ue_scheduled_groups(u),m.per_ue_service_trp_occurrences(u),signal,interference,nr);
    end
    for t=1:numel(s.trajectory)
        h=s.trajectory(t);
        if isfield(h,'admm_last')
            fprintf('outer=%d obj=%.9g change=%.3g admm=%d\n',h.outer,h.objective,h.precoder_change,h.admm_converged);
            disp(h.admm_last)
        end
    end
end
for g=1:2
    r=a.r(g);fprintf('group=%d nu=%.9g UE=%s noise=identity\n',g,r.resourceWeightNu,mat2str(r.ueIds'));
    for k=1:2
        mask=find(r.servingMask(:,k));
        h=reshape(r.H(k,:,:,:),size(r.H,2),[]);
        e=reshape(r.H_true(k,:,:,:)-r.H(k,:,:,:),size(r.H,2),[]);
        c=reshape(r.C(k,:,:,:),size(r.C,2),[]);
        fprintf('UE=%d home=%d candidates=%s Hnorm=%.6g Enorm=%.6g Cnorm=%.6g\n',r.ueIds(k),r.home_trp(k),mat2str(mask'),norm(h,'fro'),norm(e,'fro'),norm(c,'fro'));
    end
end
end
