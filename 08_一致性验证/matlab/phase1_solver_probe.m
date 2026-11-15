function phase1_solver_probe()
% Short bounded repeat with tighter inner solves; writes compact results only.
root='D:\GitHubProjects\CF-MultiTRP-WA-Robust-Optimization';
addpath(fullfile(root,'04_研究方法_WA运行实现','matlab'));
p='C:\Users\sypxx\Documents\Codex\2026-11-08\d-f-3gb-1-qq-2\outputs';
a=load(fullfile(p,'input.mat'));
cfg=a.solver;cfg.outer_iterations=8;cfg.admm_iterations=100;cfg.cg_iterations=100;
names={'matrix_pf_nominal','matrix_pf_statistical'};
for j=1:numel(names)
    cfg.branch=names{j};tic;
    sol=chasi_matrix_wmmse(a.r,cfg,{a.r.initial});
    met=chasi_evaluate(a.scene,sol,a.traffic);
    fprintf('%s time=%.3f reason=%s R5=%.9g mean=%.9g power=%.9g\n',names{j},toc,sol.reason,met.experience_bottom5_mbps_valid_ues,met.experience_mean_mbps_valid_ues,sum(sol.power));
    fprintf('trueSE=%s estimatedSE=%s\n',mat2str(sol.true_per_ue_rate',7),mat2str(sol.estimated_per_ue_rate',7));
    for t=1:numel(sol.trajectory)
        h=sol.trajectory(t);
        if isfield(h,'admm_last') && ~isempty(h.admm_last)
            fprintf('outer=%d obj=%.8g inner=%d primal=%.3g dual=%.3g ep=%.3g ed=%.3g cg=%.3g\n',h.outer,h.objective,h.admm_last.inner,h.admm_last.primal,h.admm_last.dual,h.admm_last.eps_primal,h.admm_last.eps_dual,h.admm_last.cg_relative);
        end
    end
end
end
