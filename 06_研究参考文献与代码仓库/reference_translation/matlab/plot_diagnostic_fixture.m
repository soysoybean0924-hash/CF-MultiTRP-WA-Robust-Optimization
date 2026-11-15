function plot_diagnostic_fixture()
% Requires a real run of run_diagnostic_fixture; no synthetic plotting.
root=fileparts(fileparts(mfilename('fullpath')));
result=fullfile(root,'reproduction','matlab','diagnostic_result.mat');
assert(isfile(result),'Run run_diagnostic_fixture in MATLAB first.');
s=load(result); out=s.out;
input=load(fullfile(root,'reproduction','inputs','diagnostic_fixture.mat'));
f=figure('Color','w');
plot(0:numel(out.objective)-1,out.objective,'o-','LineWidth',1.5);hold on;
plot(0:numel(input.python_objective_trace)-1, ...
    input.python_objective_trace,'x--','LineWidth',1.2);
grid on; xlabel('WMMSE outer iteration');
ylabel('Sum spectral efficiency (bit/s/Hz; prelog omitted)');
title('Synthetic diagnostic fixture, not paper reproduction');
legend('MATLAB','Original Python','Location','southeast');
savefig(f,fullfile(root,'reproduction','matlab','diagnostic_objective.fig'));
exportgraphics(f,fullfile(root,'reproduction','matlab','diagnostic_objective.png'));
f=figure('Color','w');
for t=1:size(out.primal,1)
    n=out.innerCount(t);
    semilogy(1:n,out.primal(t,1:n),'LineWidth',1.2); hold on;
end
grid on; xlabel('ADMM inner iteration'); ylabel('Primal residual ||mu-q||_F');
title('Per-outer-loop ADMM residuals (separate subproblems)');
savefig(f,fullfile(root,'reproduction','matlab','diagnostic_primal.fig'));
exportgraphics(f,fullfile(root,'reproduction','matlab','diagnostic_primal.png'));
end
