# 茶思：面向非理想测量的用户中心 Cell-Free 多 TRP 协同优化研究

本目录为 D 盘目标工程。实际研究代码已按 `01` 问题、`02` 背景、`03` 场景、`04` WA、`05` 鲁棒分类；`06` 保存隔离参考实现，`07` 保存结果证据，`08` 保存验证脚本，`09` 保存独有历史实现。`TARGET_CODE_AUDIT.csv` 列出每个 .m/.py 的相对路径、主要函数和哈希，`SOURCE_CODE_AUDIT.csv` 列出七处来源代码的同内容映射或未迁入原因。根目录 `CODE_CLASSIFICATION_REPORT.md` 给出各一级目录的真实代码数和调用链。没有发现 C/C++ 源码。

2026-11-06 来源覆盖与完整运行验收记录见 `08_一致性验证/ACCEPTANCE_REPORT.md`。代表组、六分支和四种鲁棒模式的限次 MATLAB/Python 数值对照与磁盘恢复已运行；代表组到程序真实停止点的 MATLAB 全程对照仍未完成，不能称整个运行验收通过。

MATLAB R2026b：在根目录运行 `startup; run_experiment('input')` 生成独立冻结输入，或 `startup; run_experiment('basic')` 做一步基础求解（可能耗时）。有限步三层诊断可调用 `chasi_run_layered_small(fullfile(pwd,'03_研究场景构造','inputs','structured_main_v1'),fullfile(pwd,'05_研究创新点_鲁棒施加','config','study.json'),'fixed',1,1)`。正式全量实验本轮不运行。所有运行结果应写新目录；原始快照不改。

可直接执行的限次例子：在 MATLAB 中切换到本目录，运行 `startup; run_experiment('layered_diagnostic')`。该入口读取明确标识的 `structured_main_v1` 冻结输入，只用一组、一步，不启动正式实验。六分支的一组限次入口在 `08_一致性验证/matlab/validate_six_branch_entry_small.m`。首步 MATLAB/Python 数值差异已追到 CG 阈值并修正，见 `08_一致性验证/REPORT.md`。全组长迭代、完整恢复和正式扫描尚未验收；历史服务器“converged”不是理论收敛或整个研究完成的证明。

诊断图：在具备 Matplotlib 的 Python 环境执行 `05_研究创新点_鲁棒施加/run_server_diagnostic_plots.ps1`，脚本设置本地 `SRS_CJK_FONT`，从本工程 07 中已归档的 state.json 读取，并输出新目录。另有 `python/plot_power_scan.py` 用于功率诊断图。MATLAB 批处理绘图在本机 connector 崩溃，故根入口不自动绘图；可在图形环境单独调用 `chasi_plot_layered_diagnostic.m`。
