# 研究方法 WA 运行实现

`chasi_matrix_wmmse.m` 为 MATLAB 复矩阵 WMMSE 与集中式 ADMM 主实现；接收滤波、MSE 权重、预编码二次型、TRP 功率投影、残差与停止条件在其中。`run_srs_joint.m`、`run_srs_joint_batch.m` 为服务器原入口，原路径与代码哈希契约需迁移适配，不可直接作为本工程统一入口。根目录 `run_experiment.m` 提供新输入和一步基础求解。Python `native_solver.py`、`campaign.py` 用于服务器数值对照。参考仓库是标量功率 WA 起点；当前复矩阵实现是项目扩展，不能称原仓库直接实现。

`run_srs_joint.m` 与 `run_srs_joint_batch.m` 已按本工程目录改写路径/代码哈希，`status` 从 D 盘运行通过。完整六分支求解尚未正式运行；最小六分支检查见 08。
