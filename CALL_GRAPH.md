# 调用流程

```mermaid
flowchart LR
 A[01 问题与约束] --> B[02 系统规则]
 B --> C[03 场景/信道/SRS]
 C --> D[冻结 MAT 或带哈希二进制]
 D --> E[04 复矩阵 WA / ADMM]
 D --> F[05 来源校准]
 F --> G[05 服务/beam/power 保护]
 G --> E
 E --> H[真实信道逐用户 SE 与功率]
 H --> I[诊断绘图与导出]
```

基础 MATLAB 入口为 `run_experiment.m`；三层有限步诊断为 `chasi_run_layered_small.m`。Python 分支仅为服务器对照/后处理，不从 MATLAB 包装调用。
