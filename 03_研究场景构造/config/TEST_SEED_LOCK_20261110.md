# 独立测试 seed 预先锁定

开发输入种子：场景 19、测量 1701、初始化 901。

独立测试输入种子：场景 **41**、测量 **2701**、初始化 **1901**，业务 seed 随场景为 41。测试配置为 `srs_overlap_jakes_v3_test_seed41.json`，SHA-256 `CF762E69E6D9D051759AD6CBCF5EBF14EACB40BA2306B72740FE7FAB0EB06ECA`。除三个 seed 外，其模型、SRS、资源、功率、Rank 与停止阈值和开发配置相同。此 seed 在读取任何测试性能前声明，不能据结果筛选。

运行 `matlab -batch "run('run_srs_v3_test_seed41.m')"`；它在新输入存在时校验哈希并从每支 `latest.mat` 续跑，第一轮目标为每支 20 步，**不能将 20 步上限视为收敛**。测试结果不得反过来用于调参；若后来修改设计，需另设完全独立的最终测试 seed。
