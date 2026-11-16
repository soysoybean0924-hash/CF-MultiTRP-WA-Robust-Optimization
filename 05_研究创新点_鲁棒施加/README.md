# 研究创新点：鲁棒施加

现行 v6 MATLAB 主线中，sum-SE 与 PF 的目标、结果口径、老化条件预测、SRS/CSI-RS 分源误差及与已有工作的逐项对比，见 [SRS 分源误差与鲁棒设计说明](docs/SRS_ERROR_SOURCE_ROBUST_COMPARISON.md)。该说明明确区分已接入的 PF 正式入口与尚未接入的 sum-SE 正式实验。

“创新点”仅是研究分类，不代表新颖性或收益已证实。MATLAB `chasi_layer1_calibrate.m` 对资源组开发/测试划分做经验均方球校准；`chasi_layer_protection.m` 按 Python `layered.py` 的 fixed/mapped/uniform/variable 规则计算强度，将额外 `policyCov` 加入设计协方差并进入接收滤波、MSE 权重和预编码二次项、按 `chasi_layer_service.m` 做每 UE 最多一次服务替换，并收紧有效功率上限。物理 C 保存在 `C_stat`，不能把来源份额再加成物理 C。`chasi_run_layered_small.m` 是有限步诊断入口，尚非 Python 全套可恢复正式 campaign。Python 分支保存服务器原型及 CDF 后处理。强度、消融、触发、功率扫描仍须在完成等价验证后接入 MATLAB 正式入口。


2026-11-15 存档说明：本文提到的 Python 文件已从 main 移至 archive/python-before-matlab-only-main；main 的 v6 正式链使用 MATLAB。
