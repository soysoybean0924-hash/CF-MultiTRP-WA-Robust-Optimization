# 静态审查记录（验证状态：NOT_RUN）

| 项 | 原场景 `chasi_generate_scene.m` | v2 默认 | 判定 |
|---|---|---|---|
| 几何与调度 | 7 站、21 TRP、组调度 | 调用同一函数，仅取几何、调度、服务掩码 | 复用；不足双 TRP 时显式补最近 TRP |
| 全 UE 长期增益 | `link_gain` 生成 `gains[K,L]` | `scene.largeScaleGain` 原样取用 | 已移除上一版中位数距离幂律替换 |
| 默认路损 | `0.025*(max(d,10)/100)^(-3.2)*sector` | 同一 `gains` | 默认无随机阴影 |
| 扇区 | `0.08+0.92*max(cos(offset),0)^4` | 包含在原 `gains` | 保留 |
| UMi/UMa 可选 | LOS、路损、阴影，并除以 `noiseWattPerRb*rbTotal` | 配置 `channelMode` 透传原生成器 | 仍需单独验证；不混合默认实验 |
| 接收相关 | `0.35^abs(i-j)` | 配置默认为 0.35 | 保留默认矩阵 |
| 快衰落 | 原路径先生成零均值相关高斯真值，随后有 CSI-RS/SRS 混合估计；远链路可能截断 | v2 同 gain 的新零均值相关高斯 `H0` 和独立 innovation | 估计过程重建；v2 不截断远链路，故不能称完全相同传播样本 |
| 噪声 | 求解器接收噪声 `I`；SRS 噪声按 home gain/SNR | 同归一化，`scene.largeScaleHomeMedian/10^(pilotSnrDb/10)` | 使用项目等效单位；未按物理 SRS RE 能量标定 |
| CSI-RS | 旧生成器在其估计器内使用 | v2 丢弃旧 H/C，不融合 CSI-RS | 明确 SRS-only |

旧生成器为了取得几何和增益仍计算旧 H/C，v2 不读取这些旧估计。
`distance_power_law_ablation` 仅作为明确选择的简化传播对照；不能与
默认 `source_large_scale_gain` 混入同一配对结果。旧结果与新结果的差
也包含估计路径、SRS 发送集合和远链路处理差异，不能只归因于鲁棒优化。

## 条件统计及 WA 适用条件

单 TRP、单发射列上把全部 SRS 用户信道堆为 `x~CN(0,Q)`，
`y=Mx+n`，`n~CN(0,σ²I)`。后验均值
`Q Mᴴ(MQMᴴ+σ²I)⁻¹y`，后验误差块矩阵
`Q-Q Mᴴ(MQMᴴ+σ²I)⁻¹M Q`。同导频用户共享 `y`，其交叉块通常非零。
对独立 innovation 的单次使用快照，均值乘 `a`，误差协方差为
`a²C0+(1-a²)Q`。这是快照条件模型，不是完整多时刻 Jakes 轨迹。

不同 TRP 的先验 `Q_l` 与接收噪声互相独立，观测矩阵按 TRP 分块，
故同 UE 跨 TRP 条件误差块为零。跨用户交叉块不能据此置零；数据
符号独立时它不出现在单 UE 接收协方差的二阶项。于是当前模型下
`R_yy,k=I+Σ_i F_ki F_kiᴴ+Σ_l p_l C_kl`，其中 `F_ki` 对同一用户
数据跨服务 TRP 相干求和。固定接收器 `A` 和权重 `Omega` 时，
预编码二次项另含 `D_l=Σ_k tr(A_k Omega_k A_kᴴ C_kl)`；两端求解器
均有这两处 `C` 路径。`E` 是条件期望 MSE，`-log₂ det(E)` 是
统计 WMMSE 代理；未证明其等于 `E[SE(H_true)]`。

非理想求解路径 `branch_channels` / MATLAB `prepareResources` 选择
`H_hat`，仅理想分支选 `H_true`。真实信道在求解器末尾作独立轨迹
诊断，不能用于参数选择。`C=0` 下鲁棒接收和二次项代数退化为名义
更新；待运行的固定状态对拍会检查实际实现。
