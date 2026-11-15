# Multi-TRP 分层实验独立包（本地代码未运行；服务器待验证）

本包基于 `srs_structured_20261103` 的 Python 求解器副本。只复用服务器原有冻结二进制，不写原结果目录。复制的 `native_solver.py` 保持原精度判定；新 `diagnose_nominal.py` 在独立目录记录 ADMM 原始/对偶残差、阈值、逐资源组 CG 精度，先定位旧 nominal 第一步不精确原因，不降低容差。`layered.py`、`preflight.py`、`plot_layered.py` 均**待服务器运行验证**。

## 物理和数学边界

现有冻结输入只有一套 scene/measurement/initial seed (19/1701/901)，每组一份真实信道。`source_cov` 四项为 residual、pilot、noise、age；总物理后验 C 单独保留。第一层对不同资源组作开发/测试划分，以单次误差能量的经验分位形成逐链路均方球半径；来源份额只作解释与映射，不再次相加到 C。该覆盖是**单场景组留出**，非独立信道/SRS 种子的物理覆盖保证；跨来源相关项没有作为独立集合使用。

第二层使用相同 Layer 1 输出。两个组均通过矩阵波束附加 PSD 设计项 P、有限服务替换评分/门槛，以及收紧的可行功率球起作用。预编码 `V` 的方向与幅度仍由原矩阵 WMMSE–ADMM 更新，功率是 `V` 的加权平方范数；所谓功率保护是附加保留余量，物理每 TRP 20 W 上限仍在。P 是**启发式设计项**，服务门槛和功率余量是本轮设计选择；这三个操作不是由均方球严格推导的鲁棒对应约束。`uniform` 共用一个无量纲等级；`variable` 使用配置中的三项因子。若出现收益，只能解释为整体原型对照。

第三层固定同一 `variable_factors` 和三变量操作，`fixed` 用开发配置的固定等级，`mapped` 根据四项来源份额的加权平均、上下界映射等级；`mapping` 公式及权重为本轮设计选择，原《误差建模and鲁棒设计.md》只提出与 SRS SINR、周期、复用密度相关的方向，未给可直接复现的完整系数。现有输入未导出每 UE 的真实 SRS SINR，因此横轴可用配置 `pilotSnrDb`，不能写成实测 SRS SINR。当前单档不能形成质量扫描。

所有诊断图在 `STATUS.json` 标记 `formal_result=false`。`SE<0.01` 是速率代理阈值，不能充当秒级掐头去尾体验约束。零速率包含未调度用户。正式主图所需独立 seed、严格体验约束、可比理想/名义/非协同基线和新 SRS 档输入见 [逐图数据依赖.md](逐图数据依赖.md)。`srs_reference_papers.zip` 在指定本地范围与 Downloads 未找到，本轮未据该包声称论文复现。

## ① Windows PowerShell：上传 ZIP（生成 ZIP 见交付目录）

```powershell
$zip = 'C:\Users\sypxx\Documents\CF mMIMO TRP\layered_multitrp_20261104.zip'
scp -P 22 $zip user@10.112.48.104:/data/xiwai/chasi_pilot_mmse_upload/
```

## ② Linux：解包、独立环境、预检与第一层

```bash
cd /data/xiwai/chasi_pilot_mmse_upload
mkdir -p layered_multitrp_20261104
cd layered_multitrp_20261104
unzip -n ../layered_multitrp_20261104.zip
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-server.txt
export SRS_CJK_FONT=/data/xiwai/chasi_pilot_mmse_upload/postprocess_20261104/NotoSansSC.otf
ROOT=/data/xiwai/chasi_pilot_mmse_upload/srs_structured_run_20261103_k9h1jp/srs_structured_linux_20261103
INPUT="$ROOT/output_srs_joint_main_v1/frozen_binary"
CFG=config/study.json
.venv/bin/python python/preflight.py --input "$INPUT" --config "$CFG" --out results/preflight
.venv/bin/python python/diagnose_nominal.py --input "$INPUT" --out results/nominal_diagnostic
.venv/bin/python python/layered.py calibrate --input "$INPUT" --config "$CFG" --out results/layer1
```

先审阅 `results/preflight/preflight.json`、`results/nominal_diagnostic/nominal_admm.json`、`results/layer1/layer1.json`。旧 nominal 的诊断如显示 CG/ADMM 不精确，需要再修复并重新审查精度；本包未静默放宽容差。

## ③ Linux：小规模机制试跑（各命令独立）

```bash
cd /data/xiwai/chasi_pilot_mmse_upload/layered_multitrp_20261104
ROOT=/data/xiwai/chasi_pilot_mmse_upload/srs_structured_run_20261103_k9h1jp/srs_structured_linux_20261103
INPUT="$ROOT/output_srs_joint_main_v1/frozen_binary"
LAYER=results/layer1/layer1.json
CFG=config/study.json
.venv/bin/python python/layered.py run --input "$INPUT" --config "$CFG" --layer "$LAYER" --mode uniform --max-groups 2 --max-outer 2 --out results/smoke/layer2_uniform
.venv/bin/python python/layered.py run --input "$INPUT" --config "$CFG" --layer "$LAYER" --mode variable --max-groups 2 --max-outer 2 --out results/smoke/layer2_variable
```

试跑状态若为 `inner_inexact`、`failed` 或 `limit_unconverged`，按状态保存，不能称收敛增益。当前 runner 是单进程集中式顺序 ADMM，外层每步保存最终 `state.json`/`checkpoint.npz`；长运行中断后本包**尚不支持续算**，须检查预算和输出目录后单项重跑。

## ④ Linux：第二层、第三层主对照（先于补充）

```bash
cd /data/xiwai/chasi_pilot_mmse_upload/layered_multitrp_20261104
ROOT=/data/xiwai/chasi_pilot_mmse_upload/srs_structured_run_20261103_k9h1jp/srs_structured_linux_20261103
INPUT="$ROOT/output_srs_joint_main_v1/frozen_binary"
LAYER=results/layer1/layer1.json
CFG=config/study.json
.venv/bin/python python/layered.py run --input "$INPUT" --config "$CFG" --layer "$LAYER" --mode uniform --out results/main/layer2_uniform
.venv/bin/python python/layered.py run --input "$INPUT" --config "$CFG" --layer "$LAYER" --mode variable --out results/main/layer2_variable
.venv/bin/python python/plot_layered.py --kind layer2 --states results/main/layer2_uniform/state.json results/main/layer2_variable/state.json --out results/figures/layer2
.venv/bin/python python/layered.py run --input "$INPUT" --config "$CFG" --layer "$LAYER" --mode fixed --out results/main/layer3_fixed
.venv/bin/python python/layered.py run --input "$INPUT" --config "$CFG" --layer "$LAYER" --mode mapped --out results/main/layer3_mapped
.venv/bin/python python/plot_layered.py --kind layer3 --states results/main/layer3_fixed/state.json results/main/layer3_mapped/state.json --out results/figures/layer3
```

首先生成的两组图文件名为 `layer2_diagnostic_*.png`、`layer3_diagnostic_*.png` 和 `layer3_strength_*.png`，**属于机制原型诊断图**。满足独立测试与严格约束前不改名为正式论文主图。

## ⑤ Linux：三项补充实验（各自单独运行）

```bash
cd /data/xiwai/chasi_pilot_mmse_upload/layered_multitrp_20261104
ROOT=/data/xiwai/chasi_pilot_mmse_upload/srs_structured_run_20261103_k9h1jp/srs_structured_linux_20261103
INPUT="$ROOT/output_srs_joint_main_v1/frozen_binary"
LAYER=results/layer1/layer1.json
CFG=config/study.json
for v in 0 0.05 0.1 0.15 0.2; do tag=$(printf '%s' "$v" | tr '.' 'p'); .venv/bin/python python/layered.py run --input "$INPUT" --config "$CFG" --layer "$LAYER" --mode mapped --level "$v" --out "results/strength/$tag"; done
for v in beam service power; do .venv/bin/python python/layered.py run --input "$INPUT" --config "$CFG" --layer "$LAYER" --mode mapped --ablation "$v" --out "results/ablation/$v"; done
for v in off fixed proposed; do .venv/bin/python python/layered.py run --input "$INPUT" --config "$CFG" --layer "$LAYER" --mode mapped --trigger "$v" --out "results/trigger/$v"; done
```

强度扫描的 `--level` 是无量纲基准等级；各变量实际等级、服务替换次数与有效功率上限写在每个 `state.json.mechanism`。三个触发组复用同一原始候选规则和每 UE 最多一次替换预算，分别关闭、使用固定门槛、使用拟议门槛。消融只关闭额外保护，原始服务关系与物理功率可行性仍保留。

## ⑥ Linux：功率上限扫描；SRS 质量扫描输入缺口

```bash
cd /data/xiwai/chasi_pilot_mmse_upload/layered_multitrp_20261104
ROOT=/data/xiwai/chasi_pilot_mmse_upload/srs_structured_run_20261103_k9h1jp/srs_structured_linux_20261103
INPUT="$ROOT/output_srs_joint_main_v1/frozen_binary"
for cap in 10 15 20; do .venv/bin/python python/layered.py run --input "$INPUT" --config config/study.json --layer results/layer1/layer1.json --mode mapped --power-cap "$cap" --out "results/power_cap/${cap}W"; done
```

此入口扫**功率上限**，不是功率最小化。SRS 质量扫描不能从同一冻结输入重命名产生。需在服务器具备 MATLAB、原 MATLAB 生成链、原场景依赖与独立输出目录后，按 `config/study.json.quality_scan_db` 为每档重新生成观测、LMMSE、C 和来源统计；当前 ZIP 没有可在纯 Python 端重生成 17 组输入的完整 MATLAB 依赖，因此未提供会伪造扫参的命令。

## ⑦ Windows PowerShell：下载

```powershell
$dest = 'C:\Users\sypxx\Documents\CF mMIMO TRP\layered_multitrp_20261104_server_results'
New-Item -ItemType Directory -Force -Path $dest | Out-Null
scp -P 22 -r user@10.112.48.104:/data/xiwai/chasi_pilot_mmse_upload/layered_multitrp_20261104/results/. $dest
```

回传首轮 `preflight.json`、`nominal_admm.json`、`layer1.json` 的摘要与两个 smoke `state.json` 后，再依据真实错误修复。所有原始检查点、manifest 和 SHA256SUMS 保持只读。

2026-11-15 存档说明：本文提到的 Python 文件已从 main 移至 archive/python-before-matlab-only-main；main 的 v6 正式链使用 MATLAB。
