# 当前交付与验收状态（2026-11-06）

D 盘目标工程已按研究内容放入真实代码。当前 `TARGET_CODE_AUDIT.csv` 记录 964 个 .m/.py/.c/.cpp/.h 文件；其中 09 历史隔离区 736 个、08 验证区 82 个。旧 target 对照的 13 个文件已从 03/04/08 活动或当前验证路径重分类到 09，内容和来源均保留。活动入口 `startup.m` 仅加入 01–05 直接 MATLAB 目录。逐文件代码功能、哈希与分类见 `TARGET_CODE_AUDIT.csv`、`HISTORY_09_CLASSIFICATION.csv`。

## 来源

七处来源 `SOURCE_CODE_AUDIT.csv` 的 1,932 条代码逐条复核：1,900 条有同 SHA256 目标文件，32 条有逐文件未迁入理由，原文件哈希异常和无解释条目均为 0。非代码来源 14,976 条另按配置、数据、资料、日志及其他归档清单；未迁入条目保留原位。详见 `SOURCE_COVERAGE_VERIFIED.csv`、`SOURCE_NONCODE_COVERAGE.csv` 与 `UNMIGRATED_CATEGORIES.md`。原工程、F 完整快照、压缩包、F staging 和隔离历史代码均保留，本轮没有删除。

## 已运行

同一 `structured_main_v1` 冻结输入，代表组 3 步 MATLAB/Python V 最大差 8.42e-14；两语言各自的磁盘 checkpoint 1+2 步与连续运行差 0。六分支两步最大 V 差 9.70e-11；uniform/variable/fixed/mapped 四模式两步最大 V 差 8.94e-14，服务掩码逐元素相同。三层 fixed MATLAB runner 的磁盘恢复 1+1 步与连续 2 步 V、真实 SE 差 0。用冻结包原 ADMM/CG 上限各 100 的 8 步对照 V 最大差 4.05e-10、目标差 2.74e-12，停止原因相同；均通过原预设数组容差。详细配置和日志见 `08_一致性验证/ACCEPTANCE_REPORT.md`。

修正了 MATLAB CG 相对阈值、ideal_pf 分支标识、人工 policyCov 在二次项重复相加、mapped 单组来源风险取值，以及三层终态 checkpoint 再次更新问题。修正前诊断结果保留，模型参数和历史数值未调改。活动路径扫描及 MATLAB/Python 实际解析位置均未发现旧工程或 F staging 运行依赖。

## 未完成

原内层上限下，Python 代表组在第 206 次外层更新记录程序收敛；MATLAB 同组仅运行并与 Python 对照到第 8 步。第 206 步的 MATLAB 全程停止位置和最终值尚未验收，因此“多步至真实停止的跨语言对照”未通过。前 8 步 MATLAB 计时约 391 秒，线性推算到 206 步约 2.8 小时，仅为容量估计。六分支及四模式未做长轨迹、多组统计或正式扫描；MATLAB 图形后处理、SRS 独立质量档输入、正式逐用户 CDF 和秒级协议时序仍缺或未验证。服务器历史 converged 仅表示当时程序状态。

同一冻结组派生的 2 UE/2 TRP 子例在 MATLAB/Python 均于第 180 次外层更新触发同一停止条件，最终目标与真实 SE 接近；但 V 最大逐元素差 4.14e-7，超过既定数组容差。差异定位到 ADMM 第 42 次内迭代附近，紧接 rho 降低后的 CG 更新，未通过修改参数掩盖。子例不能替代完整组的第 206 步 MATLAB 对照。详见 `08_一致性验证/ACCEPTANCE_REPORT.md`。

默认 Python 3.14 具备 NumPy，缺 SciPy 和 Matplotlib；它们分别影响一个 MAT 对照脚本和可选绘图，不影响已运行的核心求解入口。旧 target 对照仅在 09 单独调用时需额外 targetRoot。没有 Git 提交或推送。



