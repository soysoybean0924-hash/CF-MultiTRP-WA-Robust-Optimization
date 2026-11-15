# 研究背景搭建

定义用户中心 Cell-Free/多 TRP 系统的站点、天线、资源、功率、噪声和 SRS 时序规则。输入为第三部分配置 JSON；输出为 `chasi_system_rules.m` 的系统规则结构。该文件不生成随机场景；第三部分据此构造布局与信道。当前模型含简化传播与指数老化快照，不是完整 38.901/Jakes。

当前主线实际配置函数 `matlab/chasi_config.m` 已从 03 调整至此；`legacy_cf_system` 是早期搜索配置，不在默认路径。
