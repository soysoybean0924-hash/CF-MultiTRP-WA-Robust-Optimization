# 研究场景构造

MATLAB 主入口 `prepare_srs_joint_input.m`：场景、信道、SRS 复用、LMMSE/老化、冻结 MAT 与 SHA256。`chasi_load_frozen_binary.m` 按原 MATLAB 列优先 little-endian float64 格式读服务器冻结输入并逐文件校验。输入身份严格区分 `inputs/structured_main_v1` 与仅登记 manifest 的 `inputs/srs_v2`；后者不具备运行数据，不能混用。Python `srs_joint.py`、`frozen_input.py` 是服务器观测与二进制对照。输出供第四、五部分使用。
