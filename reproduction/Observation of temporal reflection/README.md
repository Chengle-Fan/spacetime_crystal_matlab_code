# 光子时间界面的时间反射与宽带频移：MATLAB 复现

本目录是对 H. Moussa *et al.* 论文 **“Observation of temporal reflection and broadband frequency translation at photonic time interfaces”** 的独立 MATLAB 复现：

- Hady Moussa, Gengyu Xu, Shixiong Yin, Emanuele Galiffi, Younes Ra’di and Andrea Alù, *Nature Physics* **19**, 863–868 (2023)
- DOI：[10.1038/s41567-023-01975-y](https://doi.org/10.1038/s41567-023-01975-y)
- Nature 文章页：[s41567-023-01975-y](https://www.nature.com/articles/s41567-023-01975-y)
- 本地主文：[Observation of temporal reflection and.pdf](<./Observation of temporal reflection and.pdf>)
- 本地补充材料：[supplementary_information.pdf](./supplementary_information.pdf)
- Nature 官方 Source Data 本地归档：[source_data.zip](./source_data.zip)
- 已解压的官方数据：[data/Time_Interface_Source_Data](./data/Time_Interface_Source_Data/)

补充材料 PDF 和 `source_data.zip` 均来自论文 DOI/Nature 在线内容。本目录没有把作者未公开的 ADS 工程、示波器采集设置或实验数据“反推”为伪原始数据。

## 复现范围与结果标签

这里的“复现”分为五类，图题和下表会明确区分：

1. **官方数据复画**：直接读取 Nature Source Data，允许的处理只包括文中说明或源码中显式保存的时间门、基线、窗函数、坐标换算和显示对齐。
2. **S1–S17 解析模型**：实现补充材料的单元 ABCD、Bloch 色散/阻抗、守恒波数频率映射、两种微观时间边界以及时间板四路径干涉。
3. **离散 Kirchhoff/MNA 电路**：对真实节点电压和支路电流列 KCL/支路方程，以梯形法求解开关传输线；它不是连续介质 FDTD。
4. **合成波包/方程预测**：用于展示时序反转、负极性、固定波数频移、有限上升沿、时间 Fabry–Pérot 等机制；没有公开原始数据的面板不会标成实验复画。
5. **不可唯一数值重建**：实验照片、PCB/Gerber、完整 BOM、缺失的端接器 S 参数、未公开的门函数/损耗补偿及五次重复测量等。

所有代码均为 MATLAB R2019b+ 基础语法，不要求额外工具箱。当前工作区已在 MATLAB R2026a 下运行自动验收；参见 [`output/validation_summary.txt`](./output/validation_summary.txt)。

## 快速开始

在 MATLAB 中进入本目录后，一键运行全部官方数据图、解析图、MNA 演示和验收：

```matlab
cd('/Users/chandler/spacetime_crystal_matlab_code/reproduction/Observation of temporal reflection')
summary = run_all_reproductions('all');
```

不读取三组大型示波器 CSV、也不运行 30 单元 MNA 场图的快速模式：

```matlab
summary = run_all_reproductions('quick');
```

在当前机器的无界面运行中，`quick` 通常约 20 秒；优化后的三个官方数据脚本合计约 9 秒。本次完整入口（含 30 单元 MNA）实测约 29 秒。总耗时仍取决于 MATLAB 启动、CSV 读取、绘图后端和机器性能，这些时间只作量级参考而非性能验收门槛。

只做数值验收：

```matlab
report = otr_validate();
```

当前机器也可从终端无界面运行：

```bash
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('/Users/chandler/spacetime_crystal_matlab_code/reproduction/Observation of temporal reflection'); run_all_reproductions('all');"
```

`run_all_reproductions` 默认在首个错误处停止，并写入运行清单。需要批量诊断时可使用：

```matlab
summary = run_all_reproductions( ...
    'IncludeSourceData', true, ...
    'IncludeCircuit', true, ...
    'RunValidation', true, ...
    'CloseGeneratedFigures', true, ...
    'ContinueOnError', true);
```

## 论文参数与歧义处理

默认值集中在 [`otr_parameters.m`](./otr_parameters.m)；MNA 的电路参数也可通过 `otr_build_mna(options)` 覆盖。

| 参数 | 本复现采用值 | 论文出处/说明 |
|---|---:|---|
| 单元数 | 30 | 主文和 SI |
| 单元物理长度 `d` | **0.2080 m** | Methods、Extended Data Fig. 1 与 SI Fig. S1 均给出 0.2080 m；主文写约 **20.1 cm**，两者冲突。本复现采用重复出现且位数更明确的 0.2080 m |
| 总长度 | 6.24 m | `30 × 0.2080 m`，用于 Fig. 1e 相位到波数的换算 |
| 无载线阻抗 `Z0` | 50 Ω | Methods；设计点 `Zoff≈50 Ω`、`Zon≈25 Ω` |
| 有效介电常数 `epsEff` | 8.36 | Methods 的微带有效介电常数；不要与基板标称介电常数混淆 |
| 基板 | 1.52 mm Rogers TMM 13i | Methods；标称相对介电常数 `12.85±0.35` |
| 设计电长度 | 72° @ 100 MHz | SI S1 / Extended Data Fig. 1 |
| 归一化并联电纳 `b` | 2.753 @ 100 MHz | SI S1；理想 Bloch 模型使用这一独立给定值 |
| 开关并联负载 | 82 pF | Methods |
| ON 电阻 | 1 Ω | Methods 的理想电压控制反射式开关 |
| 负载寄生串联电阻 | 4 Ω | Methods |
| 接地过孔寄生电感 | 8 nH | Methods |
| 最近邻寄生耦合 | 2 pF | Methods，可在离散 MNA 中开关启用 |
| 次近邻寄生耦合 | 1 pF | Methods，可在离散 MNA 中开关启用 |
| 介质损耗角正切 | 0.0019 | Methods |
| 传播衰减 `A` | **0.5 dB/m @ 100 MHz** | Methods 的 ADS 参数；`fig08` 以 `0.5·ln(10)/20` 转为导体损耗 Np/m，以 `sqrt(f/100 MHz)` 延拓，并另加 `beta·tan(delta)/2` 介质损耗 |
| 实验界面时间 | 约 3 ns（10%–90%） | 主文和 SI S9；3 ns 点的 TR 峰值约为理想突变的 90% |
| 控制线相对传播时间 | 主信号线的约 1/80 | 主文；SI S8 表述为短 40 倍以上且从中心馈电再减半 |
| 时间板持续时间 | 15、25、35 ns | 主 Fig. 3c/d；定义为函数发生器控制信号的 ON 时间，不等于慢控制电压边沿间的几何宽度 |

### 两处必须保留的文献冲突

**单元长度。** 主文正文写“约 20.1 cm”，Methods、SI 和 Extended Data 均写 `0.2080 m`。相位反演和电路模型统一使用 `0.2080 m`，因此样品长度为 `6.24 m`；README 不把两者伪装成舍入误差。

**上升沿标注。** Extended Data/SI Fig. S7 的图内曲线标签可读为约 `0、2、3、8、16 ns`，而同一图的文字说明写“进一步增大到 `8 ns` 和 `12 ns`”。因此存在 **8/12/16 ns** 的内部不一致。`fig07_finite_switching_and_homogeneity` 扫描连续的 `0–16 ns`，只把文中无歧义的 `3 ns、约 90%` 作为实验锚点，不声称生成了作者的 12 ns 或 16 ns ADS 原始曲线。MNA 调度器允许用户显式选择任意有限上升沿。

另有两套速度参数需要区分：S1–S3 的设计模型按 `theta(100 MHz)=72°` 建立；MNA 将 Methods 的 `epsEff=8.36` 转为单元 `L/C`。两者回答不同层级的问题，不应强行拟合为同一条曲线。

## 官方 Source Data 覆盖情况

`source_data.zip` 的实际内容如下。空目录也保留，因为它们直接限定了哪些论文面板不能由公开数据重建。

| 论文面板 | ZIP 中的内容 | 本复现处理 |
|---|---|---|
| Fig. 1e | `1e_raw/`：ON/OFF 实测与 ADS 相位文本 | `fig01_source_data` 以 `beta=-phase/L` 复画色散 |
| Fig. 1f | `1f_raw/tek0000.csv` 与 `V_in_sim.txt`、`V_out_sim.txt` | `fig01_source_data` 复画示波器/ADS 波形；仅移动实测 TIME 轴以对齐入射峰 |
| Fig. 2a | `2a/` **空目录** | 只能由 S1–S3、S16 理论预测；无实测系数原始数据 |
| Fig. 2b | `2b/` **空目录** | 只能由 S1–S3、S17 理论预测；无实测系数原始数据 |
| Fig. 2c | `2c_raw/tek0002.csv` | `fig02_source_data` 显式门控、FFT，复现红移 |
| Fig. 2d | `2d_raw/tek0000.csv` | `fig02_source_data` 显式门控、FFT，复现蓝移 |
| Fig. 2e | `2e/freq_conv.csv` | `fig02_source_data` 复画 30–70 MHz 扫频 |
| Fig. 2f | `2f/freq_conv_reverse.csv` | `fig02_source_data` 复画 20–34 MHz 反向扫频 |
| Fig. 3c/d | `3cd_raw/tau_15.csv`、`tau_25.csv`、`tau_35.csv` | `fig03_source_data` 只复画原始 CH1/CH2/开关监测波形；不足以唯一得到已处理的 Fig. 3d |
| Fig. 3e | `3e/` **空目录** | `fig11_slab_wavepackets` 只给出解析/合成预测，无五次实验点和范围 |
| Extended Data / SI | ZIP 中无对应目录 | 仅做解析或定性预测，并在图中注明无原始数据 |

## 模型分层

### 1. S1–S17 解析层

理想加载单元的归一化 ABCD 矩阵由 `otr_abcd_unit_cell` 严格按 S1 组装；`otr_bloch_dispersion` 使用

```text
cos(beta*d) = cos(theta) - b*sin(theta)/2                 (S2)
ZB = ±B*Z0/sqrt(A^2-1)                                   (S3)
```

并显式选择被动第一通带的 Bloch 分支。`realistic` 选项把理想电纳替换为 `(1+4) Ω – 8 nH – 82 pF` 串联支路的频率响应。

S13–S17 的核心关系为：

```text
omega2/omega1 = sqrt(C1/C2) = Z2/Z1                      (S13)
C2*v(0+) = C1*v(0-)        （接入电容，电荷连续）         (S14)
v(0+) = v(0-)              （移除电容，电压连续）         (S15)
```

在 `Zoff=50 Ω、Zon=25 Ω` 的理想设计点，S16/S17 给出：

```text
Ton = 0.375,  Ron = -0.125
Toff = 1.5,   Roff = -0.5
```

两个负反射系数解释了时间反射的负极性。`otr_temporal_slab` 再显式保留 `TT、RR、RT、TR` 四条因果路径；时间板后的总前向波是 `TT+RR`，总反向波是 `RT+TR`。

S1–S17 的可执行覆盖边界如下：

| 方程 | 状态 | 对应代码 |
|---|---|---|
| S1 | 精确数值实现 | `otr_abcd_unit_cell` |
| S2–S3 | 精确数值实现与分支检查 | `otr_bloch_dispersion`、`fig09_tlm_design` |
| S4–S5 | 文献公式可读，但不伪造完整校正 | 论文仅给连接器在 10–100 MHz 约 `S≈0.665`；缺少实测频率相关 S 参数，官方波形脚本只做透明显示/对齐 |
| S6–S9 | 门控、FFT 和固定实波数映射可执行；完整系数反演不可唯一执行 | `otr_gate_signal`、`otr_fft_spectrum`、`fig02_source_data`、`otr_frequency_map`；缺少逐实验的传播距离、损耗和相位补偿输入 |
| S10–S12 | 通过其时域初值结论数值验证，不另做符号拉普拉斯演算 | `fig04_boundary_conditions`、MNA 状态跳变 |
| S13 | 守恒波数频率映射 | `otr_frequency_map`、`fig04_boundary_conditions` |
| S14–S15 | 精确事件边界 | `otr_temporal_interface`、`otr_simulate_circuit`、`otr_validate` |
| S16–S17 | 精确散射矩阵和符号 | `otr_temporal_interface`、`otr_temporal_slab`、`fig10_scattering_spectra` |

### 2. 离散 Bloch 层

`otr_frequency_map` 不是简单地对所有频率乘固定阻抗比，而是在目标状态第一通带求解

```text
Re[beta_after(fout)] = Re[beta_before(fin)]
```

这能体现有限单元长度和串联 RLC 的色散，并检验 OFF→ON 红移、ON→OFF 蓝移以及往返映射。`fig09_tlm_design` 还复现设计参数图、首个带隙约束、Bloch 阻抗和固定波数频率比。

### 3. Kirchhoff/MNA 电路层

`otr_build_mna` 构造真正的离散电路，而非把连续介质方程换一个名字：

- 30 个基线串联 `L/R` 单元、31 个节点 `C/G`；
- 输入端为 Thevenin 电压源与 50 Ω 内阻的等效 Norton stamp，输出为 50 Ω 负载；
- 最近邻 2 pF 和次近邻 1 pF 是节点间互电容；
- 状态量为节点电压、线路电感电流，以及 RLC 模式下的开关支路电流/电容电压；
- 每一行均对应 KCL 或器件支路方程。

求解的描述系统是

```text
C(w) * xdot + G(w) * x = B * vSource(t)
```

其中 `w` 是每个开关的连接权重，使用梯形法分段积分。`otr_switch_schedule` 强制理想事件落在时间网格上，防止事件被悄悄延后一个步长。

两种物理模式：

- `physicsMode='ideal-capacitance'`：默认模式。ON 事件执行节点电荷共享 `Coff*v- = Con*v+`，线路电感电流连续；OFF 事件保持电压连续并从线路系统移除支路电荷，严格实现 S14/S15。有限上升沿被离散成可审计的增量连接事件。
- `physicsMode='rlc'`：使用固定最大拓扑，ON 时连接串联 `(Ron+Rpar)-Lpar-Cload` 支路；OFF 时隔离并清零被移除支路的内部状态。它适合检查固定 ON 段、寄生和有限连接权重，但不能替代作者未公开的 ADS 开关宏模型。

`otr_circuit_energy` 独立计算电容/电感储能、源功、50 Ω/线路/支路损耗与开关事件功，验算

```text
E(t)-E(0) = integral(Psource-Pdissipation)dt + sum(Eplus-Eminus).
```

### 4. 合成波包层

`otr_synthetic_wavepacket`、`fig05_ideal_time_reflection` 和 `demo_circuit_time_interface` 用非对称波包显示：大/小特征的接收顺序反转、反射负极性、空间波数守恒和整体频移。`fig06`、`fig07`、`fig08`、`fig11` 则分别给出时间板干涉、有限上升沿/不同步趋势、反向时间板/视频泄漏和宽带波包预测。它们是机制复现，不是缺失实验数据的替代品。

### 5. 为什么 FDTD 只能作基准，尤其不能用于 OFF 边界

仓库当前实际存在的 [`tmm`](../../tmm/) 和 [`fdtd`](../../fdtd/) 公共代码只用于交叉检查传播方向或连续介质极限；PWE 类色散检查由本目录的离散 Bloch 求解器完成。README 不链接一个当前不存在的 `pwe/` 目录，本目录的论文复现也不依赖公共模块生成结论。

特别是不能直接用 [`fdtd/fdtd1d_db.m`](../../fdtd/fdtd1d_db.m) 中常见的 `D/B` 连续时间边界描述本文的 **OFF** 事件。本文在 OFF 时把已充电的并联电容从线路系统切除：电压/电场 `v/E` 连续，而线路所见电荷/位移 `q/D` 不连续，即 S15。强行套用普适 `D/B` 守恒会得到错误的 OFF 反射系数。需要 OFF 物理时应使用 `otr_temporal_interface(...,'voltage')` 或离散 MNA；FDTD 仅是基准，不是证据层。

## 主文逐图映射

| 主文图 | 可复现现象 | MATLAB 文件 | 证据级别与限制 |
|---|---|---|---|
| Fig. 1a | 均匀时间开关产生前向折射与后向时间反射 | `fig05_ideal_time_reflection`、`demo_circuit_time_interface` | 合成解析 + 独立 MNA；不重画原论文艺术示意图 |
| Fig. 1b | 制作的 30 单元 TLM 照片 | 无 | 照片不是数值现象；本地 PDF 中可查看 |
| Fig. 1c/d | 空间反射与时间反射、时间序反转、动量守恒 | `fig05_ideal_time_reflection` | 时间界面机制；空间界面只作为概念对照 |
| Fig. 1e | ON/OFF 实测与 ADS 色散、群速度、固定 `beta` 频移 | `fig01_source_data` | 官方 Source Data 复画 |
| Fig. 1f | 非对称脉冲的 TR 顺序反转、负极性、TT | `fig01_source_data`；`demo_circuit_time_interface` | 前者为官方示波器/ADS 数据，后者为独立 Kirchhoff/MNA 机制复现 |
| Fig. 2a | OFF→ON 的 `T/R` 幅度与相位、负 TR | `fig10_scattering_spectra` | S1–S3、S16 理论；ZIP 的 `2a/` 为空，不能复原实验阴影 |
| Fig. 2b | ON→OFF 的新型电压连续边界及 `T/R` | `fig10_scattering_spectra` | S1–S3、S17 理论；ZIP 的 `2b/` 为空 |
| Fig. 2c | 约 60 MHz 到约 34 MHz 的宽带红移 | `fig02_source_data` | 官方示波器数据、显式门控和 FFT |
| Fig. 2d | 约 33.6 MHz 到约 50 MHz 的宽带蓝移 | `fig02_source_data` | 官方示波器数据、显式门控和 FFT |
| Fig. 2e/f | 30–70 MHz 红移扫频和 20–34 MHz 反向蓝移 | `fig02_source_data`、`otr_frequency_map` | 官方频率表 + 守恒实波数理论 |
| Fig. 3a/b | 空间板无限散射与时间板仅四条因果路径 | `otr_temporal_slab`、`fig06_temporal_slab_theory` | 解析机制 |
| Fig. 3c | `tau=15/25/35 ns` 的两次 TR、间隔随 `tau` 改变 | `fig03_source_data` | 官方原始示波器波形，只去稳健 DC 偏置并归一化逻辑通道 |
| Fig. 3d | 时间 Fabry–Pérot 反射零点随 `tau` 移动 | `fig06_temporal_slab_theory` | 方程预测；官方 `3cd_raw` 不足以唯一重建论文的门控/补偿光谱 |
| Fig. 3e | 固定 `k=1.8 rad/m`、`f1=28.66 MHz` 下随 `tau` 连续调谐 | `fig11_slab_wavepackets` 生成的 `fig11b_slab_duration_sweep` | 波包平均后的预测；ZIP 的 `3e/` 为空，无五次实验均值和极差 |

## Extended Data / SI 逐现象映射

主文 Extended Data Fig. 1–9 与 SI Fig. S1–S9 内容对应；下表同时列出 SI 章节。

| Extended/SI | 现象 | MATLAB 文件 | 可复现程度 |
|---|---|---|---|
| ED1 / Fig. S1 / S1 | 单元 T 网络、S1–S3、色散、Bloch 阻抗、设计参数图 | `fig09_tlm_design`、`otr_abcd_unit_cell`、`otr_bloch_dispersion` | 理想模型数值复现；实际器件可选 `realistic` |
| ED2 / Fig. S2 / S2 | PCB 布局、完整接线与测量装置 | 无数值图；MNA 只实现其电气抽象 | PCB/Gerber、照片和完整仪器连接不可由方程生成 |
| SI S3 | T 接头和阻抗失配补偿 S4–S5 | `fig01_source_data` 仅显示官方已给波形 | 缺少频率相关连接器实测数据，不能唯一重做全补偿 |
| ED3 / Fig. S3 / S5 | 接电容时电荷连续、拆电容时电压连续，频率减半/加倍 | `fig04_boundary_conditions`、`otr_temporal_interface`、MNA | S10–S17 的初值结论和边界严格检验 |
| ED4 / Fig. S4 / S6 | 60 MHz 波包约 0.55 倍下变频 | `fig02_source_data` | 官方示波器数据复画 |
| ED5 / Fig. S5 / S7 | 三块级联样品、`tau=15 ns`、17.5 MHz FWHM 波包、约 38 MHz 共用反射零点 | `fig11_slab_wavepackets` | 解析包络预测；不是作者 ADS 原始工程 |
| ED6 / Fig. S6 / S8 | 开关不同步仍近似空间均匀、动量近似守恒 | `fig07_finite_switching_and_homogeneity` | 相位相干度趋势；不是未公开 ADS 波形 |
| ED7 / Fig. S7 / S9–S10 | 单元切换速度、3 ns 实验锚点、有限上升沿抑制 TR | `fig07_finite_switching_and_homogeneity`；MNA 可设 `riseTime` | 一阶 Fourier-overlap 趋势 + 论文锚点；8/12/16 ns 冲突按前述处理 |
| ED8 / Fig. S8 / S11 | ON–OFF–ON 反向时间板、第二路径损耗、零点不完全 | `fig08_inverted_slab_and_leakage` | 含 `A=0.5 dB/m` 的方程预测；无公开实测原始数据 |
| ED9 / Fig. S9 / S12 | 视频泄漏及 `V(t)-V(t-T)` 补偿 | `fig08_inverted_slab_and_leakage` | 合成泄漏演示；无公开泄漏参考采集 |

**Table S1 / Extended Data Table 1 边界：** 本地 PDF 中表格作为图像嵌入；其文本层只能提取标题，不能可靠提取元件行。代码只采用 Methods 正文明确给出的 82 pF、1 Ω、4 Ω、8 nH、2 pF、1 pF 等参数，不臆造 Table S1 的器件型号、封装或数量。完整表格可直接在本地 PDF 最后一页/SI 相应页人工查看。

## 文件索引

### 入口、数据处理与验收

| 文件 | 作用 |
|---|---|
| [`run_all_reproductions.m`](./run_all_reproductions.m) | 全量/快速一键运行，逐任务记录 PASS/FAIL，生成 `run_manifest.mat` |
| [`otr_setup.m`](./otr_setup.m) | 从任意工作目录定位复现目录、官方数据和 `output/`；若存在则运行仓库 `startup_stm.m` |
| [`otr_parameters.m`](./otr_parameters.m) | 论文参数、冲突后采用值及派生量的唯一默认入口 |
| [`otr_read_numeric_table.m`](./otr_read_numeric_table.m) | 无工具箱读取官方数值文本/CSV |
| [`otr_read_scope_csv.m`](./otr_read_scope_csv.m) | 解析 Tektronix 示波器 CSV 的时间和通道 |
| [`otr_gate_signal.m`](./otr_gate_signal.m) | 显式端点基线、余弦渐消和门函数，门限随结果保存 |
| [`otr_fft_spectrum.m`](./otr_fft_spectrum.m) | 校准单边/双边 FFT，并保留负频率检查 |
| [`otr_save_figure.m`](./otr_save_figure.m) | 将图输出到本目录 `output/` |
| [`otr_validate.m`](./otr_validate.m) | 对 S1–S3、频率映射、S14–S17、时间板、FFT、官方数据和 MNA 做确定性验收 |

### 物理核心

| 文件 | 作用 |
|---|---|
| [`otr_abcd_unit_cell.m`](./otr_abcd_unit_cell.m) | S1 理想/实际串联 RLC 单元 ABCD |
| [`otr_bloch_dispersion.m`](./otr_bloch_dispersion.m) | S2/S3 Bloch 波数、阻抗、通带与分支选择 |
| [`otr_frequency_map.m`](./otr_frequency_map.m) | 固定实 `beta` 的 OFF↔ON 宽带频率映射 |
| [`otr_temporal_interface.m`](./otr_temporal_interface.m) | S14/S16 电荷连续 ON 与 S15/S17 电压连续 OFF 的散射矩阵 |
| [`otr_temporal_slab.m`](./otr_temporal_slab.m) | OFF–ON–OFF 或反向时间板及四条路径 |
| [`otr_synthetic_wavepacket.m`](./otr_synthetic_wavepacket.m) | 空间均匀理想切换的合成解析波包 |
| [`otr_build_mna.m`](./otr_build_mna.m) | 构造离散 Kirchhoff/MNA 线路、端口、互电容和开关支路 |
| [`otr_switch_schedule.m`](./otr_switch_schedule.m) | OFF–ON、ON–OFF、时间板/反向板、理想或平滑有限边沿调度 |
| [`otr_simulate_circuit.m`](./otr_simulate_circuit.m) | 梯形法求解 `C(w)xdot+G(w)x=b(t)` 并执行/记录事件边界 |
| [`otr_circuit_energy.m`](./otr_circuit_energy.m) | 储能、损耗、源功、开关功和离散能量平衡 |

### 出图入口

| 命令 | 主要输出/现象 |
|---|---|
| `fig01_source_data` | Fig. 1e/f 官方色散、群速度与波形 |
| `fig02_source_data` | Fig. 2c–f 官方红移/蓝移与扫频 |
| `fig03_source_data` | Fig. 3c/d 归档中的原始 `tau=15/25/35 ns` 示波器波形 |
| `fig04_boundary_conditions` | ED3/SI S3 的 ON 电荷连续与 OFF 电压连续 |
| `fig05_ideal_time_reflection` | 时序反转、负极性和波数守恒 |
| `fig06_temporal_slab_theory` | 四路径时间 Fabry–Pérot、反射零点、固定 `k` 调谐 |
| `fig07_finite_switching_and_homogeneity` | 有限开关时间和空间不同步趋势 |
| `fig08_inverted_slab_and_leakage` | 反向时间板损耗和视频泄漏补偿预测 |
| `fig09_tlm_design` | ED1/SI S1 单元设计与参数图 |
| `fig10_scattering_spectra` | 主 Fig. 2a/b 的理论散射幅相 |
| `fig11_slab_wavepackets` | ED5/SI S5 宽带包络和主 Fig. 3e 持续时间扫描预测 |
| `demo_circuit_time_interface` | 30 单元 Kirchhoff/MNA 时空图、固定 `k` 频移、事件/能量自检 |

任一函数均可单独运行并返回可审计结果结构，例如：

```matlab
r2 = fig02_source_data();
rc = demo_circuit_time_interface();
```

## 输出目录

所有新结果写入 [`output/`](./output/)；不会修改仓库公共 `pwe/tmm/fdtd` 文件。主要产物为：

- `fig01_source_data` 至 `fig11_slab_wavepackets`：同名 `.png` 和结果 `.mat`；
- `fig11b_slab_duration_sweep.png`：固定波数下的时间板持续时间扫描；
- `circuit_time_interface.png/.fig`：MNA 时空场、探针时序与能量图；
- `circuit_fixed_k_frequency_translation.png/.fig`：离散梯形线路的固定 `k` 频移；
- `demo_circuit_time_interface.mat`：MNA 关键时序、色散和验收指标；
- `validation_report.mat` 与 `validation_summary.txt`：机器可读和文本验收报告；
- `run_manifest.mat`：一键运行后每项任务的状态、用时和错误消息。

再次运行会更新同名输出。`.mat` 中保留处理门、频率轴、系数和限制说明，便于追踪“画出来的线”来自官方数据还是模型。

## 定量验收

### 自动物理/数值不变量

`otr_validate` 的硬性条件包括：

- S1 元素误差 `<1e-13`，理想无源单元行列式误差 `<1e-12`；
- S2 恒等式误差 `<1e-12`，50 MHz 时 `ZB,on` 距 25 Ω 小于 0.1 Ω；
- 理想和实际固定 `beta` 映射残差 `<1e-10 rad/m`，往返频率误差 `<1 Hz`；
- S16/S17 四个系数及符号误差 `<1e-14`；
- 时间板四路径分解误差 `<1e-12`，首个反射零点数值/解析差小于一个频率网格；
- FFT 能正确给出 `+50 MHz` 实信号峰和 `-70 MHz` 复信号峰；
- 官方扫频表尺寸、单调性和有限数值通过检查；
- MNA ON 电荷边界、OFF 电压边界、事件残差和能量相对残差分别满足源码中的严格阈值，能量阈值为 `<5e-8`。

当前 `validation_summary.txt` 的代表性结果为：

| 检查 | 当前结果 |
|---|---:|
| S1 最大元素误差 / 行列式误差 | `1.11e-16 / 1.11e-16` |
| S2 最大误差 | `1.80e-16` |
| `ZB,on(50 MHz)` | `24.94896 Ω` |
| 理想 OFF→ON 映射，输入 `[30,50,70] MHz` | `[16.74758,27.77041,38.57156] MHz` |
| 实际串联 RLC 映射 | `[17.09726,28.25621,39.05636] MHz` |
| `[Ton,Ron,Toff,Roff]` | `[0.375,-0.125,1.5,-0.5]` |
| `tau=15 ns` 理想板首个零点 | 理论 `33.33333 MHz`，数值 `33.33500 MHz` |
| 小型 MNA 能量相对残差 | `9.40e-14` |

### 官方波形的可核对数值

本目录现有 `fig02_source_data.mat` 使用显式门控得到：

| 情形 | 入射峰 | TR 峰 | 时间折射峰 | 论文报告 |
|---|---:|---:|---:|---:|
| OFF→ON 红移 | 59.797 MHz | 34.377 MHz | 33.498 MHz | 60.0 → 34.5 / 33.6 MHz |
| ON→OFF 蓝移 | 33.484 MHz | 49.587 MHz | 50.203 MHz | 33.6 → 49.5 / 50.1 MHz |

对官方 CSV 自带理论列，红移扫频的 `[TR,TT]` RMSE 为 `[2.33,2.01] MHz`，蓝移扫频为 `[3.84,3.76] MHz`。这些是公开测量点相对理论的偏差，不应通过隐藏拟合抹掉。

`fig01_source_data` 在 10–70 MHz 线性段得到 OFF/ON 实测群速度约 `100.48/58.11 Mm/s`，对应 ADS 约 `106.67/56.74 Mm/s`。由两条实测色散做固定 `beta` 检查时，60 MHz OFF 映射到约 36.39 MHz ON；这与 Fig. 2 波包门控峰值不是同一种估计器，二者不要求逐点相等。

30 单元 MNA 演示当前给出：

- 61 个基础状态（31 节点电压 + 30 线路电流）；
- 开关事件 `105 ns`，落网格误差为 0；
- Methods 参数所对应的离散低 `k` 频率比 `0.5732`；它不同于理想设计目标 0.5，是 `epsEff` 推导的线路电容与实际 82 pF 负载共同作用的结果；
- 非对称探针波形的带符号时间反转相关系数 `-0.9140`；
- 能量平衡相对残差 `2.92e-14`，ON 电荷残差为 0，OFF 电压跳变为 0 V。

## 无法由公开材料唯一复现的边界

下列内容缺少必要输入，因此本目录有意不声称“实验级精确复现”：

- Fig. 1b/ED2 的实验照片、PCB 铜皮/过孔几何、Gerber、完整 BOM 与器件批次公差；
- Table S1 的机器可读元件行；本地 PDF 文本层只有表题；
- 作者的 Keysight ADS 工程、开关宏模型、瞬态/S 参数求解设置；公开数据只含部分 ADS 导出曲线；
- Fig. 2a/b 的实验散射系数数组；官方 ZIP 中两个目录为空；
- 连接器完整频率相关 S 参数、每次实验的有效传播距离 `xr/xt`、ON/OFF 损耗曲线与相位校正，因此不能唯一重算 S4–S9 的所有补偿；
- Fig. 3d 所需的唯一入射/TR 门、独立视频泄漏参考周期、连接器与传播损耗补偿；`3cd_raw` 只能支持原始波形复画；
- Fig. 3e 的五次测量、均值、最大/最小范围；官方 `3e/` 为空；
- ED6–ED9/SI S6–S9 的 ADS/示波器原始数组；当前图只验证文中趋势和方程；
- 实际各开关的延时分布、上升/下降沿不对称、内部接地放电细节、测量噪声与 512 次平均的单次记录。

因此，本目录的验收目标是：官方数据面板数值可追溯、S1–S17 的可执行部分满足方程、离散 MNA 满足 Kirchhoff 边界和能量平衡、缺失数据面板只做明确标注的预测；不是对论文排版、照片、噪声或未公开 ADS 波形的像素级仿制。

## 与仓库公共代码的关系

仓库当前已有 `tmm/`、`fdtd/`，而本工作区没有独立的 `pwe/` 目录；PWE/Bloch 类周期色散功能在 `otr_abcd_unit_cell` 与 `otr_bloch_dispersion` 中自包含实现。这些公共代码可用于交叉验证转移矩阵或传播极限，但它们不属于 Nature 官方 Source Data，也不代替本目录的微观开关边界。最终论文复现应以本目录的官方数据脚本、S1–S17 解析函数和 Kirchhoff/MNA 电路为主；公共代码的结果只能作为独立基准。
