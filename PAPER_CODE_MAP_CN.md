# 20 份文献—代码覆盖矩阵

状态含义：

- **复现**：代码直接对应论文主要数值模型或目标图；
- **覆盖**：核心方程已经实现，但参数/平台不是逐图复现；
- **部分**：只覆盖理想化子问题，不能用来声称复现实验；
- **独立分支**：需要新的状态变量或数值框架，不应硬接到现有求解器。

库中的 20 份 PDF 对应 17 项独立工作。第 1、11、19 项是一篇工作的预印本与重复正式版本；第 3、15 项是预印本与正式发表版本。

| # | 库中文件/工作 | 主要模型或科研问题 | v2 对应 | 状态与下一步 |
| ---: | --- | --- | --- | --- |
| 1 | `1803.08731v1.pdf`，Topology of photonic time-crystals | 时间 Zak 相、单胞次序、时间畴壁态 | `demo09`、`demo10` | **部分**：已做有限响应与畴壁本征态匹配；尚需严格时间 Zak Wilson 环 |
| 2 | 2015 APL dynamic transmission line | 离散加载传输线中的真实 $k$ 带隙 | `demo03` 仅给连续介质基线 | **独立分支**：建立 LC 单胞、Bloch 电路矩阵和实验损耗模型 |
| 3 | `2507.02223v1.pdf` | 非合成 PTC 放大、相位和时间拓扑态实验 | `demo09`、`demo10` | **部分**：理想二元 PTC；尚未包含实验波形、线路色散和校准 |
| 4 | Broadband coherent wave control through photonic collisions | 时间界面的双输入相干控制 | `demo11_coherent_time_interface` | **覆盖**：理想全局时间界面；可继续加入实验开关律和脉冲带宽 |
| 5 | Conformal Spatiotemporal Modulation Enabled Geometric Frequency Combs | 共形缩放、移动边界、几何频率梳 | 无 | **独立分支**：需要事件驱动移动界面/尺度映射，周期 ST-PWE 不适用 |
| 6 | Dirac mass induced by optical gain and loss | 合成网格、非 Hermitian Dirac 质量 | `demo08` 只提供通用 FHS 算法 | **独立分支**：实现论文的离散演化映射与局域几何相 |
| 7 | Externally Modulated Low-Pass Transmission Line | 时变低通传输线的 Floquet 电路传播 | 无 | **独立分支**：实现离散电报方程、时变电容与谐波平衡 |
| 8 | Metasurface-based realization of photonic time crystals | 时变超表面、表面/辐射 Floquet 谐波 | 无 | **独立分支**：需要表面阻抗边界和开放辐射通道 |
| 9 | Observation of temporal reflection and broadband frequency translation | 电路时间界面、频率平移、不同开关微观条件 | `demo05`、`demo11`、`temporal_interface_matrix_jump` | **部分**：已显式区分跃迁律；实验电路的电荷/电压约束仍需逐装置建模 |
| 10 | Observation of temporal topological boundary states in a momentum bandgap | 两光纤环合成网格、动量 Dirac 方程、时间边界态 | 无 | **独立分支**：不能直接套连续介质 Zak 相；应实现论文离散步进方程 |
| 11 | `optica-5-11-1390.pdf` | 与第 17 项相同的正式论文 | 同第 1 项 | **重复文件** |
| 12 | Theory and applications of photonic time crystals: a tutorial | PWE/TMM、色散、有限尺寸、缺陷、无序、非线性、平台与应用 | 全包；尤其 `demo03`、`demo09`–`demo12` | **基础章节覆盖**：色散、空间有限、各向异性、无序、非线性仍是研究扩展 |
| 13 | PRL 128, 186802, Topological Space-Time Crystal | 量子紧束缚 Sambe 空间、时空平移对称 | 无 | **独立分支**：这是量子/有效晶格模型，不是标量 Maxwell PTC |
| 14 | PRX 4, 021017, Surface Impedance and Bulk Band Geometric Phases | 静态一维 PC 的 Zak 相和反射相 | `demo07_zak_phase_stpwe` 为方法延伸 | **部分**：未逐图复现表面阻抗；可作为空间拓扑基准 |
| 15 | `s41467-025-66154-4.pdf` | 第 3 项的正式发表版本 | 同第 3 项 | **预印本/正式版对应** |
| 16 | Spacetime-topological events | 非 Hermitian 离散量子行走、时空 winding、因果耦合 | 无 | **独立分支**：需实现四步 Floquet 映射和空间/时间/时空 winding |
| 17 | Spatiotemporal plane wave expansion method | 任意一维时空周期介质 ST-PWE、论文图 2 | `demo01`、`demo02`、`demo06` | **复现**：现包最完整的论文级复现 |
| 18 | Time-reflection of microwaves by a fast optically controlled time-boundary | 光控微波时间反射、传输线实验 | `demo05`、`demo11` | **部分**：已覆盖理想波分裂；尚缺线路色散、开关速度和测量链 |
| 19 | `Topological aspects of photonic time crystals.pdf` | 与第 11 项相同的正式 PDF | 同第 1 项 | **重复文件** |
| 20 | Temporal multilayer structures for designing higher-order transfer functions | 任意时间多层、透明与放大滤波器 | `demo04_temporal_multilayer_tmm` | **复现/覆盖**：主要设计类别与表格参数已实现 |

## 覆盖结论

现包的优势集中在同一物理层级：

- 一维、各向同性、局域、非色散 Maxwell 介质；
- 空间均匀的时间界面/时间晶体；
- 一维时空周期连续介质。

文献库后半部分出现了五种新的物理层级：

1. 真实传输线电路；
2. 时变超表面开放系统；
3. 合成光子网格/离散量子行走；
4. 共形移动界面；
5. 时间色散、非线性和有限空间尺寸。

这些应当建立为平行模块，并通过共同的“参数—求解—诊断—验证”接口连接，而不是把它们全部写成 `epsilon(x,t)` 传给同一个函数。
