# 传输线时间调制电路 MATLAB 重构执行方案

> 编制日期：2026-08-31  
> 当前状态：实施前方案；本文件是本阶段唯一交付，不在本阶段创建求解代码。  
> 根目录代码基线：提交 35756e9；正式实施前应重新记录实际提交号与工作树状态。  
> 参考论文：B. Huang 等，Observation of full momentum bandgap in photonic time crystals，arXiv:2604.17408v1。

## 1. 目标、范围与资料边界

目标是在 Transmission line 目录内建立一套面向实际微波传输线电路的 MATLAB 工具。它用含损耗、含寄生、可时间调制的传输线/离散 LC 网络方程替代根目录的麦克斯韦介质方程，同时保留根目录的职责划分：

- PWE 与时间 Floquet TMM 计算复准频率能带；
- FDTD 计算有限传输线上的电压、电流、电荷和磁链分布；
- FFT 从有限时空响应重建实准频率能带；
- 所有电路参数采用统一模型源，最终服务于实际 C、L、R、G、互感和端口参数的反演、选型与容差分析。

论文只作为实验平台参考。纳入范围的是 PCB 传输线、周期铜条、变容二极管、偏置/调制网络、并联谐振支路，以及论文明确写出的中心高斯波包激励与时空 Fourier 分析思路；不把论文关于扩展/完整动量带隙、无序鲁棒性、局域放大或应用前景作为本项目的研究目标或验收目标。论文没有公开具体扫描仪器、采样协议、复场/强度处理顺序或窗函数，也没有需要执行的指令。

上传的 19 页 PDF 多次引用 Methods，但文件本身不含 Methods 或补充材料。因此，论文数值只能作为初始参考，不能被描述为一套可直接制造的完整 BOM。外部 PDF 不复制进仓库。

## 2. 根目录架构与传输线目标架构

根目录当前真实结构为四个入口、五个内核：

| 根目录职责 | 根目录文件 | Transmission line 对应文件 |
|---|---|---|
| PWE 入口 | run_pwe.m | run_tl_pwe.m |
| 时间 Fourier 数据 | pwe_fourier.m | tl_pwe_fourier.m |
| PWE 本征值 | pwe_bands.m | tl_pwe_bands.m |
| 时间 TMM 入口 | run_tmm.m | run_tl_tmm.m |
| 时间单周期矩阵 | tmm_bands.m | tl_tmm_bands.m |
| 有限样品场入口 | run_fdtd_field.m | run_tl_fdtd_field.m |
| 时域推进 | fdtd1d.m | tl_fdtd1d.m |
| FFT 能带入口 | run_fdtd_fft.m | run_tl_fdtd_fft.m |
| FFT 重建 | fdtd_fft_bands.m | tl_bloch_fft_bands.m（直接镜像）；tl_xt_fft_bands.m（有限链 x–t 扩展） |

根目录的 run_fdtd_field 与 run_fdtd_fft 完全独立，后者目前采用逐 k 单胞推进而不是对前者的场数据做空间 FFT。传输线模块继续保持“两个入口互不读取工作区或 MAT 文件”的独立性，但作一项有意扩展：

- 主 FFT 路径采用“受论文测量结果启发”的有限链 x–t 二维 FFT 目标。run_tl_fdtd_fft 自己启动一套专用宽带 FDTD，随后调用 tl_xt_fft_bands；具体实验扫描流程须由后续测量协议确定；
- 另设 tl_bloch_fft_bands，以逐 k Bloch 单胞推进和同相位采样提供无限体基准。它用于验证与暗模诊断，不增加日常入口。

这既保持根目录的入口职责，又能直接对接未来逐单元时空测量。

## 3. 物理模型

### 3.1 麦克斯韦量与传输线量的对应

| 一维电磁量 | 传输线量 | 说明 |
|---|---|---|
| 电场 E | 节点电压 V | 主要电观测量 |
| 磁场 H | 支路电流 I | 与论文的局域磁场观测最接近 |
| 电位移 D | 电荷 Q=C V | 时间界面连续状态 |
| 磁感应 B | 磁链 Φ=L I | 时间界面连续状态 |
| 介电参数 ε | 并联电容 C 或 C′ | 可由变容管调制 |
| 磁导参数 μ | 串联电感 L 或 L′ | 包含铜条/耦合结构 |
| 电导损耗 | 并联漏导 G 或 G′ | 介质、偏置网络等损耗 |
| 磁/串联损耗 | 串联电阻 R 或 R′ | 铜损、电感 ESR 等 |

内部推进必须使用 Q/Φ，而不是直接推进 V/I。这样在理想时间突变处，自由电荷和外加冲击电压均不存在时，Q 与 Φ自然连续；V=Q/C 和 I=Φ/L 可以发生物理跃变。这对应根目录使用 D/B 而不是直接使用 E/H 的原则。

### 3.2 含损耗、含时变参数的电报方程

连续传输线写成

\[
-\frac{\partial V}{\partial x}
=R'I+\frac{\partial \Phi}{\partial t},
\qquad
-\frac{\partial I}{\partial x}
=G'V+\frac{\partial Q}{\partial t},
\]

\[
\Phi=L'I,\qquad Q=C'V .
\]

当 C 随时间变化时，

\[
\frac{\partial Q}{\partial t}
=C'\frac{\partial V}{\partial t}
+\frac{\partial C'}{\partial t}V .
\]

因此不得把并联电流错误地写成只有 C′∂V/∂t；采用 Q 状态可自动保留泵浦项。若 L 随时间变化，同理必须保留 ∂(L′I)/∂t。

连续参数与每个数值/物理单元长度 Δx 的集中参数关系为

\[
L_s=L'\Delta x,\quad R_s=R'\Delta x,\quad
C_p=C'\Delta x,\quad G_p=G'\Delta x .
\]

代码必须明确 Δx 是“数值剖分长度”还是“真实单胞长度 a”。真实 LC 链的正弦型晶格色散是物理效应，不能在网格收敛中被错误消除。

### 3.3 有限离散链的保守状态

对节点 n 和连接 n、n+1 的支路，基础更新方程为

\[
\dot q_n=i_{n-1/2}-i_{n+1/2}-G_pv_n-i_{0,n}+i_{\mathrm{src},n},
\]

\[
\dot\phi_{n+1/2}
=v_n-v_{n+1}-R_si_{n+1/2},
\]

\[
v_n=\frac{q_n}{C_{\mathrm{eff},n}(t)},\qquad
\boldsymbol\phi=\mathbf L_{\mathrm{series}}\boldsymbol i .
\]

第一版可用有效串联电感；第二版必须允许最近邻互感 S，使 Lseries 成为对称正定矩阵。有限链中预先分解该矩阵，避免每个时间步重复求逆。

### 3.4 两种论文实验拓扑

SSPP 单元：

- 相邻铜条等效为串联/耦合电感网络；
- 每个节点通过过孔、变容二极管和滤波/偏置寄生形成时变并联电容；
- 最小状态为节点电荷 q 与串联支路磁链 φ。

拓扑依据为论文第 5 页及 Fig. 2a–c（PDF 第 12 页）。

CROW 单元：

- 在 SSPP 单元上增加并联谐振支路；
- 理想近似为并联 L0；
- 完整模型中 L0 与 200 pF 隔直电容串联；本项目另加入待标定的电感 ESR 以描述损耗，论文没有给出其数值。

拓扑依据为论文第 7 页及 Fig. 3a–c（PDF 第 13 页）。

完整 CROW 辅助状态可写为

\[
\dot\phi_{0,n}=v_n-R_0i_{0,n}-v_{b,n},\qquad
\dot q_{b,n}=i_{0,n},
\]

\[
i_{0,n}=\frac{\phi_{0,n}}{L_0},\qquad
v_{b,n}=\frac{q_{b,n}}{C_{\mathrm{block}}}.
\]

这样无需借用论文中未给出的有限隔直电容修正式。

论文只给出了 k=0 和 k=π/a 的频率锚点，没有给出下列完整色散。作为本方案的 provisional model，在无损、理想隔直、最近邻 LC 链以及互感符号按下式定义时，暂定静态解析基准为

\[
\omega_{\mathrm{SSPP}}^2(k)=
\frac{4\sin^2(ka/2)}
{C_0\,[L_s+2S\cos(ka)]},
\]

\[
\omega_{\mathrm{CROW}}^2(k)=
\frac{1}{L_0C_0}+\omega_{\mathrm{SSPP}}^2(k).
\]

在 k=π/a 处得到论文给出的 Ls−2S 组合。必须检查 Ls>2|S|，确保电感矩阵正定。有限 Cblock、更多互感、端口和寄生由完整状态模型求解，不强行套用上述简式。

### 3.5 变容二极管与调制波形

初始线性模型采用外部规定的正周期电容：

\[
C_{\mathrm{eff}}(t)=C_0+\Delta C\cos(2\pi f_mt+\varphi).
\]

进入器件标定阶段后改为

\[
C_{\mathrm{eff}}(t)
=C_{\mathrm{par}}
+C_{\mathrm{var}}\!\left[
V_{\mathrm{DC}}+
V_m|H_{\mathrm{BF}}(f_m)|
\cos(2\pi f_mt+\varphi+\angle H_{\mathrm{BF}})
\right],
\]

并由实测 C–V 表、封装寄生和调制链频响生成一个周期的数值波形及其谐波。PWE/TMM/FDTD 使用同一波形回调。

第一阶段只处理“由泵浦预先规定、与信号幅值无关”的线性时变电容。信号引起的变容管非线性、自泵浦耗尽和饱和属于后续非线性扩展，不能与线性 Floquet 能带混用。

## 4. 论文实验参数的可追踪种子

代码内部一律保存 SI 值；MHz 文本被解释为普通频率 f，角频率必须显式转换为 Ω=2πf。论文虽使用 Ω/ω 符号，却在正文和图中用 MHz 标注，不能省略 2π。下表页码同时对应该 PDF 的页码和页面印刷号。

| 参数 | SSPP 初值 | CROW 初值 | 来源位置 | status/备注 |
|---|---:|---:|---|---|
| 论文标记的铜条间距 w | 4 mm | 12 mm | 第 5、7 页 | paper-given；不能直接当作 a |
| 直流偏置 | 14 V | 10.8 V | 第 5、7 页 | paper-given |
| 偏置点变容值 Cvar | 12 pF | 17 pF | 第 5、7 页 | paper-given |
| 有效静态电容 C0 | 19.8 pF | 27.8 pF | 第 5、7 页 | paper-given/effective |
| 调制幅值 ΔC | 2.38 pF | 5 pF | 第 5、8 页 | paper-given |
| 相对调制深度 ΔC/C0 | 0.1202 | 0.1799 | 由上两行计算 | paper-derived |
| 调制频率 fm | 575/675/725 MHz | 700 MHz | 第 5–6、8 页 | paper-given |
| 并联电感 L0 | 无 | 9.2 nH | 第 7 页 | paper-given |
| L0 支路隔直电容 | 无 | 200 pF | 第 7 页 | paper-given |
| 静态频率锚点 | fr0=371 MHz | fcol=315 MHz | 第 5、8 页 | paper-given |
| 耦合电感参考 | 未给 | Ls≈165 nH | 第 8 页定义及 Fig. 3e（第 13 页） | paper-derived/reference；非样品 BOM |
| 论文定义的 Q0 | 不使用 | 20.7 | Fig. 3e（第 13 页） | paper-given/reference；非常规损耗 Q |

论文把 C0 称为包含变容管、滤波网络影响和放大器输出电阻影响的“有效电容”。因此 C0−Cvar 不能未经标定就解释成一只实际并联电容，C0 也不能替代器件级 C–V/阻抗模型。

由 SSPP 的 C0 和 fr0 可得到

\[
L_s-2S=
\frac{4}{(2\pi f_{r0})^2C_0}
\approx37.18\ {\rm nH}.
\]

这只确定一个组合参数，不能分别确定 Ls 与 S。CROW 的 L0=9.2 nH、C0=27.8 pF 给出 fcol≈314.7 MHz，可作为单位和公式的首个回归。

以下量在上传 PDF 中缺失，必须标为待标定，禁止赋予看似精确的默认值：

- 晶格常数 a；论文没有明确说明 a=w；
- SSPP 的 Ls、S，CROW 样品的实际 Ls、S，以及更远邻互感；
- 铜损、变容管 ESR、电感 ESR、并联漏导和偏置网络损耗；
- 源阻抗、负载阻抗、端口耦合和探针传递函数；
- 变容管型号、C–V/ESR 曲线、封装寄生和所需泵浦电压；
- BF、放大器、隔直/偏置网络的完整参数与相位延迟；
- FR4 厚度、介电常数、损耗角正切、铜厚和完整几何；
- 样品确切单元数、总长度和照片中各单元接口的用途；
- 激励中心频率、带宽、绝对幅值、测量仪器链、采样率和空间扫描细节；
- 泵浦电压幅值、沿各单元的调制幅相与空间同步性。

论文定义 Q=ωcol C0 Z0，其中 Z0 指自由空间阻抗。它不是由串联电阻得到的常规谐振器品质因数，不能用于反推 R，也不能把该 Z0 当作实验端口的 50 Ω。

## 5. 统一约定与模型接口

### 5.1 数学与单位约定

- 内部单位：H、F、Ω、S、m、s、rad/s；
- 相位：exp(ikx−iωt)；
- Floquet 乘子：λ=exp(−iωT)，ω=i log(λ)/T；
- 时间第一 Floquet 区：Re(ω)/Ω∈[−1/2,1/2)；
- 空间第一布里渊区：ka∈[−π,π)；
- Im(ω)>0 表示净增长，Im(ω)<0 表示净衰减；
- 图上同时给出工程单位 MHz 和无量纲 ω/Ω、ka/π；
- 所有输入字段名带物理含义，禁止用同一字段混放 Hz 与 rad/s。

### 5.2 单一参数源

tl_build_model.m 是四条流水线唯一的电路参数源。建议返回：

~~~text
model.kind                              % 'sspp' 或 'crow'
model.cell.a
model.series.Ls, Rs, mutualS
model.shunt.C0, Gp, Cpar
model.resonator.L0, R0, Cblock
model.bulk                              % 严格周期单胞，供 PWE/TMM/Bloch
model.finite                            % 每单元容差、无序、泵浦误差和端口
model.modulation.fmHz, OmegaRadPerSec
model.modulation.bulkFcn(t)
model.modulation.finiteFcn(cellIndex,t)
model.ports.Zsource, Zload
model.provenance.<parameter>            % value/source/page/status/uncertainty
model.derived                           % fcol、边界频率、Zref、稳定性上界等
~~~

来源状态固定为 paper-given、paper-derived、figure-estimate、datasheet、measured、fit 或 assumption。论文未报告的不确定度必须写成 unknown/not-reported，禁止伪造数值。任何 assumption 必须由入口脚本显式开启；不允许核心函数静默补齐缺失硬件参数。builder 返回不可变的配置快照；拟合与选型函数只通过显式 overrides 重新调用 builder，不能旁路这一参数源。

### 5.3 公共 API 与结果形状

实施前冻结以下调用关系：

~~~matlab
model     = tl_build_model(physicalCfg);
fourier   = tl_pwe_fourier(model,pweCfg);
pwe       = tl_pwe_bands(fourier,model,kScan,pweCfg);
tmm       = tl_tmm_bands(model,kScan,tmmCfg);
field     = tl_fdtd1d(model,fdtdCfg);
fftBands  = tl_xt_fft_bands(observation,referenceObservation,fftCfg);
bulkBands = tl_bloch_fft_bands(model,kScan,blochCfg);
~~~

- kScan 一律为 1×Nk；复频和模态权重一律为 Nmode×Nk；原始本征矢另带 harmonic/state 维；
- PWE/TMM 必须返回 k、omegaRaw、omegaFolded、omegaSelected、selectedRawIndices、replicaOrder、selectionConfidence、模式条件数和配置快照；
- FDTD 必须显式区分节点与支路的 Yee 错位：

~~~text
field.node.x, field.node.t, field.node.V, field.node.Q
field.branch.x, field.branch.tHalf
field.branch.IHalf, field.branch.PhiHalf
field.branch.IAtNodeTime          % 采用文档化中心化后才允许返回
~~~

- observation 是一个带 data、x、t、cellIndex、channelId、channelOffset 和 observableName 的结构，禁止把节点坐标/整数时间与支路坐标/半时间数据混传；
- referenceObservation 来自同一源和端口、关闭调制或移除周期负载的独立参考运行，用于构造 source-support mask；FFT 内核不得隐式调用 PWE/TMM 生成掩码。

## 6. 目标文件结构

~~~text
Transmission line/
├── TRANSMISSION_LINE_IMPLEMENTATION_PLAN.md   # 本阶段交付
├── README.md                                  # 后续：模型、单位、入口和限制
├── tl_build_model.m                           # 统一电路与来源追踪
├── run_tl_pwe.m
├── tl_pwe_fourier.m
├── tl_pwe_bands.m
├── run_tl_tmm.m
├── tl_tmm_bands.m
├── run_tl_fdtd_field.m
├── tl_fdtd1d.m
├── run_tl_fdtd_fft.m
├── tl_xt_fft_bands.m                          # 主：有限链 x–t FFT
├── tl_bloch_fft_bands.m                       # 辅：逐 k 无限体 FFT
├── tl_fit_parameters.m                        # 静态/损耗/动态参数反演
├── tl_select_components.m                     # 标称值、寄生与容差回算
└── validate_tl_suite.m                        # Base MATLAB 一键验收
~~~

文件保持扁平，全部加 tl_ 前缀，避免 MATLAB 路径上与根目录同名函数相互遮蔽。从 Transmission line 目录运行，不使用隐藏 addpath、startup 或根目录数值内核依赖。生成的 MAT、FIG、PNG 和拟合中间结果只写入已忽略的 output/。

## 7. 各数值模块实施要求

### 7.1 PWE

将每个 Bloch 波数下的保守状态写成

\[
\dot{\boldsymbol s}=\mathbf A(k,t)\boldsymbol s .
\]

采用

\[
\mathbf A(k,t)=\sum_p\mathbf A_p(k)e^{-ip\Omega t},\qquad
\boldsymbol s(t)=e^{-i\omega t}\sum_m\boldsymbol s_m e^{-im\Omega t},
\]

建立

\[
\sum_n\left[
i\mathbf A_{m-n}(k)-m\Omega\delta_{mn}\mathbf I
\right]\boldsymbol s_n
=\omega\boldsymbol s_m
\]

的 block-Toeplitz Floquet 本征问题。实施要点：

1. tl_pwe_fourier 只在无重复端点网格上采样与 k 无关的 1/C、1/L、损耗和谐振器周期系数，覆盖场截断所需的全部差阶 [−2M,2M]；tl_pwe_bands 再用单胞 incidence/Bloch 矩阵为每个 k 组装 Ap(k)。Fourier 层不得在没有 kScan 的情况下暗中采样完整 A(k,t)。
2. 保留 Nt≥4M+1 的最低检查；方波/实测尖锐波形必须同时对 M 和 Nt 加密。
3. R/G 会使问题非 Hermitian；不得照搬根 PWE 对无损 ε 卷积矩阵的 Hermitian 正定结论。应改为检查 C、L 正值、R/G 非负和电感矩阵正定。
4. SSPP、理想 CROW、有限 Cblock CROW 的状态维数不同，不能硬编码“每个 k 只有两条带”。保留 raw 谱、折叠谱、中心谐波权重和 V/I 可观测权重。
5. 无损普通点可用模态重叠追踪；含 R/G 的非 Hermitian 问题使用左右本征矢的双正交重叠、频率邻近和 invariant-subspace 共同追踪，并输出条件数。简并或 exceptional point 处允许分支身份未定义，禁止强行连续编号。
6. 零频环流或端口暗模保留在 raw 输出中，但不伪装成可观测电压带。去除 Floquet 副本后另返回 omegaSelected、selectedRawIndices、replicaOrder 和 selectionConfidence，供入口绘图及与有限维 TMM 公平比较。
7. 输出还包含 k、omegaRaw、omegaFolded、左右本征矢、中心谐波/V/I 权重、截断/采样元数据和完整模型快照。

### 7.2 时间 Floquet TMM

这里的 TMM 是“固定 Bloch k 的时间单周期状态转移矩阵”，不是工程上沿空间级联的 ABCD 矩阵。

1. 对方波或分段常数电容，在每层构造 Aj(k)，计算 Pj=expm(Aj Δtj)。
2. 单周期矩阵按真实时间顺序左乘，U=PN…P2P1；Q/Φ 状态跨时间界面无需额外跃迁矩阵。
3. 对正弦或实测 C(t)，用中点分层，并把时间切片数加倍直至复准频率收敛。
4. 由 eig(U) 得 λ，再用 ω=i log(λ)/T 折叠；强增长/衰减时同时监控 U 的条件数和乘子动态范围。
5. 含 R/G 时不假定 det(U)=1，也不强制增长/衰减共轭。无调制且无损时检查能量守恒；有调制且无损时检查辛/体积结构、乘子互易性以及 ΔE=Wpump；一般情况使用完整源—泵浦—损耗—端口能量账本。
6. 跨 k 使用与 PWE 一致的左右模态/子空间匹配；若只把每个 k 当作无序本征值集合，也必须在输出中声明，不得伪造分支编号。
7. 输出 k、U、lambda、omegaRaw、omegaFolded、本征矢、条件数、切片收敛信息与输入快照。

TMM 先于 PWE 完成，作为任意周期波形离散化与 PWE 截断收敛的独立基准。

### 7.3 有限链 FDTD

tl_fdtd1d 在真实离散 LC 链上推进 q/φ，并恢复 V/I：

1. 无损 SSPP 是最小可验证内核；随后依次加入 Rs/Gp、互感 S、CROW 的 L0/R0、有限 Cblock 和端口。
2. 电感/电阻损耗采用中心或梯形处理，避免显式欧拉产生伪增长；泵浦切换必须落在整数时间节点。
3. 稳定性由完整电路最大本征频率控制，包括整个泵浦周期的最小 C、完整互感矩阵、Cblock 和寄生高频极点。无损 staggered/leapfrog 基准强制 dt·ωmax≤1.8（严格稳定边界为 2）并用 dt/2、dt/4 验证；若改用其他积分器，必须先登记其离散放大矩阵和稳定域。物理 Floquet 增长不能被误判为数值不稳定；未来若调制 L，还要单独规定 L 界面与半时间格的对齐。
4. 边界支持 matched、periodic、open、short；实验主算例显式给出 Thevenin/Norton 源、源内阻和负载，不把根目录 sponge 当成 50 Ω 端接。
5. 支持复高斯窄带波包、宽带脉冲和可移动单元源；源位置、调制启动相位和端口参考面进入结果元数据。
6. 按 5.3 节的 node/branch 两套坐标和时间格记录 V、I、Q、Φ、空间 ROI、探针、能量、累计 R/G 耗散、泵浦做功、源功和端口流出功；历史数组可用 single，内部推进保持 double。
7. 电流 I 是论文磁场 Hz 的首选电路代理；电压 V 是电场代理。代理关系必须标注，不能声称等于 PCB 全波场。

### 7.4 FFT 能带

主路径由 run_tl_fdtd_fft 独立创建专用 bandCfg：

~~~matlab
bandField = tl_fdtd1d(model,bandCfg);
referenceField = tl_fdtd1d(referenceModel,referenceCfg);
observation = tl_make_observation( ...
    bandField.branch.IAtNodeTime,bandField.branch.x, ...
    bandField.node.t,bandCfg.cellMap);
referenceObservation = tl_make_observation( ...
    referenceField.branch.IAtNodeTime,referenceField.branch.x, ...
    referenceField.node.t,referenceCfg.cellMap);
fftBands = tl_xt_fft_bands( ...
    observation,referenceObservation,fftCfg);
~~~

tl_make_observation 可实现为入口局部函数，不必新增公开文件。referenceModel/referenceCfg 由入口显式定义为同源、同端口的无调制或无周期负载参考，具体选择写入结果；不得由 FFT 内核猜测。

要求：

1. 默认对有符号复 I 或 V 先做变换，再取功率。直接对 |I|² 变换会引入 DC 和二倍频；论文图注 FFT(|Hz|²) 的处理顺序不明确，因此只作为可选对照，不能作为默认科学路径。
2. 时间门使用无重复端点并覆盖整数个调制周期；空间 ROI 覆盖整数个单胞。分别使用空间窗和周期型时间 Hann 窗，并返回 coherent gain、energy gain、ENBW 和原生 bin 间隔。
3. 抗混叠检查先于折叠：recordEvery·dt 的 Nyquist 频率必须覆盖全部保留载频/Floquet 边带，单胞采样必须覆盖目标 k，时间窗必须与 T 公度。Floquet 折叠不能挽救采样阶段已经发生的混叠。
4. 先保留 signed k 与物理频率的 raw 谱，再把等价 bin 的功率求和到 ka∈[−π,π) 和 Re(ω)/Ω∈[−1/2,1/2)，不能只重排索引。
5. 零填充只细化显示网格，不计入真实分辨率。返回未归一化功率、全局归一化谱、逐 k 显示谱，以及由独立 referenceObservation 得到的源支持掩码；无参考支持的逐列归一化噪声不参与 ridge。
6. ridge 必须根据源参考谱、显著度和跨窗口连续性独立提取，不能用 PWE/TMM 目标值反选最近峰。
7. 主 x–t FFT 只验收 Re(ω)。bulk Im(ω) 优先由 PWE/TMM/Bloch 逐 k模式验证；若从有限链估计增长，先投影到固定 k 和单一 ridge/模态，并称为“有限链观测净增长”，因为它还受群速度流出、端口、拍频和多模混合影响。任何谱线宽都不解释为 Im(ω)。
8. 多节点单胞按 cellIndex 做 Bloch DFT，并使用 channelId/channelOffset 保留各通道谱；非相干功率和只作附加总览，不能覆盖通道分辨结果。

辅助 tl_bloch_fft_bands 对每个 k 推进单胞内全部状态并在 pT 同相位采样，复用根目录的奇偶频轴、Hann 窗和 Floquet 折叠约定。它用于无限体回归和有限链误差拆分，主有限链响应图仍使用 x–t FFT。

### 7.5 参数反演与实际元件选型

tl_fit_parameters.m 只用 Base MATLAB，可用对数参数化加 fminsearch 保证 L/C/R/G 为允许值。按可辨识性分阶段拟合：

1. 静态无损种子：fcol 固定 L0C0；布里渊区边界固定 Ls−2S；完整色散曲率分离 Ls 与 S。
2. 阻抗分离：仅有色散通常只能确定 LC 组合；必须加入单元输入阻抗、V/I 比或端口 S 参数才能分别确定 L 与 C。
3. 损耗分离：用静态传播衰减、谐振线宽和时域衰减联合拟合 Rs、Gp、变容管 ESR 与电感 ESR；只凭论文定义的 Q 不反演电阻。
4. 动态标定：用 C–V 数据和 BF 传递函数先生成 Ceff(t)，再用边带幅相或时域测量校正 ΔC、相位和高次谐波。
5. 端口/寄生：最后拟合源/负载、焊盘、过孔、隔直电容和探针响应，避免端口效应污染体参数。

tl_select_components.m 输出 target、nearest available、tolerance、ESR/Q、SRF、封装寄生和来源。器件筛选至少检查：

- 全泵浦周期内 Ceff>0，变容管始终满足允许偏置与反向电压；
- 电感与电容的自谐振频率高于所有纳入模型的载频与 Floquet 边带；
- 标称值与容差回算后仍满足目标静态色散、阻抗和调制深度；
- 源、负载与偏置网络不会在工作带内形成未建模强谐振；
- 固定随机种子的容差/寄生 Monte Carlo 报告 5%、50%、95% 分位，不只报告一个最佳样品。

## 8. 分阶段实施与阶段出口

| 阶段 | 主要工作 | 阶段出口 |
|---|---|---|
| P0 基线冻结 | 记录提交、路径、MATLAB 版本；建立参数来源表 | 无代码歧义，论文值与缺失值分级完成 |
| P1 统一模型 | 推导 SSPP/CROW 状态矩阵；实现 tl_build_model | 静态解析锚点和单位审计通过 |
| P2 TMM | 实现方波精确层与正弦/实测波形切片 | 常参数解析解、界面连续和切片收敛通过 |
| P3 PWE | Fourier 差阶、block-Toeplitz、本征筛选 | 与收敛 TMM 的复频交叉通过 |
| P4 FDTD | 无损链→损耗→互感→CROW→真实端口 | 能量/耗散闭合、dt/链长/边界收敛通过 |
| P5 FFT | 专用宽带 FDTD、x–t FFT、Bloch 辅助模式 | ridge 与 PWE/TMM 在分辨率内一致 |
| P6 参数反演 | 导入实测 C–V、S 参数、静态色散与衰减 | 输出带置信范围的 R/L/C/G/S |
| P7 器件落地 | 标称值选择、寄生、容差和泵浦约束 | 形成可采购候选表与回算报告 |
| P8 文档验收 | 四入口、README、validate_tl_suite、Code Analyzer | R2026a 全通过；R2020a 状态如实记录 |

依赖顺序为 P0→P1→P2→P3→P4→P5→P6→P7→P8。没有静态实测数据时，P1–P5 只能完成合成解析基准、37.18 nH/314.7 MHz 两个可识别论文锚点，以及显式 assumption 模型的求解器验证；不能称为完整论文参数验证，更不能宣称已经确定实际电路 R/L/C 值。

## 9. 预注册验收标准

### 9.1 路径、依赖与数据

- 从 Transmission line 目录运行四个入口；which 解析到本目录 tl_* 文件；
- MATLAB R2020a+、仅 Base MATLAB；本机 R2026a 实跑，未实测版本不得写成已通过；
- 无 addpath/startup、无根目录内核运行时依赖、无隐藏缓存；
- 非法零/负/NaN/Inf 元件、非正定互感矩阵、Hz/rad·s⁻¹ 混用、未对齐时间界面一律报错；
- 每个物理参数具有来源位置、单位、固定 status 和不确定度；论文未报告的不确定度明确写 unknown/not-reported；
- 从本目录执行 matlab -batch "validate_tl_suite"；任一断言失败必须使 batch 非零退出；
- 配置与结果保存不可变快照；生成物只进入 output/，外部 PDF 和原始测量数据默认不入库，只登记路径与校验和。

### 9.2 解析与求解器交叉

- 每个 benchmark 预先登记正的 ωscale、Qscale、Φscale、可观测权重阈值、简并点处理和比较模式集合，禁止用待比较理论曲线事后筛选模式；
- 常参数无损连续线恢复 vp=1/√(L′C′)、Zc=√(L′/C′)；无互感 SSPP 恢复离散 LC 解析色散，理想 CROW 恢复中心频率公式；
- 频率误差统一用 |Δω|/ωscale，避免 k=0/零模的相对误差零分母；TMM 常参数误差不高于 10⁻¹⁰，PWE 常系数误差不高于 10⁻⁸；
- 时间界面分别计算

\[
e_Q=\frac{\|Q^+-Q^-\|_2}
{\max(\|Q^-\|_2,Q_{\mathrm{floor}})},\qquad
e_\Phi=\frac{\|\Phi^+-\Phi^-\|_2}
{\max(\|\Phi^-\|_2,\Phi_{\mathrm{floor}})},
\]

  其中 floor 由预注册全局尺度给出；两者均不高于 10⁻¹⁰；
- 被动、无泵浦、R/G≥0 时要求 max Im(ω)/ωscale≤10⁻¹⁰；
- PWE–TMM 比较先按中心谐波/V/I 权重形成模式集合，再用双正交模态或 invariant-subspace 做一一匹配；简并点比较子空间/集合距离，不强制分支标签；
- Re(ω) 使用 minℓ|ΔRe(ω)+ℓΩ|/Ω 的循环距离，Im(ω) 使用 |ΔIm(ω)|/Ω 的普通绝对距离。截断/切片收敛后，两者分别满足中位数≤0.01、最大值≤0.05。

### 9.3 FDTD、端口与能量

- dt、dt/2、dt/4 使用同一物理网络、源、ROI 和共同记录时刻。复探针及场图在共同 ROI 上用 \(\|u-u_{\rm ref}\|_2/\max(\|u_{\rm ref}\|_2,u_{\rm floor})\) 比较，误差≤1%，且第二次加密误差至少下降；活跃频点的绝对相位误差≤0.01 rad，不使用“相位相对误差”；
- 封闭 periodic/open/short、无源、无调制、无损基准要求 \(|\mathcal E(t)-\mathcal E(0)|/\max(\mathcal E(0),\mathcal E_{\rm floor})\le0.1\%\)；
- 一般算例统一核对

\[
\Delta\mathcal E=W_{\mathrm{source}}+W_{\mathrm{pump}}
-D_{R/G}-W_{\mathrm{ports}},
\]

  残差除以各能量项绝对值之和与 \(\mathcal E_{\rm floor}\) 的较大者后不高于 1%。matched 端口的流出能量必须计入；Wpump 由保守状态和参数变化/界面功计算，用于捕获遗漏的 \(\dot C V\) 项；
- 合成均匀线在预注册通带、入射参考高于噪声阈值时，以端口方向波分解定义 Γ=reflected/incident，并要求 20log10|Γ|<−40 dB；实际周期线只报告频率依赖反射，不把全带 −40 dB 当作要求；
- 改变链长、记录 ROI、源位置和端口距离后，体模观测量在预注册的窗口不确定度内稳定。

### 9.4 FFT 与理论

- raw signed-k 谱、折叠谱、源支持掩码和未归一化功率全部保留；
- source-support mask 必须由同源同端口的独立无调制/无周期负载参考产生；参考低于预注册阈值的 bin 不归一化、不提 ridge；
- 设记录窗长度 Tobs、空间 ROI 长度 LROI，使用实际窗口返回的 ENBWt 与 ENBWx：

\[
\Delta\omega_{\mathrm{eff}}
=\mathrm{ENBW}_t\frac{2\pi}{T_{\mathrm{obs}}}
+|v_g|\,\mathrm{ENBW}_x\frac{2\pi}{L_{\mathrm{ROI}}}.
\]

  其中 vg 使用预注册的保守群速度上界。Bloch 逐 k 模式不含第二项。Bloch ridge 循环误差不超过时间有效分辨率；有限链 ridge 不超过 Δωeff，最大不超过 2Δωeff；
- 至少比较两种链长、两种记录时间、两种 ROI 和两个源位置；两次配置在公共支持区的 ridge 差不得超过两者 Δωeff 之和，并据此区分端口/Fabry–Pérot 峰与体带；
- zero padding 不参与上述门槛；
- FFT 线宽不作为 Im(ω) 验收。只有 Bloch 固定 k/单模投影且 |γref|Tfit≥2 时比较增长率，门槛为 max(0.05|γref|,1/Tfit)；近零增长只报告绝对误差。有限链增长始终标为观测净增长；
- 声明扫描区内若发生浮点溢出，该配置验收失败。必须缩短窗口、采用稳定缩放/QR 或明确收窄并重新验证范围；不得删除溢出样本后仍宣布整区通过。

### 9.5 实际参数与器件

- 论文 SSPP 种子回算 Ls−2S≈37.18 nH；
- 论文 CROW 种子回算 fcol≈314.7 MHz；
- 拟合使用加权训练数据并保留独立验证数据；静态色散验证 RMS 门槛预注册为 max(1%, 实验有效分辨率)，同时报告残差结构，禁止靠增加不可辨识寄生在训练集硬凑 1%；
- 最终参数表同时列 target、fit confidence、selected part、tolerance、parasitic 和数据来源；
- 容差扫描固定随机种子，至少 1000 个样本并检查分位数随样本数稳定，报告带边、净增长、特性阻抗和最大器件电压/电流的 5%、50%、95% 分位区间；
- 若缺少 C–V、S 参数或损耗测量，结果只能标为 simulation seed，不能标为 experimental component value。

## 10. 需要用户/实验侧提供的数据

在进入 P6 前，建议按优先级收集：

1. PCB 单胞完整尺寸、层叠、材料、铜厚、过孔和准确晶格常数 a；
2. 变容二极管型号及工作频率下的 C(V)、ESR(V)、封装寄生；
3. 电感、隔直电容和滤波网络的阻抗/S 参数、Q 与 SRF；
4. 未调制有限链或单胞的复 S11/S21，最好包含多个长度；
5. 静态逐单元复电压或电流场，用于提取 signed k 和衰减；
6. 泵浦在每个单元的幅值、相位、上升沿和空间同步误差；
7. 源、负载、连接器、探针和仪器的参考面/阻抗；
8. 动态调制下的原始复时空数据，而不仅是与论文图类似的强度图。

若只能先得到部分数据，拟合顺序保持不变，并冻结不可辨识参数；禁止让优化器用一个未测寄生去补偿另一个未测寄生。

## 11. 主要风险与控制

| 风险 | 后果 | 控制 |
|---|---|---|
| 把 ε/μ 简单替换为 C/L | 漏掉 R/G、寄生和 \dot C V | 以 q/φ 重推全部状态矩阵 |
| 把 w 当成 a | k 轴与电感反演错误 | a 保持必填待测字段 |
| 把论文 Q 当损耗品质因数 | 得到虚假 R | R 只由衰减/线宽/阻抗标定 |
| 只拟合色散 | L/C、R/G 不可辨识 | 联合阻抗、衰减和 V/I 数据 |
| 对电流强度 I² 直接 FFT | 产生 DC/二倍频伪带 | 默认先变换复 I/V，再取功率 |
| 忽略有限 Cblock | CROW 多出未建模极点 | 使用辅助状态直接建模 |
| 互感符号或矩阵非正定 | 非物理伪增长 | 正定审计与静态解析回归 |
| 把真实晶格色散当数值误差 | 错误“网格收敛” | 分开物理单胞与数值子划分 |
| 端口/有限长度峰当体带 | FFT ridge 偏离 | 链长、ROI、端口、Bloch 基准收敛 |
| 强泵浦数值溢出 | FFT 伪影或丢失衰减支 | 缩短窗、记录缩放、用 TMM/PWE 求复频 |
| 变容管非线性/泵浦耗尽 | 线性模型失效 | 先限定小信号；另立非线性阶段 |
| 论文缺少 Methods | 无法唯一复现 BOM | 所有缺失值显式标定，不作论文复现声明 |

## 12. 完成定义

以下条件全部满足，才认为“从麦克斯韦模板重构为实验传输线仿真”完成：

1. 四个入口均从 Transmission line 目录独立运行；所有物理电路参数只来自 tl_build_model，扫描、离散、源、ROI、窗口和绘图参数由各入口显式建立 solver config；
2. SSPP 与 CROW 的理想、损耗、互感和有限隔直电容模型都有明确状态方程；
3. PWE、TMM、Bloch FFT、有限链 FDTD/x–t FFT 的交叉验收通过；
4. FDTD 输出电路 V/I，并通过探针/全波校准映射到实验可测量；I 只作为 Hz 的电路代理，同时记录真实端口、源和调制元数据；
5. 实测 C–V、S 参数和衰减数据进入分级反演，输出带不确定度的 R/L/C/G/S；
6. 候选实际器件通过 SRF、Q/ESR、偏置、泵浦幅值、容差和寄生回算；
7. validate_tl_suite 非零失败时 MATLAB batch 返回失败，README 如实记录适用范围和未验证项；
8. 不把论文研究结论当验收目标，不把 simulation seed 误称为实际电路值。

建议正式编码从 P1 的 tl_build_model 与静态 SSPP/CROW 解析回归开始，再依次完成 TMM、PWE、FDTD 和 FFT；不要从绘图或论文动态图反向拼装参数。
