# Transmission line：时变传输线电路仿真

当前模型与根仓库同步为 **V3.1 / 3.1.0**。

本目录把根目录的麦克斯韦求解流程重构为离散传输线/LC 网络流程：PWE 与时间 Floquet TMM 求复准频率，有限链 FDTD 求电压/电流场，逐 `k_c` 有限宽度高斯波包与固定电压探针重建实准频率响应。内部单位统一为 SI，相位约定为

\[
\exp(ikx-i\omega t),\qquad \lambda=\exp(-i\omega T),
\]

所以 `Im(omega)>0` 表示增长。

参考论文只用于 SSPP/CROW 实验拓扑和文中明确给出的参数种子，不用于复现其研究结论。论文没有给出晶格常数、完整电感/互感、损耗、端口、器件型号和测量链；因此默认示例必须显式设置 `allowAssumptions=true`，其结果只能称为 simulation seed，不能称为实验 BOM。

## 文件与入口

- `tl_build_model.m`：SSPP/CROW 的唯一物理参数源和来源快照。
- `tl_pwe_fourier.m`、`tl_pwe_bands.m`：时间 Fourier 系数与 block-Toeplitz PWE。
- `tl_tmm_bands.m`：固定 Bloch 波数的时间单周期矩阵。
- `tl_fdtd1d.m`：有限离散链的 Q/Φ leapfrog 推进。
- `tl_fdtd_gaussian_k_scan.m`：逐 `k_c` 构造有限宽度复高斯电压波包并调用完整有限链 FDTD。
- `tl_fdtd_fft_bands.m`：多固定 V 探针时间 FFT、物理频率保留和显式 Floquet 功率折叠。
- `tl_bloch_fft_bands.m`：逐 k、同调制相位采样的无限周期体辅助基准，不作为实验式主图。
- `tl_xt_fft_bands.m`：一次端口宽带激励的有限链 signed x–t FFT 诊断工具；横轴语义不同于高斯源中心扫描。
- `tl_fit_parameters.m`：基于静态色散、阻抗、衰减和锚点的受约束反演接口。
- `tl_select_components.m`：用户器件目录筛选、SRF 检查和容差回算。
- `run_tl_pwe.m`、`run_tl_tmm.m`、`run_tl_fdtd_field.m`、`run_tl_fdtd_fft.m`：四个互相独立的示例入口。
- `validate_tl_suite.m`：本机回归与交叉验证。

在本目录运行，例如：

```matlab
run_tl_pwe
run_tl_tmm
run_tl_fdtd_field
run_tl_fdtd_fft
```

批量验证：

```text
/Applications/MATLAB_R2026a.app/bin/matlab -batch "validate_tl_suite"
```

核心函数只使用 Base MATLAB，不依赖根目录函数、工作区缓存或 MAT 文件。

## 统一模型

最小建模方式是：

```matlab
physicalCfg = struct('topology','sspp','allowAssumptions',true);
model = tl_build_model(physicalCfg);
```

`topology` 可为 `sspp` 或 `crow`。常用覆盖字段为 `a`、`Ls`、`mutualS`、`Rs`、`C0`、`Cpar`、`deltaC`、`Gp`、`fmHz`、`Zsource`、`Zload`、`cellCount`；CROW 另有 `L0`、`R0` 和 `Cblock`。总节点电容定义为 `C(t)=C0+Cpar+deltaC*f(t)`：`C0` 是名义调制支路，`Cpar` 是独立提取的加性寄生。论文所报 `C0` 被视为有效总值，因此论文种子取 `Cpar=0`，另设非零 `Cpar` 时必须避免重复计入寄生。每单元泵浦误差通过 `phaseByCell` 和 `amplitudeScaleByCell` 指定。

内部状态不是直接的 V/I，而是整数时间节点的电荷 Q 和半整数时间支路的磁链 Φ。这样在理想电容时间界面，Q/Φ 连续而 V/I 可以跃变，不会遗漏时变电容的泵浦项。`model.snapshot` 不含函数句柄，可直接随结果归档。

论文种子的两个首要单位回归是：SSPP 的 `Ls-2*S≈37.18 nH`，以及在默认 `Cpar=0` 时 CROW 的 `1/(2*pi*sqrt(L0*C0))≈314.7 MHz`。

## PWE 与 TMM

```matlab
k = linspace(-pi/model.cell.a,pi/model.cell.a,181);
pweCfg = struct('Mtime',5,'Nt',256);
fourier = tl_pwe_fourier(model,pweCfg);
pwe = tl_pwe_bands(fourier,model,k,pweCfg);

tmmCfg = struct('temporalSlices',1024);
tmm = tl_tmm_bands(model,k,tmmCfg);
```

PWE 返回原始 Floquet 副本、折叠谱、中心谐波/可观测权重、选中分支和左右本征矢条件数。TMM 对方波使用精确常值层，对正弦波使用中点 `expm` 切片。两者都不在简并点强行赋予物理分支身份；提高 `Mtime/Nt` 或 `temporalSlices` 后应做集合意义的收敛比较。

## 有限链 FDTD

```matlab
fdtdCfg = struct();
fdtdCfg.dt = model.modulation.period/256;
fdtdCfg.nSteps = 8*256;
fdtdCfg.boundaryType = 'matched';  % open/matched/periodic/short
fdtdCfg.modulationEnabled = true;
fdtdCfg.source = struct('type','thevenin','node',1, ...
    'impedance',50,'waveformFcn',@(t) exp(-((t-1e-9)/0.2e-9).^2));
field = tl_fdtd1d(model,fdtdCfg);
```

`IHalf0/PhiHalf0` 定义在 `t=-dt/2`。结果显式分开：

```text
field.node.x, field.node.t, field.node.V, field.node.Q
field.branch.x, field.branch.tHalf
field.branch.IHalf, field.branch.PhiHalf, field.branch.IAtNodeTime
```

`IAtNodeTime` 是相邻两个半时间电流的中心平均，允许与节点时间数据共同绘图或做 x–t FFT。内部推进始终为 double，`precision='single'` 只压缩历史数组。稳定性强制 `dt*omegaMaximum<=1.8`。逐 `k_c` 扫描只复用带完整模型快照和边界拓扑的 `stabilityAudit`，不会省略任何有限链推进。`open` 表示电路开路端，并非无反射空间边界；`matched` 只是按给定端口阻抗加终端电导。`field.energy` 给出源功、泵浦功、R/G 耗散、端口耗散和数值账本残差。

`I` 只是局域磁场的电路代理，`V` 只是电场代理；两者不等于 PCB 的三维全波场。

## FFT 输入约束

主入口 `run_tl_fdtd_fft` 不再把一次端口脉冲的空间 FFT 波数或无限周期体的 Bloch `k` 冒充实验式逐列激发。它对每个高斯源中心 `k_c` 独立构造具有指定强度 FWHM 的复电压波包，用波包中心处的瞬时局域单胞模态补全 Q/Φ（CROW 还补全谐振器状态），调用 `tl_fdtd1d` 推进完整有限链，只把多个固定节点的复电压 `V(t)` 交给 FFT。若 `phaseByCell/amplitudeScaleByCell` 使初始电容非均匀，代码仍用每节点实际电容令 `Q_i=C_iV_i`，保证指定的 `V(x,0)` 精确成立，同时警告其余伴随状态只是中心局域窄带近似。无损模型的群速度在邻近 k 用电路能量内积连续跟踪同一模态；有损模型会明确标记这只是右本征矢连续性启发式，例外点附近应改用端口激励或双正交分析。半时间电流包络再按离散 leapfrog 群速度回退半步；SSPP 的 `k_c=0` 模式用 `zeroKPropagationDirection` 指定的 `k→0` 阻抗极限并核查定向功率，若有损/过阻尼模型没有可分辨传播方向则拒绝伪造“单向”初值：

```matlab
scan = tl_fdtd_gaussian_k_scan(model,scanCfg);
bands = tl_fdtd_fft_bands( ...
    scan.probeSignals,scan.time,scan.k,fftCfg);
```

时间窗强制采用 `[start,end)` 的无重复端点、至少两个整数调制周期和周期型 Hann 窗。程序先保留完整 Nyquist 区的 `rawOmega/rawPower`，再在未移位整数 DFT bin 上把相差整数倍 Ω 的功率显式累加到 `[-Omega/2,Omega/2)`；奇偶采样数都不依赖 reshape 假设。多个固定 V 探针的功率非相干相加，降低单点位于场节点导致的漏支。

结果同时返回 `foldedPower`、`rawPower`、未归一化 `rawColumnPower`、`activeKMask`、全局/逐列归一化功率和窗/采样元数据。所有 FFT 功率使用周期 Hann 的相干振幅校正；格点复指数的峰值不随零填充因子改变，连续谱积分则必须乘相应的频率 bin 宽（x–t 谱乘 k–ω bin 面积）。逐列归一化只用于观察峰位；零填充只改变绘图间距，原生归一化频率分辨率仍是 `1/analysisPeriodCount`。横轴是有限波包的中心 `k_c`，每列含宽度约由 FWHM 决定的 Bloch 波数组合，因此只能称为有限源、有限链、有限窗响应，不能当成无限体精确本征值，也不能从线宽读取 `Im(omega)`。

`tl_xt_fft_bands` 仍可分析一次宽带端口激励的 signed x–t 数据；其可选 `maximumPhysicalOmega` 只限制进入 Floquet 折叠的 raw 频率 bin，完整 raw 谱仍保留。`tl_bloch_fft_bands` 仍可作无限周期体数值回归。两者的横轴与边界条件必须分别标注，不能与上述主图混用。科研使用应分别改变波包 FWHM、探针位置、链长/边界距离、`dt`、分析周期数和零填充因子，确认峰位对前五项收敛且不把零填充误作分辨率提升。

## 参数反演

```matlab
fitData.staticDispersion = struct( ...
    'kRadPerM',kMeasured, ...
    'frequencyHz',fMeasured, ...
    'uncertaintyHz',sigmaF);
fitCfg.parameterNames = {'Ls','mutualS'};
fit = tl_fit_parameters(physicalCfg,fitData,fitCfg);
```

可用数据块和字段见 `help tl_fit_parameters`。只有色散时通常只能稳健确定 LC 组合；要分别确定 L/C 需要 V/I 比、输入阻抗或 S 参数。要分离 `Rs/Gp/R0`，需要传播衰减、谐振线宽和时域衰减。`isLocallyIdentifiable=false` 时，不应引用单个最佳值作为实验参数。最终还应保留独立验证集。

## 器件目录与筛选

```matlab
catalog.inductors = struct( ...
    'partNumber','example-L','value',39e-9, ...
    'toleranceFraction',0.05,'esrOhm',0.4,'srfHz',4e9, ...
    'package','0402','source','manufacturer datasheet');
catalog.capacitors = struct( ...
    'partNumber','example-C','value',20e-12, ...
    'toleranceFraction',0.05,'esrOhm',0.1,'srfHz',6e9, ...
    'package','0402','source','manufacturer datasheet');
report = tl_select_components(model,catalog,struct());
```

`value` 对电感用 H，对电容/变容管用 F；`toleranceFraction=0.05` 表示 ±5%。筛选只把满足 `srfHz >= srfSafetyFactor*maximumFrequencyHz` 的候选纳入最近标称值比较，并至少运行 1000 个固定种子容差样本。`C0` 与 `Cpar` 都是去嵌后的等效量，不能直接用一只最近值电容替代；如果输入的 `C0` 已包含焊盘、偏置和封装寄生，就必须保持 `Cpar=0`。

## 进入实验参数阶段前仍需的数据

需要补充准确晶格常数与 PCB 层叠、变容管 C(V)/ESR/封装寄生、L/C 的 RF 阻抗与 SRF、未调制样品复 S11/S21、多长度衰减、逐单元复场，以及每单元泵浦幅相/上升沿。缺少这些数据时，本目录可以验证求解器和合成电路，但不能唯一确定实际 R/L/C/G/S。
