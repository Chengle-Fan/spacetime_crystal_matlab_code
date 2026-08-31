# Transmission line：时变传输线电路仿真

本目录把根目录的麦克斯韦求解流程重构为离散传输线/LC 网络流程：PWE 与时间 Floquet TMM 求复准频率，有限链 FDTD 求电压/电流场，Bloch 或有限链 x–t FFT 重建实准频率能带。内部单位统一为 SI，相位约定为

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
- `tl_bloch_fft_bands.m`：逐 k、同调制相位采样的无限体 FFT 基准。
- `tl_xt_fft_bands.m`：有限链 signed x–t FFT、Floquet 折叠和参考支持 ridge。
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

`topology` 可为 `sspp` 或 `crow`。常用覆盖字段为 `a`、`Ls`、`mutualS`、`Rs`、`C0`、`deltaC`、`Gp`、`fmHz`、`Zsource`、`Zload`、`cellCount`；CROW 另有 `L0`、`R0` 和 `Cblock`。每单元泵浦误差通过 `phaseByCell` 和 `amplitudeScaleByCell` 指定。

内部状态不是直接的 V/I，而是整数时间节点的电荷 Q 和半整数时间支路的磁链 Φ。这样在理想电容时间界面，Q/Φ 连续而 V/I 可以跃变，不会遗漏时变电容的泵浦项。`model.snapshot` 不含函数句柄，可直接随结果归档。

论文种子的两个首要单位回归是：SSPP 的 `Ls-2*S≈37.18 nH`，以及 CROW 的 `1/(2*pi*sqrt(L0*C0))≈314.7 MHz`。

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

`IAtNodeTime` 是相邻两个半时间电流的中心平均，允许与节点时间数据共同绘图或做 x–t FFT。内部推进始终为 double，`precision='single'` 只压缩历史数组。稳定性强制 `dt*omegaMaximum<=1.8`。`field.energy` 给出源功、泵浦功、R/G 耗散、端口耗散和数值账本残差。

`I` 只是局域磁场的电路代理，`V` 只是电场代理；两者不等于 PCB 的三维全波场。

## FFT 输入约束

有限链 FFT 需要两次独立运行：被测调制链和同源、同端口的无调制（或明确移除周期负载的）参考链。输入结构必须含：

```text
data, x, t, cellIndex, channelId, channelOffset, observableName
```

`data` 为 `Nt×Nsample`。每个通道必须在同一组连续单胞中各采一次；`x=cellIndex*a+channelOffset+共同原点`。时间数据采用无重复终点，并满足 `Nt*dt` 是整数个调制周期。默认先变换有符号复 V/I，再取功率；不要默认对 `abs(I).^2` 做 FFT。

```matlab
fftCfg = struct('cellPeriod',model.cell.a, ...
    'temporalPeriod',model.modulation.period);
bands = tl_xt_fft_bands(observation,referenceObservation,fftCfg);
```

结果同时保留 `rawPower/rawOmega`、第一 Floquet 区 `foldedPower/foldedOmega`、分通道功率、独立参考功率、`sourceSupportMask`、全局/逐列显示归一化和不借助 PWE/TMM 的连续 ridge。零填充只改变显示间距；真实分辨率见 `nativeOmegaResolution` 和 `nativeKResolution`。有限链 ridge 只解释 `Re(omega)`，谱线宽不是 `Im(omega)`。

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

`value` 对电感用 H，对电容/变容管用 F；`toleranceFraction=0.05` 表示 ±5%。筛选只把满足 `srfHz >= srfSafetyFactor*maximumFrequencyHz` 的候选纳入最近标称值比较，并至少运行 1000 个固定种子容差样本。有效 `C0` 包含变容管、偏置、焊盘和其他寄生，不能直接用一只最近值电容替代。

## 进入实验参数阶段前仍需的数据

需要补充准确晶格常数与 PCB 层叠、变容管 C(V)/ESR/封装寄生、L/C 的 RF 阻抗与 SRF、未调制样品复 S11/S21、多长度衰减、逐单元复场，以及每单元泵浦幅相/上升沿。缺少这些数据时，本目录可以验证求解器和合成电路，但不能唯一确定实际 R/L/C/G/S。

