# Spacetime Media MATLAB — V3.1 / 时空介质 MATLAB V3.1

一维光学时间晶体的数值仿真与理论计算工具：平面波展开法 (PWE)、时间传输矩阵法 (TMM)、有限样品时域有限差分法 (FDTD)，以及逐 k 高斯波包激发、固定探针采样的 FFT 能带重建。四个入口脚本扁平化位于仓库根目录，无需任何路径设置。

A MATLAB toolkit for 1D scalar photonic time crystals: PWE/TMM bands, finite-sample FDTD field evolution, and per-k Gaussian-wavepacket FDTD with fixed-probe FFT quasifrequency reconstruction. Four flat entry scripts are provided at the repository root.

---

## 四个执行入口 / The Four Entry Scripts

| 脚本 | 方法 | 内容 |
|------|------|------|
| `run_pwe.m` | 平面波展开 (PWE) | 可选方波/正弦调制的空间均匀时间晶体；实部显示完整 Floquet 谱的指定副本，虚部仅显示两条代表根 |
| `run_tmm.m` | 时间传输矩阵法 (TMM) | 两层方波时间晶体；只绘制复准频率的实部与虚部 |
| `run_fdtd_field.m` | 时域有限差分 (FDTD) | 以用户指定中心 k 或波长的高斯初值激发，只绘制所选 E/D 空间—时间场分布 |
| `run_fdtd_fft.m` | 逐 k FDTD–FFT | 对每个中心 k 构造高斯波包并运行有限空间 FDTD；保存代表性 E(x,t)，由固定 E 探针重建第一时间 Floquet 区能带 |

在仓库根目录直接运行：

```matlab
run_pwe      % PWE：设置时间晶体参数，绘制 Re(ω) 与 Im(ω)
run_tmm      % TMM：设置方波时间晶体参数，绘制 Re(ω) 与 Im(ω)
run_fdtd_field  % 有限空间样品：只绘制用户选择的 |E| 或 |D|
run_fdtd_fft    % 逐 k 高斯波包 FDTD：绘制代表性场图与探针 FFT 能带
```

`run_pwe` 和 `run_tmm` 各生成一张包含实部/虚部的图。两个 FDTD 入口彼此独立：`run_fdtd_field` 研究一个用户所选中心 k 的有限样品场演化；`run_fdtd_fft` 不读取前者数据，而是遍历中心 k，对每个 k 使用相同包络构造一组新的高斯 E/H 初值，调用同一个 `fdtd1d` 有限空间内核，并在样品内部固定位置记录 E(t)。FFT 入口默认生成一张代表性 E(x,t) 图和一张能带图。

## 统一归一化 / Normalization

- 真空光速、介电常数、磁导率均归一化：**c0 = eps0 = mu0 = 1**，默认 **mu_r = 1**。
- 时间周期 T，调制角频率 Ω = 2π/T。
- 长度与时间为无量纲单位；无量纲波数常用 k/Ω（即 k·T/(2π)）。

## 相位约定 / Phase Conventions

- 介电函数展开：**ε(t) = Σ_m ε_m · exp(−i·m·Ω·t)**。
- 准频率：λ = exp(−i·ω·T)，ω = i·log(λ)/T（主对数，Re(ω) 折进第一时间布里渊区）。
- 纯时间晶体保持空间波数 k 守恒且不折叠；只把 Re(ω) 折入 [−Ω/2, Ω/2)。

## 方法适用范围 / Method Scope

| 方法 | 适用范围 | 输出 |
|------|---------|------|
| PWE | 空间均匀、时间周期、无损、严格正 ε、μ=1 的普通介质 | 完整 Floquet 本征谱；`run_pwe` 用实部纵轴范围显示副本，虚部只绘制两条代表根 |
| TMM | 空间均匀、无损、严格正 ε/μ 的两层方波时间晶体 | 两条精确复准频率能带的实部与虚部 |
| FDTD | 可选有限空间/全空间区域与有限开启时窗；方波或正弦调制 | 初值激发后的复电场 E(x,t) 或电位移 D(x,t) |
| FFT | 每个中心 k 独立进行有限空间高斯波包 FDTD；在样品内多个固定点采样 E(t) | 第一时间 Floquet 区内的有限源、有限样品、有限窗实准频率响应谱 |

限制：标量一维、非色散、无损且 ε/μ 严格为正。PWE 支持方波/正弦调制；TMM 只支持两层方波。有限样品 FDTD 内部始终推进时间界面连续的 D/B，激发只由复数 E/H 初值给出。默认 `fixed-distant` 是经过传播距离检查的固定远端截断，并非数学开放边界；可选海绵也不是 PML，且在 k-gap 增益问题中可能散射出被放大的额外 k 分量。`sampleCellCount` 与 `cellSize` 只定义有限样品的几何长度，不代表空间晶格。

`run_fdtd_field` 可设置有限的 `modulationStartTime/modulationEndTime`、方波开启相位、`sample/global` 调制区域、`normalized/um-fs` 单位，以及 `normalized-k/wavelength` 两种激发输入。脉冲宽度明确采用强度 FWHM。`fdtd1d` 可只记录 E 或 D、只保存指定空间 ROI，并允许历史数组采用 single；内部推进仍保持 double。对解析上已知材料下界的方波/正弦入口，可使用认证下界避免重复的全时空 CFL 预扫描，实际推进时仍检查材料未违反该下界。

FFT 入口使用真实有限空间数据：`run_fdtd_fft` 通过共享函数 `fdtd_gaussian_k_scan` 对每个 `k_c` 构造具有指定强度 FWHM 的单向复高斯波包，调用 `fdtd1d` 推进完整空间 E/D/B/H 状态，只把样品内部多个固定点的 E(t) 送入 `fdtd_fft_bands`。V3.1 对 `k_c=0` 要求用 `zeroKPropagationDirection` 明确选择方向，磁场半步包络使用 Yee 离散群速度，不再因 `sign(0)=0` 退化成双向分裂。`topology_aspect/reproduce_fdtd_band.m` 也调用同一组根目录函数验证论文参数，不再复制逐 k 初值与推进代码。分析窗采用 `[start,end)` 的整数周期、无重复端点采样和周期型 Hann 窗；先计算完整物理频率谱，再把相差整数倍 Ω 的功率显式累加到 `[-Omega/2,Omega/2)`。多个探针的功率非相干相加，降低单点恰落在场节点造成的漏支风险。零填充只细化绘图网格，原生分辨率仍约为 `1/analysisPeriodCount`。返回值同时保留折叠功率、raw 物理频率功率、全局/逐列归一化功率、有效列掩码和窗/采样元数据；raw/折叠功率采用 Hann 相干增益校正，格点复指数的峰值不随零填充改变，做频率积分时仍须乘 bin 宽。主图逐 k 列归一化，只用于看峰位。横轴是高斯源的中心 `k_c`，每列还含有限波包谱宽，因此图是源加权的有限样品响应，不能冒充无限体精确本征值，也不能由谱宽读取 `Im(omega)`。科研使用必须检查 FWHM、探针位置、样品/边界距离、dx、dt 和分析周期数收敛。详见各入口脚本的逐段中文注释与 [repair.md](repair.md)。

传输线模板见 [`Transmission line/`](Transmission%20line/README.md)，对应论文图 2/3 的传输线复现见 [`Full_momemtum/`](Full_momemtum/README.md)，光学时间晶体复现见 [`topology_aspect/`](topology_aspect/README.md)。论文目录不复制求解器；运行时会验证共享模板的实际解析路径并随结果归档。

## 版本与旧版恢复 / Version & Legacy

- 当前版本 **V3.1 / 3.1.0**（见 `VERSION.txt`）。总核查与修正记录见 `repair.md` 第 15 节；它优先于第 1–14 节中的历史实现描述。
- V3 曾移除旧版 `reproduction/` 和其他历史目录；当前 `topology_aspect/` 是随后按指定论文重新建立的独立复现目录，不是旧版目录恢复。
- 旧版全部代码可从 git 标签 **`legacy-pre-v3-refactor`** 恢复。
- 一键核查：`summary = validate_v31_suite();`；加入两套论文有限样品冒烟测试可用 `validate_v31_suite(struct('runPaperSmoke',true))`。
- 最小环境：**MATLAB R2020a+**（Base MATLAB，无 Toolbox 依赖）。四个入口脚本均从仓库根目录运行，无 `addpath`/startup 依赖；本次实测于 **R2026a Update 3**，R2020a 为声明下限（本机未安装，未实测）。

## 参考文献 / References

- J. Park and B. Min, "Spatiotemporal plane wave expansion method for arbitrary space-time periodic photonic media," *Optics Letters* **46**, 484–487 (2021). [DOI: 10.1364/OL.411622](https://doi.org/10.1364/OL.411622)
- D. Ramaccia, A. Alù, A. Toscano, and F. Bilotti, "Temporal multilayer structures for designing higher-order transfer functions using time-varying metamaterials," *Applied Physics Letters* **118**, 101901 (2021). [DOI: 10.1063/5.0042567](https://doi.org/10.1063/5.0042567)
