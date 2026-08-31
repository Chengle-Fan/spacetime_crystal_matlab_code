# Spacetime Media MATLAB — V3 / 时空介质 MATLAB V3 版

一维光学时间晶体的数值仿真与理论计算工具：平面波展开法 (PWE)、时间传输矩阵法 (TMM)、有限样品时域有限差分法 (FDTD)，以及逐 k、有限时间原胞同相位采样的 FFT 能带重建。四个入口脚本扁平化位于仓库根目录，无需任何路径设置。

A MATLAB toolkit for 1D scalar photonic time crystals: PWE/TMM bands, finite-sample FDTD field evolution, and per-k finite stroboscopic FDTD–FFT quasifrequency reconstruction. Four flat entry scripts are provided at the repository root.

---

## 四个执行入口 / The Four Entry Scripts

| 脚本 | 方法 | 内容 |
|------|------|------|
| `run_pwe.m` | 平面波展开 (PWE) | 可选方波/正弦调制的空间均匀时间晶体；只绘制复准频率的实部与虚部 |
| `run_tmm.m` | 时间传输矩阵法 (TMM) | 两层方波时间晶体；只绘制复准频率的实部与虚部 |
| `run_fdtd_field.m` | 时域有限差分 (FDTD) | 以用户指定中心 k 或波长的高斯初值激发，只绘制所选 E/D 空间—时间场分布 |
| `run_fdtd_fft.m` | 逐 k FDTD–FFT | 对每个 k 推进有限个时间原胞并同相位采样，只绘制第一时间 Floquet 区能带 |

在仓库根目录直接运行：

```matlab
run_pwe      % PWE：设置时间晶体参数，绘制 Re(ω) 与 Im(ω)
run_tmm      % TMM：设置方波时间晶体参数，绘制 Re(ω) 与 Im(ω)
run_fdtd_field  % 有限空间样品：只绘制用户选择的 |E| 或 |D|
run_fdtd_fft    % 独立逐 k 有限采样：只绘制 FFT Floquet 能带
```

`run_pwe` 和 `run_tmm` 各生成一张包含实部/虚部的图。两个 FDTD 入口彼此独立：`run_fdtd_field` 研究一个用户所选中心 k 的有限样品场演化；`run_fdtd_fft` 不读取前者数据，而是遍历 k，对每个 k 只记录有限个时间原胞的同相位电场样本，再沿原胞序列做 FFT。两个脚本都只生成一张图。

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
| PWE | 空间均匀、时间周期、无损、严格正 ε、μ=1 的普通介质 | 两条复准频率能带的实部与虚部 |
| TMM | 空间均匀、无损、严格正 ε/μ 的两层方波时间晶体 | 两条精确复准频率能带的实部与虚部 |
| FDTD | 可选有限空间/全空间区域与有限开启时窗；方波或正弦调制 | 初值激发后的复电场 E(x,t) 或电位移 D(x,t) |
| FFT | 每个守恒实 k 独立推进；有限个时间原胞在相同调制相位采样 E | 第一时间 Floquet 区内的有限窗实准频率重建谱 |

限制：标量一维、非色散、无损且 ε/μ 严格为正。PWE 支持方波/正弦调制；TMM 只支持两层方波。有限样品 FDTD 内部始终推进时间界面连续的 D/B，激发只由复数 E/H 初值给出。默认 `fixed-distant` 是经过传播距离检查的固定远端截断，并非数学开放边界；可选海绵也不是 PML，且在 k-gap 增益问题中可能散射出被放大的额外 k 分量。`sampleCellCount` 与 `cellSize` 只定义有限样品的几何长度，不代表空间晶格。

`run_fdtd_field` 可设置有限的 `modulationStartTime/modulationEndTime`、方波开启相位、`sample/global` 调制区域、`normalized/um-fs` 单位，以及 `normalized-k/wavelength` 两种激发输入。脉冲宽度明确采用强度 FWHM。`fdtd1d` 可只记录 E 或 D、只保存指定空间 ROI，并允许历史数组采用 single；内部推进仍保持 double。对解析上已知材料下界的方波/正弦入口，可使用认证下界避免重复的全时空 CFL 预扫描，实际推进时仍检查材料未违反该下界。

FFT 入口同样使用有限数据：每个 k 只有 `temporalCellCount` 个同相位采样点，观测跨度为 `(temporalCellCount-1)*T`，原生分辨率约为 `Delta[Re(omega)/Omega]=1/temporalCellCount`。周期型 Hann 窗会展宽谱峰，零填充只细化绘图网格。等相位采样把相差整数倍 Ω 的 Floquet 副本自动映射到 `[-Omega/2,Omega/2)`；不进行空间 FFT、空间折叠或有限样品 ROI 处理。主图仍按每个 k 独立归一化，只表示实准频率峰位；函数另返回未归一化列功率、后半窗主导增长率和 Courant/FFT/窗元数据。该增长率只近似 k-gap 的正增长支，不能由谱宽恢复完整的复频率两支。空间推进采用 `kYee=2*sin(k*dx/2)/dx`，科研使用时必须检查 dx、dt 与时间原胞数收敛。详见各入口脚本的逐段中文注释与 [repair.md](repair.md)。

论文复现见 [`topology_aspect/`](topology_aspect/README.md)。该目录不复制模板求解器；复现入口会验证并直接调用本根目录的 TMM/FDTD 函数，同时在返回值中记录实际解析路径。

## 版本与旧版恢复 / Version & Legacy

- 当前版本 **3.0.1**（见 `VERSION.txt`）。当前用户收窄后的 PWE、TMM、FDTD 专项实现记录见 `repair.md` 第 8–10 节；它们优先于该文件前半部分保留的历史验收要求。
- V3 曾移除旧版 `reproduction/` 和其他历史目录；当前 `topology_aspect/` 是随后按指定论文重新建立的独立复现目录，不是旧版目录恢复。
- 旧版全部代码可从 git 标签 **`legacy-pre-v3-refactor`** 恢复。
- 最小环境：**MATLAB R2020a+**（Base MATLAB，无 Toolbox 依赖）。四个入口脚本均从仓库根目录运行，无 `addpath`/startup 依赖；本次实测于 **R2026a**，R2020a 为声明下限（本机未安装，未实测）。

## 参考文献 / References

- J. Park and B. Min, "Spatiotemporal plane wave expansion method for arbitrary space-time periodic photonic media," *Optics Letters* **46**, 484–487 (2021). [DOI: 10.1364/OL.411622](https://doi.org/10.1364/OL.411622)
- D. Ramaccia, A. Alù, A. Toscano, and F. Bilotti, "Temporal multilayer structures for designing higher-order transfer functions using time-varying metamaterials," *Applied Physics Letters* **118**, 101901 (2021). [DOI: 10.1063/5.0042567](https://doi.org/10.1063/5.0042567)
