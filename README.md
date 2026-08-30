# Spacetime Media MATLAB — V3 / 时空介质 MATLAB V3 版

一维光学时间晶体的数值仿真与理论计算工具：平面波展开法 (PWE)、时间传输矩阵法 (TMM)、时域有限差分法 (FDTD) 和场数据二维 FFT。四个入口脚本扁平化位于仓库根目录，无需任何路径设置。

A MATLAB toolkit for 1D scalar photonic time crystals: plane-wave expansion (PWE), temporal transfer-matrix method (TMM), finite-difference time-domain (FDTD), and field-history FFT reconstruction. Four flat entry scripts are provided at the repository root.

---

## 四个执行入口 / The Four Entry Scripts

| 脚本 | 方法 | 内容 |
|------|------|------|
| `run_pwe.m` | 平面波展开 (PWE) | 可选方波/正弦调制的空间均匀时间晶体；只绘制复准频率的实部与虚部 |
| `run_tmm.m` | 时间传输矩阵法 (TMM) | 两层方波时间晶体；只绘制复准频率的实部与虚部 |
| `run_fdtd_field.m` | 时域有限差分 (FDTD) | 初值激发有限时间晶体，只绘制空间—时间场分布并保存复场数据 |
| `run_fdtd_fft.m` | 二维 FFT | 读取有限样品复场，只绘制第一时间 Floquet 区内的局部 FFT 能带 |

在仓库根目录直接运行：

```matlab
run_pwe      % PWE：设置时间晶体参数，绘制 Re(ω) 与 Im(ω)
run_tmm      % TMM：设置方波时间晶体参数，绘制 Re(ω) 与 Im(ω)
run_fdtd_field  % 第一步：计算有限样品复场，只绘制 |E(x,t)|
run_fdtd_fft    % 第二步：读取上一步数据，只绘制 FFT 局部能带
```

`run_pwe` 和 `run_tmm` 各生成一张包含实部/虚部的图。`run_fdtd_field` 生成一张场分布图，并把第二步必需的完整复电场保存到已被 Git 忽略的 `fdtd_field_data.mat`；`run_fdtd_fft` 独立读取该文件并只生成一张 FFT 能带图。

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
| FDTD | 静态背景中的有限、空间均匀时间调制样品；方波或正弦调制 | 初值激发后的复电场 E(x,t) |
| FFT | 有限样品内部、整数个调制周期的复电场记录 | 初值加权的局部实频响应谱 |

限制：标量一维、非色散、无损且 ε/μ 严格为正。PWE 支持方波/正弦调制；TMM 只支持两层方波。FDTD 内部始终推进时间界面连续的 D/B，外边界海绵不是 PML，激发只由复数 E/H 初值给出。`sampleCellCount` 与 `cellSize` 只定义有限样品的几何长度，不代表空间晶格。`excitationRegion='bulk'` 或 `'k-gap'` 选择初始波包中心 k；单次 Gaussian 初值只能重建该 k 附近的局部谱。FFT 使用样品内部 ROI、空间/时间 Hann 窗、全谱统一 dB 归一化和有符号 k，不做空间折叠、不逐 k 放大弱列，也不能恢复 k-gap 内的 Im(ω)。详见各入口脚本头部注释与 [repair.md](repair.md)。

## 版本与旧版恢复 / Version & Legacy

- 当前版本 **3.0.1**（见 `VERSION.txt`）。当前用户收窄后的 PWE、TMM、FDTD 专项实现记录见 `repair.md` 第 8–10 节；它们优先于该文件前半部分保留的历史验收要求。
- V3 为破坏性重构：旧目录（`core/`、`tmm/`、`fdtd/`、`topology/`、`demos/`、`examples/`、`tests/`、`learn/`、`docs/`、`reproduction/`、`Simple_FDTD_NathanZechar/`、`output/`）与旧根入口（`startup_stm.m`、`stm_init.m`、`run_all_demos.m`、`stm_run_examples.m`）均已移除。
- 旧版全部代码可从 git 标签 **`legacy-pre-v3-refactor`** 恢复。
- 最小环境：**MATLAB R2020a+**（Base MATLAB，无 Toolbox 依赖）。四个入口脚本均从仓库根目录运行，无 `addpath`/startup 依赖；本次实测于 **R2026a**，R2020a 为声明下限（本机未安装，未实测）。

## 参考文献 / References

- J. Park and B. Min, "Spatiotemporal plane wave expansion method for arbitrary space-time periodic photonic media," *Optics Letters* **46**, 484–487 (2021). [DOI: 10.1364/OL.411622](https://doi.org/10.1364/OL.411622)
- D. Ramaccia, A. Alù, A. Toscano, and F. Bilotti, "Temporal multilayer structures for designing higher-order transfer functions using time-varying metamaterials," *Applied Physics Letters* **118**, 101901 (2021). [DOI: 10.1063/5.0042567](https://doi.org/10.1063/5.0042567)
