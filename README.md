# Spacetime Media MATLAB — V3 / 时空介质 MATLAB V3 版

一维光子时空晶体（photonic spacetime crystals）的数值仿真与理论计算工具：平面波展开法 (PWE)、时间传输矩阵法 (TMM)、时域有限差分法 (FDTD)。三个入口脚本扁平化位于仓库根目录，无需任何路径设置。

A MATLAB toolkit for 1D scalar photonic spacetime crystals: plane-wave expansion (PWE), temporal transfer-matrix method (TMM), and finite-difference time-domain (FDTD). Three flat entry scripts at the repo root — no `addpath`, no startup file.

---

## 三个执行入口 / The Three Entry Scripts

| 脚本 | 方法 | 内容 |
|------|------|------|
| `run_pwe.m` | 平面波展开 (PWE) | 可选方波/正弦调制的空间均匀时间晶体；只绘制复准频率的实部与虚部 |
| `run_tmm.m` | 时间传输矩阵法 (TMM) | 两层方波时间晶体；只绘制复准频率的实部与虚部 |
| `run_fdtd.m` | 时域有限差分 (FDTD) | 波包演化 + 有限样品 E-FFT 能带（对 TMM 基准，含 1/L 收敛与无样品参考）；一般时空晶体的双重折叠谱 |

在仓库根目录直接运行：

```matlab
run_pwe      % PWE：设置时间晶体参数，绘制 Re(ω) 与 Im(ω)
run_tmm      % TMM：设置方波时间晶体参数，绘制 Re(ω) 与 Im(ω)
run_fdtd     % FDTD：波包、宽带 FFT 能带、双重折叠谱
```

结果留在工作区；默认不写盘、不创建目录。`run_pwe` 和 `run_tmm` 各只生成一张包含实部/虚部的图。

## 统一归一化 / Normalization

- 真空光速、介电常数、磁导率均归一化：**c0 = eps0 = mu0 = 1**，默认 **mu_r = 1**。
- 空间周期 Λ、时间周期 T；g = 2π/Λ，Ω = 2π/T。
- 长度与时间为无量纲单位；无量纲波数常用 k/Ω（即 k·T/(2π)）。

## 相位约定 / Phase Conventions

- 介电函数展开：**ε(x,t) = Σ_{m,n} ε_mn · exp(i·n·g·x − i·m·Ω·t)**。
- 准频率：λ = exp(−i·ω·T)，ω = i·log(λ)/T（主对数，Re(ω) 折进第一时间布里渊区）。
- 能带折叠：Re(ω) 与 k 均折进第一布里渊区 [−Ω/2, Ω/2] × [−g/2, g/2]。

## 方法适用范围 / Method Scope

| 方法 | 适用范围 | 输出 |
|------|---------|------|
| PWE | 空间均匀、时间周期、无损、严格正 ε、μ=1 的普通介质 | 两条复准频率能带的实部与虚部 |
| TMM | 空间均匀、无损、严格正 ε/μ 的两层方波时间晶体 | 两条精确复准频率能带的实部与虚部 |
| FDTD | 任意 (x,t) 双周期介质（全场时域） | 场演化；FFT 能带给出实频率脊线，带隙虚部需 TMM/PWE 补充 |

限制：标量一维、非色散、无损（默认）。PWE 入口专注于光学时间晶体；在 `run_pwe.m` 顶部用 `modulationType='square'` 或 `'sinusoidal'` 选择调制，并设置 `epsHigh/epsLow`、占空比或相位、周期、Fourier 截断和 k 扫描。PWE 不调用 TMM、不计算权重、不返回本征矢，只保留 `bands.k` 和两条 `bands.omega`。材料采样器拒绝零、负数和复数 ε。TMM 只考虑两层方波调制；在 `run_tmm.m` 顶部设置两层的 ε/μ、占空比、周期和 k 扫描。`tmm_bands` 采用时间界面处连续的 `[D;B]` 状态，只返回 `result.k` 和两条 `result.omega`，不再暴露单周期矩阵、Floquet 乘子或带隙掩码。FDTD 的 FFT 能带只给出实频率脊线，无法恢复带隙内的复准频率虚部。详见各入口脚本头部注释与 [V3_REFACTOR_PLAN.md](V3_REFACTOR_PLAN.md)。

## 版本与旧版恢复 / Version & Legacy

- 当前版本 **3.0.1**（见 `VERSION.txt`）。3.0.1 修复 V3.0 验收问题（详见 `repair.md`）：`run_fdtd` 的 Case B 现在是真正的有限样品基准——在样品内部测量 ROI 内对 **E(x,t)** 做 2-D FFT，其功率脊线与精确 TMM 通带交叉验证（Method A 脊线），并辅以无样品参考（入射 k/Ω 覆盖与海绵残余背散射的量化）以及三个样品长度在公共有符号 k 网格上的 1/L 脊线外推。预注册脊线容差（由折叠网格量化间隔导出，**不是**旧周期 0.10/0.30）为验收门，外推到 L→∞ 的残差同样对该容差断言。默认出货的是**扩增 ROI（样品内每侧固定 10 个格点的窗口，ROI 随样品长度 L 增长，约占样品 98–99%）**类；**固定中心 ROI** 类——中心短 ROI 无法把体平面波模式与样品-背景杂散模式分离，属**如实记录的物理预期失败**，不作验收门——连同两组记录时长、两类 ROI/窗函数、dx/dt 网格收敛与 Case C 对 PWE 的定量对照，按仓库扁平结构约定放入一次性回归脚本（存于 `/tmp`，不入库）执行，并在控制台与本文件如实披露，不掩盖。这是 3.0.1 中**唯一的定量验收门**。Case A（波包）与 Case C（双重折叠谱）为展示用途，在 3.0.1 中不作定量验收，其原生分辨率在控制台如实报告。Case B 采用弱对比 ε(t)=1.3↔1：强对比二元时间晶体会打开动量带隙，其参数增长 Im(ω)·T 在有限时间窗谱中先于能带解析而主导谱——这是如实记录的物理限制，而非被掩盖的效应。
- V3 为破坏性重构：旧目录（`core/`、`tmm/`、`fdtd/`、`topology/`、`demos/`、`examples/`、`tests/`、`learn/`、`docs/`、`reproduction/`、`Simple_FDTD_NathanZechar/`、`output/`）与旧根入口（`startup_stm.m`、`stm_init.m`、`run_all_demos.m`、`stm_run_examples.m`）均已移除。
- 旧版全部代码可从 git 标签 **`legacy-pre-v3-refactor`** 恢复。
- 最小环境：**MATLAB R2020a+**（Base MATLAB，无 Toolbox 依赖）。三个入口脚本从仓库根目录直接完成，无 `addpath`/startup 依赖；本次 3.0.1 实测于 **R2026a**，R2020a 为声明下限（本机未安装，未实测）。所用 API（`tiledlayout`/`nexttile` 为 R2019b+，`exportgraphics` 为 R2020a，带 `print` 回退）均满足 R2020a。

## 参考文献 / References

- J. Park and B. Min, "Spatiotemporal plane wave expansion method for arbitrary space-time periodic photonic media," *Optics Letters* **46**, 484–487 (2021). [DOI: 10.1364/OL.411622](https://doi.org/10.1364/OL.411622)
- D. Ramaccia, A. Alù, A. Toscano, and F. Bilotti, "Temporal multilayer structures for designing higher-order transfer functions using time-varying metamaterials," *Applied Physics Letters* **118**, 101901 (2021). [DOI: 10.1063/5.0042567](https://doi.org/10.1063/5.0042567)
