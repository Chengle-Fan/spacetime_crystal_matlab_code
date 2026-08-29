# Spacetime Media MATLAB — V3 / 时空介质 MATLAB V3 版

一维光子时空晶体（photonic spacetime crystals）的数值仿真与理论计算工具：平面波展开法 (PWE)、时间传输矩阵法 (TMM)、时域有限差分法 (FDTD)。三个入口脚本扁平化位于仓库根目录，无需任何路径设置。

A MATLAB toolkit for 1D scalar photonic spacetime crystals: plane-wave expansion (PWE), temporal transfer-matrix method (TMM), and finite-difference time-domain (FDTD). Three flat entry scripts at the repo root — no `addpath`, no startup file.

---

## 三个执行入口 / The Three Entry Scripts

| 脚本 | 方法 | 内容 |
|------|------|------|
| `run_pwe.m` | 平面波展开 (PWE) | 一般时空晶体 ε(x,t) 的双折叠布里渊区能带；二元时间晶体 PWE/TMM 交叉验证 |
| `run_tmm.m` | 时间传输矩阵法 (TMM) | 二元时间晶体的精确 Floquet 能带、动量带隙与增长/衰减分支 |
| `run_fdtd.m` | 时域有限差分 (FDTD) | 波包演化 + 宽带 FDTD-FFT 能带（对 TMM 基准）；一般时空晶体的双重折叠谱 |

在仓库根目录直接运行：

```matlab
run_pwe      % PWE：一般时空晶体 + PTC，双折叠第一布里渊区能带
run_tmm      % TMM：二元时间晶体精确能带
run_fdtd     % FDTD：波包、宽带 FFT 能带、双重折叠谱
```

结果留在工作区；默认不写盘、不创建目录（各脚本末尾的 `doSave` 开关默认关闭）。

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
| PWE | 任意 (x,t) 双周期、无损、标量介质 | 复准频率能带（含带隙内增长/衰减分支） |
| TMM | 空间均匀、时间分段常数（二元/多层）介质 | 精确复准频率（FDTD-FFT 与 PWE 的基准） |
| FDTD | 任意 (x,t) 双周期介质（全场时域） | 场演化；FFT 能带给出实频率脊线，带隙虚部需 TMM/PWE 补充 |

限制：标量一维、非色散、无损（默认）。FDTD 的 FFT 能带只给出实频率脊线（功率谱），无法恢复带隙内的复准频率虚部——带隙结构请以 TMM/PWE 为准。详见各入口脚本头部注释与 [V3_REFACTOR_PLAN.md](V3_REFACTOR_PLAN.md)。

## 版本与旧版恢复 / Version & Legacy

- 当前版本 **3.0.0**（见 `VERSION.txt`）。V3 为破坏性重构：旧目录（`core/`、`tmm/`、`fdtd/`、`topology/`、`demos/`、`examples/`、`tests/`、`learn/`、`docs/`、`reproduction/`、`Simple_FDTD_NathanZechar/`、`output/`）与旧根入口（`startup_stm.m`、`stm_init.m`、`run_all_demos.m`、`stm_run_examples.m`）均已移除。
- 旧版全部代码可从 git 标签 **`legacy-pre-v3-refactor`** 恢复。
- 最小环境：**MATLAB R2020a+**（Base MATLAB，无 Toolbox 依赖）。

## 参考文献 / References

- J. Park and B. Min, "Spatiotemporal plane wave expansion method for arbitrary space-time periodic photonic media," *Optics Letters* **46**, 484–487 (2021). [DOI: 10.1364/OL.411622](https://doi.org/10.1364/OL.411622)
- D. Ramaccia, A. Alù, A. Toscano, and F. Bilotti, "Temporal multilayer structures for designing higher-order transfer functions using time-varying metamaterials," *Applied Physics Letters* **118**, 101901 (2021). [DOI: 10.1063/5.0042567](https://doi.org/10.1063/5.0042567)
