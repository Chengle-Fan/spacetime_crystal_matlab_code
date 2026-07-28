# Spacetime Media MATLAB Toolbox / 时空介质 MATLAB 工具箱

光子时空晶体 (Photonic Space-Time Crystals) 的数值仿真与理论计算工具。

A MATLAB toolkit for analyzing and simulating **one-dimensional spacetime-periodic media** — photonic time crystals (PTC) and spacetime crystals.

---

## 快速开始 / Quick Start

### 第一次运行 / First Run

```matlab
startup_stm        % 设置路径
test_smoke         % 自检测试
```

### 建议的学习路径 / Suggested Learning Path

**阶段 A — 时间界面 / Temporal Interfaces**
- `demos/demo05_fdtd_temporal_interface.m` — 单时间界面波分裂
- `demos/demo11_coherent_time_interface.m` — 双输入相干控制

**阶段 B — 光子时间晶体 / Photonic Time Crystals**
- `demos/demo03_ptc_bands_pwe_vs_tmm.m` — PWE 与 TMM 交叉验证
- `demos/demo09_finite_ptc_order_and_phase.m` — 有限周期响应
- `demos/demo12_ptc_convergence_audit.m` — 收敛审计

**阶段 C — 拓扑 / Topology**
- `demos/demo10_temporal_domain_wall_mode.m` — 时间畴壁模式

**阶段 D — 一般时空周期介质 / General ST Media**
- `demos/demo01_reproduce_fig2_stpwe.m` — 论文图 2 复现
- `demos/demo02_complex_frequency_and_momentum_gaps.m` — 复谱分析
- `demos/demo06_fdtd_spacetime_wavepacket.m` — 波包演化
- `demos/demo07_zak_phase_stpwe.m` — 空间 Zak 相位

运行全部示例 / Run all demos:
```matlab
run_all_demos(false)   % 快速
run_all_demos(true)    % 完整（含论文级精度）
```

---

## 代码结构 / Code Structure

| 目录 | 用途 |
|------|------|
| `core/` | ST-PWE 核心（Fourier 域求解） |
| `tmm/` | 时间传输矩阵法（精确求解时间多层结构） |
| `fdtd/` | 一维 D/B-Yee 时域有限差分 |
| `topology/` | 拓扑不变量（Zak 相位、Chern 数） |
| `demos/` | 12 个论文复现演示脚本（v2） |
| `examples/` | 8 个示例脚本（v1） |
| `tests/` | 测试脚本 |

### 版本说明 / Version Notes

本仓库合并了两个版本：

- **v1** (原始): 根目录及 `examples/`，含基础框架
- **v2** (增强): 新增 `demos/` 及 12 个复现论文图表的演示脚本，扩展 `core/`、`tmm/`，增加 PWE 收敛审计、时间畴壁、双输入相干界面等功能

---

## 文档导航 / Documentation

| 文档 | 说明 |
|------|------|
| [docs/manual.md](docs/manual.md) | **完整中文用户手册** — 参数设计、能带计算、拓扑不变量、TMM/FDTD 仿真、函数索引 |
| [docs/conventions.md](docs/conventions.md) | **物理约定与适用域** — 时间因子符号、界面连续量、模型选择边界、常见误判 |
| [docs/paper-map.md](docs/paper-map.md) | **文献覆盖矩阵** — 库中 20 份 PDF 与代码的逐项覆盖状态 |
| [docs/roadmap.md](docs/roadmap.md) | **研究路线图** — 从现有代码出发的课题扩展方向 |
| [docs/validation.md](docs/validation.md) | **Validation record** — v2 包的三层数值验证结果 (English) |

### 核心参考论文 / References

- J. Park and B. Min, "Spatiotemporal plane wave expansion method for arbitrary space-time periodic photonic media," *Optics Letters* **46**, 484–487 (2021). [DOI: 10.1364/OL.411622](https://doi.org/10.1364/OL.411622)
- D. Ramaccia, A. Alù, A. Toscano, and F. Bilotti, "Temporal multilayer structures for designing higher-order transfer functions using time-varying metamaterials," *Applied Physics Letters* **118**, 101901 (2021). [DOI: 10.1063/5.0042567](https://doi.org/10.1063/5.0042567)
- G. R. Morgenthaler, "Velocity modulation of electromagnetic waves," *IRE Trans. Microwave Theory Tech.* **6**, 167 (1958).
- T. Fukui, Y. Hatsugai & H. Suzuki, "Chern numbers in discretized Brillouin zones," *J. Phys. Soc. Jpn.* **74**, 1674 (2005).
