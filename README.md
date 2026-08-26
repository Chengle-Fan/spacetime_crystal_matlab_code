# Spacetime Media MATLAB Toolbox / 时空介质 MATLAB 工具箱

光子时空晶体 (Photonic Space-Time Crystals) 的数值仿真与理论计算工具。

A MATLAB toolkit for analyzing and simulating **one-dimensional spacetime-periodic media** — photonic time crystals (PTC) and spacetime crystals.

---

## 快速开始 / Quick Start

```matlab
startup_stm        % 设置路径
test_smoke         % 自检测试
run_all_demos(false)   % 运行全部示例（快速）
```

首次使用前先设置路径并跑通自检测试；看到 `All smoke tests passed.` 即环境就绪。

---

## 学习路径 / Learning Path

**新手从这里开始** → [learn/README.md](learn/README.md)：11 章中文教程，从 MATLAB 入门到动手研究——

> 物理背景 → 三种数值方法 → ST-PWE / 时间 TMM / FDTD 精讲 → 拓扑不变量 → 12 个 demo 导读 → 练习 → 排错。

---

## 代码结构 / Code Structure

| 目录 | 用途 |
|------|------|
| `core/` | ST-PWE 核心（Fourier 域求解） |
| `tmm/` | 时间传输矩阵法（精确求解时间多层结构） |
| `fdtd/` | 一维 D/B-Yee 时域有限差分 |
| `topology/` | 拓扑不变量（Zak 相位、Chern 数） |
| `demos/` | 12 个论文复现演示脚本（v2） |
| `examples/` | 10 个示例脚本（v1） |
| `tests/` | 测试脚本 |

### 版本说明 / Version Notes

本仓库合并了两个版本：

- **v1** (原始): 根目录及 `examples/`，含基础框架
- **v2** (增强): 新增 `demos/` 及 12 个复现论文图表的演示脚本，扩展 `core/`、`tmm/`，增加 PWE 收敛审计、时间畴壁、双输入相干界面等功能

---

## 文档导航 / Documentation

| 文档 | 说明 |
|------|------|
| [learn/README.md](learn/README.md) | **教程（主入口）** — 从零上手的 11 章中文学习路径 |
| [docs/tool-reference.md](docs/tool-reference.md) | **函数工具书** — 全部 113 个 `.m` 文件的函数签名、输入/输出参数与功能说明 |
| [docs/paper-map.md](docs/paper-map.md) | **文献覆盖矩阵** — 库中 20 份 PDF 与代码的逐项覆盖状态 |
| [docs/roadmap.md](docs/roadmap.md) | **研究路线图** — 从现有代码出发的课题扩展方向 |
| [docs/validation.md](docs/validation.md) | **Validation record** — v2 包的三层数值验证结果 (English) |

---

### 核心参考论文 / References

- J. Park and B. Min, "Spatiotemporal plane wave expansion method for arbitrary space-time periodic photonic media," *Optics Letters* **46**, 484–487 (2021). [DOI: 10.1364/OL.411622](https://doi.org/10.1364/OL.411622)
- D. Ramaccia, A. Alù, A. Toscano, and F. Bilotti, "Temporal multilayer structures for designing higher-order transfer functions using time-varying metamaterials," *Applied Physics Letters* **118**, 101901 (2021). [DOI: 10.1063/5.0042567](https://doi.org/10.1063/5.0042567)
- G. R. Morgenthaler, "Velocity modulation of electromagnetic waves," *IRE Trans. Microwave Theory Tech.* **6**, 167 (1958).
- T. Fukui, Y. Hatsugai & H. Suzuki, "Chern numbers in discretized Brillouin zones," *J. Phys. Soc. Jpn.* **74**, 1674 (2005).
