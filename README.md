# Spacetime Media MATLAB Toolbox

时空介质 MATLAB 工具箱 — 光子时空晶体 (Photonic Space-Time Crystals) 的数值仿真与理论计算工具。

## 版本说明

本仓库合并了两个版本：

### v1 (原始版本)
位于根目录及 `examples/` 目录下，包含基础框架：
- `stm_init.m` / `stm_run_examples.m` — 入口脚本
- `core/` — 核心函数（PWE、Fourier 调制 slab 等）
- `examples/` — 7 个示例
- `tmm/` — 传输矩阵方法
- `topology/` — 拓扑不变量计算
- `fdtd/` — 时域有限差分仿真

### v2 (增强版本)
新增 `demos/` 目录及增强功能：
- `startup_stm.m` / `run_all_demos.m` — 新版入口脚本
- `demos/` — 12 个复现论文图表的演示脚本
- `core/` — 新增 Fig2 复现相关函数
- `tmm/` — 新增域壁模式、有限晶体响应等函数
- 增强文档：`PAPER_CODE_MAP_CN.md`、`RESEARCH_ROADMAP_CN.md`、`START_HERE_CN.md`

## 快速开始

```matlab
% v1 入口
stm_init
stm_run_examples

% v2 入口
startup_stm
run_all_demos
```

## 参考

- [Phys. Rev. X 4, 021017](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.4.021017)
