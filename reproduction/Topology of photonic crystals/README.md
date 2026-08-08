# 复现：Topology of photonic time-crystals (Lustig et al., 2018)

本目录包含对论文 **"Topology of photonic time-crystals"** (E. Lustig, Y. Sharabi, and M. Segev, arXiv:1803.08731v1, 2018) 中图1和图3的MATLAB复现代码。

## 论文概述

该论文研究了**二元光子时间晶体（photonic time-crystal, PTC）**的拓扑性质。光子时间晶体是一种介电常数在时间上周期性调制、空间上均匀的介质。由于空间平移对称性，动量 k 是守恒量；由于时间周期性，Floquet 定理适用。

核心发现：
- **图1**展示了二元PTC的Floquet能带结构、Zak相位以及FDTD模拟的脉冲演化
- **图3**展示了时间折射和反射波之间的相对相位是拓扑可观测量的实验/数值验证

## 物理模型

### 参数
| 参数 | 符号 | 值 |
|------|------|-----|
| 介电常数1 | ε₁ | 3 |
| 介电常数2 | ε₂ | 1 |
| 调制周期 | T | 2π（Ω=1） |
| 片段持续时间 | t₁ = t₂ | T/2（等占空比） |
| 背景介电常数 | ε_bg | 2 |
| 时间反演对称点 | t=0 | 片段1的中点 |

### 关键方程

**Floquet色散关系：**
```
ω_F(k) = (i/T)·ln(λ(k))
```
其中 λ(k) 是单周期monodromy矩阵 U(k) 的特征值。当 |Tr(U)/2| ≤ 1 时 ω_F 为实数（能带），否则为复数（带隙）。

**Zak相位（Eq. 5）：**
```
θ_n^Zak = ∮_BZ i ⟨u_{n,k}| ∂/∂k |u_{n,k}⟩ dk
```
由时间反演对称性（ε(t)=ε(-t)）量子化为0或π。

**拓扑相位符号律（Eq. 6）：**
```
sgn(φ_s) = (-1)^s · (-1)^l · exp(i Σ_{m=1}^{s-1} θ_m^Zak)
```

## 文件说明

| 文件 | 内容 | 对应论文图 |
|------|------|-----------|
| `fig1a_ptc_schematic.m` | 二元PTC的介电常数ε(t)示意图 | Fig. 1(a) |
| `fig1b_band_structure.m` | Floquet能带结构 ω_F(k) + Zak相位标注 | Fig. 1(b) |
| `fig1c_fdtd_in_band.m` | FDTD模拟：能带内的脉冲演化（4个出射脉冲） | Fig. 1(c) |
| `fig1d_fdtd_in_gap.m` | FDTD模拟：带隙内的脉冲演化（2个出射脉冲+指数增长） | Fig. 1(d) |
| `fig3_relative_phase.m` | 前6个带隙的相对相位 φ = arg(E⁻/E⁺) | Fig. 3(a-f) |
| `compute_zak_phases.m` | **独立Zak相位计算**：详细的Wilson-loop诊断+Eq. (6)验证 | Fig. 1(b) + 3 |
| `run_all_reproductions.m` | 主脚本：依次运行以上所有复现 | 全部 |

## 运行方法

### 环境要求
- MATLAB R2020a 或更高版本
- Signal Processing Toolbox（仅用于`stm_stft.m`的频谱分析，非必需）

### 运行步骤

1. **单个图复现**：
   ```matlab
   % 在MATLAB中切换到本目录，然后运行：
   run('fig1a_ptc_schematic')
   run('fig1b_band_structure')
   run('fig1c_fdtd_in_band')       % 耗时较长（FDTD模拟）
   run('fig1d_fdtd_in_gap')        % 耗时较长（FDTD模拟）
   run('fig3_relative_phase')      % 计算量中等
   ```

2. **全部复现**（推荐）：
   ```matlab
   run('run_all_reproductions')
   ```

3. **输出位置**：所有图片和中间数据保存在 `output/` 子目录中。

### 注意事项
- **Fig. 1c 和 Fig. 1d 的FDTD模拟较耗时**（每次约需几分钟），请耐心等待。
- Fig. 3的`doFDTD`标志默认为`false`，仅使用TMM精确解计算相对相位。如需FDTD验证，将脚本中的`doFDTD = false`改为`true`。
- 所有脚本会自动添加到父目录的STM工具箱路径（通过`startup_stm.m`）。

## 核心计算方法

### 能带结构计算（Fig. 1b）
1. 对每个k值，计算monodromy矩阵 U(k) = Π_m exp(-i·k·G_m·t_m)
2. 对角化 U(k) 得到Floquet特征值 λ(k)
3. Floquet频率 ω_F(k) = i·ln(λ)/T
4. 带边条件：|Tr(U)/2| = 1

### Zak相位计算
1. 在每个能带区间内，跟踪U(k)的双正交本征态
2. 计算相邻k点间的本征态重叠（Wilson链）
3. Zak相位 = -arg(Π_j ⟨v_L(k_j)|v_R(k_{j+1})⟩)

### FDTD模拟（Fig. 1c-d）
1. 使用D/B-Yee蛙跳格式（D和B在时间界面上连续）
2. 介电常数在PTC窗口内按二元周期调制
3. 初始条件：高斯包络的平面波脉冲
4. 记录|E(x,t)|并可视化

### 相对相位提取（Fig. 3）
1. **TMM方法**：使用`temporal_finite_crystal_response`直接计算有限周期PTC的正向/反向输出振幅比
2. **FDTD方法**（可选）：对PTC结束后的场做空间傅里叶变换，通过定向分解提取E⁺和E⁻分量
3. 绘制每个带隙中φ(k) = arg(E⁻/E⁺)的曲线

## 代码依赖

本复现代码依赖父目录的STM工具箱：
- `tmm/temporal_crystal_bands.m` — 能带结构
- `tmm/temporal_crystal_monodromy.m` — monodromy矩阵
- `tmm/temporal_finite_crystal_response.m` — 有限周期PTC响应
- `tmm/temporal_db_to_directional.m` — D/B → 定向分量转换
- `fdtd/fdtd1d_db.m` — FDTD求解器

## 参考文献

1. E. Lustig, Y. Sharabi, and M. Segev, "Topology of photonic time-crystals," arXiv:1803.08731v1 (2018).
2. F. R. Morgenthaler, "Velocity Modulation of Electromagnetic Waves," IRE Trans. Microwave Theory Tech. 6, 167 (1958).
3. J. Park and B. Min, "Spatiotemporal plane wave expansion method for arbitrary space-time periodic photonic media," Opt. Lett. 46, 484 (2021).

## 复现状态

- [x] Fig. 1(a): PTC示意图
- [x] Fig. 1(b): Floquet能带结构 + Zak相位
- [x] Fig. 1(c): FDTD能带内脉冲
- [x] Fig. 1(d): FDTD带隙中脉冲
- [x] Fig. 3(a-f): 前6个带隙的相对相位
