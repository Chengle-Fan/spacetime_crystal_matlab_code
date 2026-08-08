# 复现：Topological aspects of photonic time crystals (Lustig et al., 2018)

本目录包含对论文 **"Topological aspects of photonic time crystals"** (E. Lustig, Y. Sharabi, and M. Segev, *Optica* 5, 1390–1395, 2018, DOI: [10.1364/OPTICA.5.001390](https://doi.org/10.1364/OPTICA.5.001390)) 的 MATLAB 复现代码。

> **注意：** 本复现针对 *Optica* 正式发表版本。原 arXiv preprint (1803.08731v1) 的图号有所不同——详见下方的图号映射表。

## 论文概述

该论文首次在**光子时间晶体（photonic time-crystal, PTC）**中引入拓扑能带理论。PTC 是一种介电常数在时间上周期性调制、空间上均匀的介质。由于空间平移对称性，动量 k 守恒；由于时间周期性，Floquet 定理适用。

核心发现：
- **图1** 展示了二元PTC的介电常数示意图和Floquet能带结构（含Zak相位 0 或 π）
- **图2** 展示了FDTD模拟的脉冲演化：带内（4个出射脉冲）vs 带隙（指数增长 + 2个出射脉冲）
- **图4** 展示了前6个动量带隙中前向/后向Floquet模式间的相对相位，验证了拓扑符号律（Eq. 6）
- **图5** 展示了两个不同拓扑的PTC之间的时间界面态——一种"时间拓扑边缘态"

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
Ω(k) = (1/T) cos⁻¹(W - X)
```
其中 |W - X| ≤ 1 时 Ω 为实数（能带），否则为复数（带隙）。

**Zak相位（Eq. 5）：**
```
θ_m^Zak = ∫_{-π/T}^{π/T} dΩ [ i ∫_0^T dt ε(t) D*_{m,Ω}(t) ∂_Ω D_{m,Ω}(t) ]
```
由时间反演对称性（ε(t)=ε(-t)）量子化为0或π。

**拓扑相位符号律（Eq. 6）：**
```
sgn(φ_s) = δ (-1)^s (-1)^l exp(i Σ_{m=1}^{s-1} θ_m^Zak)
```

## 图号映射（arXiv v1 → Optica 正式版）

| arXiv v1 | Optica 正式版 | 描述 |
|----------|--------------|------|
| Fig. 1(a) | Fig. 1(a) | PTC介电常数示意图 |
| Fig. 1(b) | Fig. 1(b) | Floquet能带结构 + Zak相位 |
| Fig. 1(c) | Fig. 2(a) | FDTD：带内脉冲演化（4个出射脉冲） |
| Fig. 1(d) | Fig. 2(b) | FDTD：带隙中脉冲演化（2个出射脉冲 + 指数增长） |
| — | Fig. 3 | 概念示意图（空间光子晶体 vs PTC）— *非计算性* |
| Fig. 3(a-f) | Fig. 4(a-f) | 前6个带隙的相对相位 |
| — | Fig. 5(a-c) | 时间拓扑边缘态（两个PTC的界面） |

## 文件说明

| 文件 | 内容 | 对应论文图 |
|------|------|-----------|
| `fig1_ptc_bands.m` | 二元PTC的ε(t)示意图 + Floquet能带结构 + Zak相位 | **Fig. 1(a,b)** |
| `fig2_fdtd_simulations.m` | FDTD模拟：带内脉冲 + 带隙脉冲 | **Fig. 2(a,b)** |
| `fig4_relative_phase.m` | 前6个带隙的相对相位 φ = arg(E⁻/E⁺) + Eq.(6)验证 | **Fig. 4(a-f)** |
| `fig5_temporal_edge_state.m` | 时间拓扑边缘态：双PTC界面 + 局域化时间峰 | **Fig. 5(a-c)** |
| `compute_zak_phases.m` | 独立Zak相位计算（Wilson-loop诊断 + Eq.(6)验证） | Fig. 1(b) + Fig. 4 |
| `run_all_reproductions.m` | 主脚本：依次运行以上所有复现 | 全部 |

### 历史文件（arXiv v1 版本，已保留供参考）
| 文件 | 原对应图 |
|------|---------|
| `fig1a_ptc_schematic.m` | v1 Fig. 1(a) |
| `fig1b_band_structure.m` | v1 Fig. 1(b) |
| `fig1c_fdtd_in_band.m` | v1 Fig. 1(c) |
| `fig1d_fdtd_in_gap.m` | v1 Fig. 1(d) |
| `fig3_relative_phase.m` | v1 Fig. 3(a-f) |

## 运行方法

### 环境要求
- MATLAB R2020a 或更高版本
- Signal Processing Toolbox（仅用于频谱分析，非必需）

### 运行步骤

1. **单个图复现**：
   ```matlab
   % 在MATLAB中切换到本目录，然后运行：
   run('fig1_ptc_bands')            % 快速（能带计算）
   run('fig2_fdtd_simulations')     % 耗时较长（FDTD，~6-20分钟）
   run('fig4_relative_phase')       % 中等计算量（TMM + Zak相位）
   run('fig5_temporal_edge_state')  % 中等（域壁模式分析）
   ```

2. **全部复现**（推荐）：
   ```matlab
   run('run_all_reproductions')
   ```

3. **输出位置**：所有图片和中间数据保存在 `output/` 子目录中。

### 注意事项
- **Fig. 2 的 FDTD 模拟较耗时**（每次约 3-10 分钟），请耐心等待。可在脚本中设置 `doFDTD = false` 以加载已保存的数据。
- **Fig. 5** 默认使用解析域壁模式分析（快速），设置 `doFDTD = true` 可运行完整的 FDTD 模拟（较慢）。
- 所有脚本会自动添加到父目录的STM工具箱路径（通过`startup_stm.m`）。

## 核心计算方法

### 能带结构计算（Fig. 1b）
1. 对每个k值，计算monodromy矩阵 U(k)
2. 对角化 U(k) 得到Floquet特征值 λ(k)
3. Floquet频率 Ω(k) = cos⁻¹(W-X)/T
4. 带边条件：|W-X| = 1

### Zak相位计算
1. 在每个能带区间内，跟踪U(k)的双正交本征态
2. 计算相邻k点间的本征态重叠（Wilson链）
3. Zak相位 = -arg(Π_j ⟨v_L(k_j)|v_R(k_{j+1})⟩)，量子化为 0 或 π

### FDTD模拟（Fig. 2）
1. 使用 D/B-Yee 蛙跳格式（D和B在时间界面上连续）
2. 介电常数在PTC窗口内按二元周期调制
3. 初始条件：高斯包络的平面波脉冲
4. 记录 |D(x,t)| 并可视化

### 相对相位提取（Fig. 4）
1. **TMM方法**：使用 `temporal_finite_crystal_response` 直接计算有限周期PTC的正向/反向输出振幅比
2. 绘制每个带隙中 φ(k) = arg(E⁻/E⁺) 的曲线
3. 验证拓扑符号律 Eq. (6)

### 时间拓扑边缘态（Fig. 5）
1. **域壁模式分析**：使用 `temporal_domain_wall_mode` 在两个不同拓扑的PTC界面处匹配增长/衰减Floquet本征态
2. 结果：在界面处产生局域化的时间振幅峰——类似空间拓扑边缘态的时间版本
3. 可选FDTD验证

## 代码依赖

本复现代码依赖父目录的STM工具箱：
- `tmm/temporal_crystal_bands.m` — 能带结构
- `tmm/temporal_crystal_monodromy.m` — monodromy矩阵
- `tmm/temporal_finite_crystal_response.m` — 有限周期PTC响应
- `tmm/temporal_domain_wall_mode.m` — 时间域壁模式分析
- `tmm/temporal_db_to_directional.m` — D/B → 定向分量转换
- `fdtd/fdtd1d_db.m` — FDTD求解器

## 参考文献

1. E. Lustig, Y. Sharabi, and M. Segev, "Topological aspects of photonic time crystals," *Optica* 5, 1390–1395 (2018). DOI: [10.1364/OPTICA.5.001390](https://doi.org/10.1364/OPTICA.5.001390)
2. arXiv preprint: [1803.08731v1](https://arxiv.org/abs/1803.08731) (2018) — *注意：图号与正式发表版不同*
3. F. R. Morgenthaler, "Velocity Modulation of Electromagnetic Waves," IRE Trans. Microwave Theory Tech. 6, 167 (1958).
4. J. Park and B. Min, "Spatiotemporal plane wave expansion method for arbitrary space-time periodic photonic media," Opt. Lett. 46, 484 (2021).

## 复现状态

- [x] Fig. 1(a,b): PTC示意图 + Floquet能带结构 + Zak相位
- [x] Fig. 2(a,b): FDTD模拟（带内 + 带隙）
- [x] Fig. 4(a-f): 前6个带隙的相对相位 + 拓扑验证
- [x] Fig. 5(a-c): 时间拓扑边缘态
