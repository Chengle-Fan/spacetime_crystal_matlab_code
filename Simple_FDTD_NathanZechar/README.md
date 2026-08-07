# Simple FDTD — Nathan Zechar

**作者：** Nathan Zechar, Wright State University, 2021  
**许可协议：** BSD 3-Clause  
**语言：** MATLAB

---

## 概述

本仓库包含一套用于教学和入门学习的**时域有限差分法（FDTD, Finite-Difference Time-Domain）**  MATLAB 代码，分别实现了**一维（1D）、二维（2D）和三维（3D）**的电磁波传播仿真。该代码基于 **Yee 元胞（Yee Cell）** 算法，直接在时域中离散化麦克斯韦旋度方程组，逐步推进电磁场。

代码以**极简风格**编写，专注于清晰展示 FDTD 核心算法，非常适合：
- 学习 FDTD 方法的基本原理
- 理解 Yee 元胞的空间排布与蛙跳（leapfrog）时间推进
- 作为更复杂电磁仿真项目的基础模板

---

## 文件结构

```
Simple_FDTD_NathanZechar/
├── license.txt                          # BSD 3-Clause 许可证
└── Simple_FDTD_NathanZechar/
    ├── license.txt                      # BSD 3-Clause 许可证（副本）
    └── Simple_FDTD_NathanZechar/
        ├── FDTD_1D.m                    # 一维 FDTD 仿真（Ex-Hy 模式）
        ├── FDTD_2D.m                    # 二维 FDTD 仿真（TMz 模式）
        └── FDTD_3D_Folder/
            ├── FDTD_3D.m                # 三维 FDTD 仿真（全 6 分量）
            ├── vars2Single.m            # 将变量转换为单精度浮点
            └── gpu2Single.m             # 将变量迁移至 GPU（gpuArray）
```

---

## 各文件详解

### 1. [FDTD_1D.m](Simple_FDTD_NathanZechar/Simple_FDTD_NathanZechar/FDTD_1D.m) — 一维 FDTD

**物理模型：** 一维空间中沿 x 方向传播的平面电磁波，场分量为 $E_x$ 和 $H_y$。

**主要参数：**

| 参数 | 含义 | 默认值 |
|------|------|--------|
| `f0` | 激励源频率 | 1 MHz |
| `Lf` | 每波长采样点数 | 10 |
| `Lx` | x 方向波长数 | 20 |
| `nt` | 时间步数 | 1000 |

**边界条件：** 未设置（理想电导体 PEC 边界，即场在边界处为零）。

**激励源：** 高斯包络调制的正弦波点源，表达式为：
$$\text{source}(t) = \sin(2\pi f_0 \cdot dt \cdot t) \cdot \exp\left(-\frac{1}{2}\left(\frac{t-20}{8}\right)^2\right)$$

**更新方程（蛙跳格式）：**

- **H 场更新：**  $H_y^{n+1/2}(i) = H_y^{n-1/2}(i) + \frac{\Delta t}{\varepsilon_0 \Delta x}\left[E_x^n(i+1) - E_x^n(i)\right]$
- **E 场更新：**  $E_x^{n+1}(i) = E_x^n(i) + \frac{\Delta t}{\mu_0 \Delta x}\left[H_y^{n+1/2}(i) - H_y^{n+1/2}(i-1)\right]$

> **注意：** 代码中 E 场和 H 场的更新系数似有互换，建议在学习和使用时注意核对。CFL 稳定条件为 $\Delta t = 0.99 \cdot \Delta x / c_0$。

**可视化：** 实时绘制 $E_x$ 场沿 x 轴的分布曲线。

---

### 2. [FDTD_2D.m](Simple_FDTD_NathanZechar/Simple_FDTD_NathanZechar/FDTD_2D.m) — 二维 FDTD

**物理模型：** 二维 TMz 模式（$E_z$, $H_x$, $H_y$），波在 xy 平面内传播。

**主要参数：**

| 参数 | 含义 | 默认值 |
|------|------|--------|
| `f0` | 激励源频率 | 1 MHz |
| `Lf` | 每波长采样点数 | 10 |
| `Lx, Ly` | x, y 方向波长数 | 各 8 |
| `nt` | 时间步数 | 1000 |

**更新方程：**

- **H 场更新（法拉第定律离散化）：**
  $$H_x^{n+1/2} = H_x^{n-1/2} - \frac{\Delta t}{\mu_0 \Delta y} \frac{\partial E_z}{\partial y}$$
  $$H_y^{n+1/2} = H_y^{n-1/2} + \frac{\Delta t}{\mu_0 \Delta x} \frac{\partial E_z}{\partial x}$$

- **E 场更新（安培定律离散化）：**
  $$E_z^{n+1} = E_z^n + \frac{\Delta t}{\varepsilon_0 \Delta x}\frac{\partial H_y}{\partial x} - \frac{\Delta t}{\varepsilon_0 \Delta y}\frac{\partial H_x}{\partial y}$$

**CFL 条件：** $\Delta t = 0.99 \cdot \frac{1}{c_0\sqrt{\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2}}}$

**激励源：** 与 1D 相同的高斯调制正弦波，位于计算域中心。

**可视化：** 使用 `imagesc` 实时显示 $E_z$ 场的二维彩色图像。

**性能监控：** 每个时间步使用 `tic`/`toc` 输出耗时。

---

### 3. [FDTD_3D.m](Simple_FDTD_NathanZechar/Simple_FDTD_NathanZechar/FDTD_3D.m) — 三维 FDTD

**物理模型：** 完整的三维 Yee 元胞，包含全部 6 个电磁场分量（$E_x, E_y, E_z, H_x, H_y, H_z$）。

**主要参数：**

| 参数 | 含义 | 默认值 |
|------|------|--------|
| `f0` | 激励源频率 | 1 MHz |
| `Lf` | 每波长采样点数 | 10 |
| `Lx, Ly, Lz` | x, y, z 方向波长数 | 各 8 |
| `nt` | 时间步数 | 1000 |
| `single` | 是否使用单精度 | 1（启用） |
| `usegpu` | 是否使用 GPU 加速 | 0（禁用） |

**Yee 元胞场分量排布：**

| 场分量 | 网格维度 | 在元胞中的位置 |
|--------|----------|----------------|
| $E_x$ | `(Nx, Ny+1, Nz+1)` | 棱边中心（x 方向） |
| $E_y$ | `(Nx+1, Ny, Nz+1)` | 棱边中心（y 方向） |
| $E_z$ | `(Nx+1, Ny+1, Nz)` | 棱边中心（z 方向） |
| $H_x$ | `(Nx+1, Ny, Nz)` | 面心（x 方向） |
| $H_y$ | `(Nx, Ny+1, Nz)` | 面心（y 方向） |
| $H_z$ | `(Nx, Ny, Nz+1)` | 面心（z 方向） |

**CFL 条件：** $\Delta t = 0.99 \cdot \frac{1}{c_0\sqrt{\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2} + \frac{1}{\Delta z^2}}}$

**激励源：** 使用 Sigmoid 包络的正弦波，与 1D/2D 的高斯包络不同：
$$\text{source}(t) = \frac{\sin(2\pi f_0 \cdot dt \cdot t)}{1 + \exp(-0.2(t-60))}$$

**可视化：** 使用 `slice` 函数显示 $E_z$ 在三个正交切片平面上的分布，颜色映射范围为 $[-0.001, 0.001]$。

---

### 4. [vars2Single.m](Simple_FDTD_NathanZechar/Simple_FDTD_NathanZechar/FDTD_3D_Folder/vars2Single.m) — 单精度转换

将所有电磁场分量和更新系数从默认的 `double` 精度转换为 `single` 精度，可节省约 **50% 的内存**，并显著提升三维仿真的计算速度。

**使用方法：** 在 `FDTD_3D.m` 中设置 `single = 1`，脚本会自动调用此函数。

---

### 5. [gpu2Single.m](Simple_FDTD_NathanZechar/Simple_FDTD_NathanZechar/FDTD_3D_Folder/gpu2Single.m) — GPU 加速

将所有电磁场分量和更新系数转换为 `gpuArray` 类型，使计算在 GPU 上执行。

**前提条件：**
- 需要 **Parallel Computing Toolbox**
- 需要支持 CUDA 的 **NVIDIA GPU**

**使用方法：** 在 `FDTD_3D.m` 中设置 `usegpu = 1`（同时确保 `single = 1`），脚本会自动调用此函数。

---

## 快速开始

### 系统要求

- **MATLAB** R2018b 或更高版本（推荐）
- （可选）Parallel Computing Toolbox + NVIDIA GPU（仅 3D GPU 加速需要）

### 运行步骤

1. **1D 仿真：**
   ```matlab
   run('FDTD_1D.m')
   ```
   将看到一个动态更新的电场一维传播曲线图。

2. **2D 仿真：**
   ```matlab
   run('FDTD_2D.m')
   ```
   将看到一个动态更新的二维电场彩色图像（圆柱波扩散）。

3. **3D 仿真：**
   ```matlab
   run('FDTD_3D_Folder/FDTD_3D.m')
   ```
   将看到一个三维切片图，显示 $E_z$ 在三个正交平面上的分布。

   - 对于更高性能，可设置 `single = 1` 启用单精度。
   - 如有 GPU，可额外设置 `usegpu = 1` 启用 GPU 加速。

---

## 关键概念

### Yee 元胞（Yee Cell）

FDTD 的核心是 Kane S. Yee 于 1966 年提出的空间离散方案：电场和磁场分量在空间上交错排布，每个电场分量被四个磁场分量环绕（反之亦然）。这种排布自然地满足了法拉第定律和安培定律的积分形式。

### 蛙跳时间推进（Leapfrog Scheme）

时间上，E 场和 H 场交替更新，相差半个时间步长 $\Delta t/2$。这种显式格式无需求解线性方程组，每一步的计算量仅为 $O(N)$（$N$ 为网格点数）。

### CFL 稳定性条件

显式时间推进的稳定性要求时间步长满足：
$$\Delta t \leq \frac{1}{c_0 \sqrt{\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2} + \frac{1}{\Delta z^2}}}$$

本代码采用 $0.99$ 的安全系数以确保稳定性。

---

## 局限性与改进方向

1. **无吸收边界条件（ABC）：** 当前代码使用默认的 PEC 边界（场在边界处为零），会引入反射。实际应用中通常需要 **Mur 吸收边界** 或 **PML（完美匹配层）**。

2. **无介质参数分布：** 当前代码仅模拟真空中传播。添加 $\varepsilon_r(x,y,z)$ 和 $\mu_r(x,y,z)$ 的分布即可模拟介质散射体、波导等结构。

3. **无近场-远场变换：** 无法计算散射截面或辐射方向图。

4. **无色散材料模型：** 不支持 Drude、Lorentz、Debye 等色散介质模型。

5. **E/H 场更新系数：** 1D 代码中 E 场和 H 场的更新系数可能存在互换，建议初学者对照标准 Yee 格式仔细核验。

---

## 参考文献

1. Yee, K. S. (1966). *Numerical solution of initial boundary value problems involving Maxwell's equations in isotropic media.* IEEE Transactions on Antennas and Propagation, 14(3), 302–307.

2. Taflove, A., & Hagness, S. C. (2005). *Computational Electrodynamics: The Finite-Difference Time-Domain Method* (3rd ed.). Artech House.

3. Sullivan, D. M. (2013). *Electromagnetic Simulation Using the FDTD Method* (2nd ed.). IEEE Press.

---

## 致谢

感谢 **Nathan Zechar**（Wright State University, 2021）编写并开源此教学代码。
