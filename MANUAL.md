# STM 工具包用户手册：时空光子晶体计算

> **STM (Space-Time Media) Toolbox**
> 一维时空调制介质的光子能带结构、拓扑不变量和时域仿真

---

## 目录

1. [概述](#1-概述)
2. [快速开始](#2-快速开始)
3. [参数设计：定义你的时空晶体](#3-参数设计定义你的时空晶体)
4. [构建系统矩阵](#4-构建系统矩阵)
5. [能带结构计算](#5-能带结构计算)
6. [拓扑不变量计算](#6-拓扑不变量计算)
7. [场分布图绘制](#7-场分布图绘制)
8. [TMM 交叉验证](#8-tmm-交叉验证)
9. [FDTD 时域仿真验证](#9-fdtd-时域仿真验证)
10. [完整科研流程示例](#10-完整科研流程示例)
11. [函数索引](#11-函数索引)

---

## 1. 概述

### 1.1 工具包能做什么

本 MATLAB 工具包用于分析和仿真 **一维时空周期性介质** —— 同时具有空间周期（`Λ`）和时间周期（`T`）的电磁介质。这种介质也称为**光子时间晶体**（Photonic Time Crystal, PTC）或**时空晶体**（Spacetime Crystal）。

工具包内包含三种独立的方法，互为验证：

| 方法 | 目录 | 适用场景 |
|------|------|---------|
| **ST-PWE** | `core/` | 傅里叶域求解，快速得到全频带结构，支持复数频率/波矢分析 |
| **TMM** | `tmm/` | 时间转移矩阵，精确求解多层时间结构，验证复数能带 |
| **FDTD** | `fdtd/` | 时域仿真，直观展示波包演化，验证理论预测 |

此外还提供拓扑计算工具（`topology/`）用于计算 Zak 相位和 Chern 数。

### 1.2 核心物理

时空调制介电常数：

$$ \varepsilon(x,t) = \sum_{m,n} \varepsilon_{mn}\, e^{i n g x - i m \Omega t} $$

其中 $g = 2\pi/\Lambda$ 为空间倒格矢，$\Omega = 2\pi/T$ 为时间调制频率。电场采用 Floquet-Bloch 展开：

$$ E(x,t) = e^{i(kx-\omega t)} \sum_{m,n} E_{mn}\, e^{i(n g x - m \Omega t)} $$

两种特征值问题对应两种禁带：

| 视角 | 特征值方程 | 禁带类型 | 物理意义 |
|------|-----------|----------|---------|
| 固定 $k$ → 求 $\omega$ | $A·v = \omega·B·v$ | **动量禁带** (momentum gap) | $\text{Im}(\omega)\neq 0$ → 时间指数增长/衰减（参量放大） |
| 固定 $\omega$ → 求 $k$ | $A·v = k·v$ | **频率禁带** (frequency gap) | $\text{Im}(k)\neq 0$ → 空间指数衰减（倏逝波） |

---

## 2. 快速开始

### 2.1 设置 MATLAB 路径

```matlab
>> stm_init();
```

将 `core/`, `tmm/`, `fdtd/`, `topology/`, `examples/`, `tests/` 全部加入路径，并在根目录下创建 `output/` 文件夹。

### 2.2 运行所有示例

```matlab
>> stm_run_examples(false)    % 快速示例（5 个脚本，几分钟内完成）
>> stm_run_examples(true)     % 完整示例（含论文级精度计算，约需 10-30 分钟）
```

每个示例的输出图片保存在 `output/` 文件夹中。
也可以单独调用：

```matlab
>> example_floquet_bands_and_fields('paper')    % 重现论文 Fig.2（论文级精度）
>> example_floquet_bands_and_fields('quick')     % 快速版本用于验证
>> example_complex_gaps()
>> example_ptc_pwe_vs_tmm()
>> example_tmm_multilayer()
>> example_fdtd_interface()
>> example_fdtd_wavepacket()
>> example_zak_phase()
>> example_chern_pump()
```

### 2.3 自检测试

```matlab
>> stm_selftest()
```

验证傅里叶系数一致性、TMM 合成正确性、单值矩阵保面积、Chern 数量子化和 FDTD 稳定性。

---

## 3. 参数设计：定义你的时空晶体

这是整个流程的起点。你需要做三件事：
1. 定义介电常数 $\varepsilon(x,t)$ 的空间分布和时空调制方式
2. （可选）定义磁导率 $\mu(x,t)$
3. 设定截断阶数 $N_{\text{space}}$ 和 $M_{\text{time}}$

### 3.1 方式 A：提供解析傅里叶系数（推荐，精度最高）

直接提供 $\varepsilon_{mn}$ 的解析表达式。以 Park-Min Fig.2 为例，介电常数在一个周期内的构型：

| $x$ 范围 | $\varepsilon(x,t)$ |
|----------|-------------------|
| $[0, \frac{3}{4}\Lambda)$ | $\varepsilon_1 = 2$（静态） |
| $[\frac{3}{4}\Lambda, \Lambda)$ | $\varepsilon_c[1 + \text{modDepth}\cdot\sin(\Omega t)]$（时空调制） |

#### 步骤 1：加载预设参数

```matlab
% 使用内置 Park-Min 材料预设
p = stm_preset_modulated_slab();
% p 包含: Lambda, c0, g, eps1, epsc, modDepth, xModStart, xModEnd, OmegaBar, Omega, T
```

> **关键参数设置建议：**
> - `modDepth` ∈ (0, 1)：调制深度越大，动量禁带越宽
> - `OmegaBar`：调制频率与空间布里渊区边界的频率之比，决定禁带位置
> - `Nspace` 和 `Mtime`：截断阶数，决定计算精度。对于正弦调制 `Mtime=1` 已经精确；对于方波调制需要更大的 `Mtime`

#### 步骤 2：编写 $\varepsilon(x,t)$ 的实空间函数

参考内置的 `stm_permittivity_modulated_slab.m`：

```matlab
function epsr = stm_permittivity_modulated_slab(x, t, p)
    % 输入：x, t 可以是标量或兼容数组
    % 返回值：与 x,t 同尺寸的介电常数分布
    xCell = mod(x, p.Lambda);
    mask = (xCell >= p.xModStart) & (xCell < p.xModEnd);
    epsr = p.eps1 * ones(size(x));
    epsr(mask) = p.epsc * (1 + p.modDepth * sin(p.Omega * t));
end
```

#### 步骤 3：提供解析傅里叶系数

参考内置的 `stm_fourier_modulated_slab.m`：

```matlab
function coeff = stm_fourier_modulated_slab(m, n, p)
    % ε_{mn} = (1/Λ)·∫₀^{Λ} (1/T)·∫₀^{T} ε(x,t)·e^{-i n g x + i m Ω t} dx dt
    % 在调制区域和静态区域分别作积分
    staticPart = stm_interval_fourier(n, 0, p.xModStart, p.Lambda);
    modPart = stm_interval_fourier(n, p.xModStart, p.xModEnd, p.Lambda);

    switch m
        case 0                             % 时间直流分量
            coeff = p.eps1*staticPart + p.epsc*modPart;
        case 1                             % 时间基频 (e^{-iΩt})
            coeff = 1i*(p.epsc*p.modDepth/2)*modPart;
        case -1                            % 时间基频 (e^{+iΩt})
            coeff = -1i*(p.epsc*p.modDepth/2)*modPart;
        otherwise
            coeff = 0;                     % 更高阶谐波为 0（正弦调制仅含 ±1 阶）
    end
end
```

#### 辅助工具：空间区间傅里叶积分

```matlab
% 文件：stm_interval_fourier.m
% (1/Λ)·∫_{x_a}^{x_b} exp(-i n g x) dx
function value = stm_interval_fourier(n, xa, xb, Lambda)
    if n == 0
        value = (xb - xa)/Lambda;
    else
        g = 2*pi/Lambda;
        value = (exp(-1i*n*g*xb) - exp(-1i*n*g*xa))/(-1i*n*g*Lambda);
    end
end
```

### 3.2 方式 B：提供实空间函数 + 数值采样（通用，灵活）

对于没有解析傅里叶系数的任意 $\varepsilon(x,t)$ 分布，用数值采样计算系数。

```matlab
% 步骤 1：定义实空间材料函数
materialFun = @(X, T) stm_permittivity_modulated_slab(X, T, p);   % 或者自定义

% 步骤 2：采样得到傅里叶系数
mOrders = -1:1;         % 需要的时间谐波阶数
nOrders = -3:3;         % 需要的空间谐波阶数
table = stpwe_sample_fourier_coefficients(materialFun, ...
    mOrders, nOrders, p.Lambda, p.T, 512, 256);
% 最后两个参数：空间采样点数 Nx，时间采样点数 Nt

% 步骤 3：从采样表查询系数
epsCoeff = @(m, n) stpwe_lookup_coefficient(table, m, n);
```

> **`stpwe_sample_fourier_coefficients` 参数说明：**
>
> | 参数 | 含义 | 建议值 |
> |------|------|--------|
> | `Nx` | 每个空间周期内的采样点数 | `max(128, 8*numel(nOrders))` |
> | `Nt` | 每个时间周期内的采样点数 | `max(128, 8*numel(mOrders))` |
>
> 采样使用均匀端点外网格（避免 Gibbs 效应）。

### 3.3 方式 C：二元（方波型）时间晶体

对于在两种值之间方波切换的简单时间晶体：

```matlab
% 使用 temporal_binary_eps_coeff.m 直接给出解析傅里叶系数
epsA = 1.0;    epsB = 4.0;    % 两种介电常数值
dutyA = 0.5;                    % εA 的占空比

epsCoeff = @(m, n) double(n == 0) * temporal_binary_eps_coeff(m, epsA, epsB, dutyA);
% 注意此处 n==0：空间均匀的时间晶体（无空间调制）
```

> **对比三种方式：**
> | 方式 | 优点 | 缺点 | 适用场景 |
> |------|------|------|---------|
> | A. 解析系数 | 精度最高，计算最快 | 需要推导解析式 | 标准波形（正弦、方波等） |
> | B. 数值采样 | 通用性好 | 需要较高采样率 | 任意复杂波形 |
> | C. 二元解析 | 对间断波形精确 | 收敛慢，需大 Mtime | 方波时间晶体 |

---

## 4. 构建系统矩阵

将傅里叶系数组装成特征值问题所需的大矩阵。

### 4.1 调用方式

```matlab
% 语法
sys = stpwe_build_system(epsCoeff, muCoeff, Nspace, Mtime, g, Omega)

% 例 1：非磁性时空晶体（使用内置 Park-Min 预设）
epsCoeff = @(m,n) stm_fourier_modulated_slab(m, n, p);
sys = stpwe_build_system(epsCoeff, [], 10, 1, p.g, p.Omega);
%                                 ↑↑↑  ↑  Nspace=10, Mtime=1

% 例 2：磁性材料
muCoeff = @(m,n) double(m==0 && n==0);   % 简单的 μ=1 均匀分布
sys = stpwe_build_system(epsCoeff, muCoeff, 12, 2, p.g, p.Omega);
```

### 4.2 参数说明

| 参数 | 含义 | 典型值 |
|------|------|--------|
| `epsCoeff(m,n)` | 返回介电常数 $mn$ 阶傅里叶系数的函数句柄 | 见 3.1 或 3.2 |
| `muCoeff(m,n)` | 返回磁导率 $mn$ 阶傅里叶系数的函数句柄 | `[]` 或 `@(m,n) double(m==0 && n==0)` |
| `Nspace` | 空间截断阶数 $n \in [-N_{\text{space}}, N_{\text{space}}]$ | **静态区 5-15**，时空区 8-20 |
| `Mtime` | 时间截断阶数 $m \in [-M_{\text{time}}, M_{\text{time}}]$ | 正弦调制 1-2，方波调制 15-30 |
| `g` | 空间倒格矢 $2\pi/\Lambda$ | 由周期 Λ 决定 |
| `Omega` | 时间调制角频率 $2\pi/T$ | 由周期 T 决定 |

> **截断阶数选择技巧：**
> - 对于**正弦调制**，$M_{\text{time}}=1$ 已经精确（仅含 ±1 阶傅里叶系数）
> - 对于**方波调制**，$M_{\text{time}}$ 需要 15-30 以获取收敛结果
> - $N_{\text{space}}$ 取 8-12 通常已足够，论文级精度取 20
> - 总矩阵大小为 $2S \times 2S$，其中 $S = (2N_{\text{space}}+1)(2M_{\text{time}}+1)$ → `Nspace=10, Mtime=1`时矩阵为 126×126

### 4.3 sys 结构体字段

| 字段 | 大小 | 含义 |
|------|------|------|
| `sys.Nspace`, `sys.Mtime` | 标量 | 截断阶数 |
| `sys.nList`, `sys.mList` | $S \times 1$ | 每个谐波对应的 $(n,m)$ 索引 |
| `sys.S` | 标量 | 谐波总数 $=(2N_{\text{space}}+1)(2M_{\text{time}}+1)$ |
| `sys.g`, `sys.Omega` | 标量 | 倒格矢和调制频率 |
| `sys.Ceps`, `sys.Cmu` | $S \times S$ | 介电/磁导率卷积矩阵 |
| `sys.G` | $S \times S$ | 空间波数矩阵 $\text{diag}(n·g)$ |
| `sys.W` | $S \times S$ | 时间频率矩阵 $\text{diag}(m·\Omega)$ |
| `sys.I`, `sys.Z` | $S \times S$ | 单位阵和零矩阵 |
| `sys.Bomega` | $2S \times 2S$ | $\begin{bmatrix}0 & C_\mu \\ C_\varepsilon & 0\end{bmatrix}$ |

---

## 5. 能带结构计算

### 5.1 固定 k → 复 ω（动量禁带分析）

这是最常用的计算模式，给定 Bloch 波矢扫描，求解对应的 Floquet 频率。

```matlab
% 设置 k 空间扫描路径
kBar = linspace(-0.5, 0.5, 201);     % 归一化波矢 kBar = k·Λ/(2π)
kValues = p.g * kBar;                 % 物理波矢

% 在每个 k 点求解
allOmega = [];  allK = [];  allWeight = [];  allFold = [];
for ik = 1:numel(kValues)
    sol = stpwe_solve_omega(sys, kValues(ik));

    % 可选：筛选有效模式
    fBar = sol.omega / (p.g * p.c0);          % 归一化频率
    keep = isfinite(fBar) ...                  % 排除无穷大
        & real(fBar) >= 0 & real(fBar) <= fMax ...  % 窗口范围
        & abs(imag(fBar)) <= imagTolerance ...       % 虚部阈值
        & sol.m0Weight >= weightThreshold;           % m=0 占比阈值

    % 存储结果
    allK = [allK; repmat(kBar(ik), sum(keep), 1)];
    allOmega = [allOmega; real(fBar(keep))];
    allWeight = [allWeight; sol.m0Weight(keep)];
    allFold = [allFold; abs(imag(fBar(keep)))];
end

% 绘制能带图（用 m0Weight 着色）
figure('Color','w');
scatter(allK, allOmega, 8, allWeight, 'filled');
xlabel('k\Lambda/(2\pi)');  ylabel('\omega\Lambda/(2\pi c)');
colormap(flipud(gray));  colorbar;
```

> **输出 `sol` 结构体字段：**
>
> | 字段 | 含义 |
> |------|------|
> | `sol.omega` | 所有本征值（复频率），$\omega$ |
> | `sol.R`, `sol.L` | 右/左特征向量矩阵，已归一化为 $L^\dagger B R = I$ |
> | `sol.m0Weight` | 各模式的 m=0 分量占比（越大表示越接近静态模式） |

#### 筛选参数建议

| 筛选参数 | 作用 | 默认值 |
|---------|------|--------|
| `fMax` | 显示的最大频率 | 取决于研究范围 |
| `imagTolerance` | 允许的 |Im(ω)| 上限 | `3e-3` 用于"准传播"模式 |
| `weightThreshold` | 最小 m0 权重 | `0.035` 用于筛选物理模式 |

### 5.2 固定 ω → 复 k（频率禁带分析）

```matlab
fBarSweep = linspace(0.38, 0.43, 121);     % 归一化频率扫描
kStore = [];  fStore = [];  kiStore = [];

for iw = 1:numel(fBarSweep)
    solK = stpwe_solve_k(sys, p.g * p.c0 * fBarSweep(iw));
    kNorm = solK.k / p.g;                   % 归一化波矢

    % 筛选感兴趣的 k 模式
    keep = isfinite(kNorm) ...
        & real(kNorm) > -0.5 & real(kNorm) < 0.5 ...
        & abs(imag(kNorm)) < 0.12;

    kStore = [kStore; real(kNorm(keep))];
    fStore = [fStore; repmat(fBarSweep(iw), sum(keep), 1)];
    kiStore = [kiStore; imag(kNorm(keep))];
end

% 绘制频率禁带
figure('Color','w');
scatter(kStore, fStore, 8, abs(kiStore), 'filled');
xlabel('Re(k)\Lambda/(2\pi)');  ylabel('\omega\Lambda/(2\pi c)');
colorbar;  title('Im(k) 指示频率禁带');
```

### 5.3 能带折叠

将扩展区频率折叠到第一时间布里渊区 $[-\Omega/2, \Omega/2)$：

```matlab
sol = stpwe_solve_omega(sys, kValue);
foldedOmega = stpwe_fold_frequency(sol.omega, sys.Omega);
% real(foldedOmega) ∈ [-Omega/2, Omega/2)，imag 保持不变
```

这在需要将傅里叶变换的结果映射到标准 Floquet 能带图中非常有用。

### 5.4 静态参考能带（调制深度 = 0）

```matlab
% 计算静态（无时间调制）情况下的能带作为参考基准
epsCoeff0 = @(n) stm_fourier_modulated_slab(0, n, p);   % 仅 m=0 分量
fStatic = stpwe_static_bands(kValues, epsCoeff0, Nspace, p.g, p.c0, 8);

% 绘制：静态能带（实线）+ 折叠带（虚线）
figure;  hold on;
for ib = 1:size(fStatic,1)
    plot(kBar, fStatic(ib,:), 'k-', 'LineWidth', 1.5);                      % 静态带
    plot(kBar, fStatic(ib,:) + p.OmegaBar, 'k--', 'LineWidth', 0.5);        % 第 1 上折叠
    plot(kBar, fStatic(ib,:) - p.OmegaBar, 'k--', 'LineWidth', 0.5);        % 第 1 下折叠
end
xlabel('k\Lambda/(2\pi)');  ylabel('\omega\Lambda/(2\pi c)');
```

### 5.5 跟踪单条能带

当需要计算某一条特定能带的拓扑量时，需要在整个 BZ 范围内"跟踪"它：

```matlab
% 设置 k 网格（注意：不含重复右端点）
kGrid = -p.g/2 + (0:160) * p.g/161;

% 选项参数
opt.omegaWindow = [0.15 0.70] * p.g * p.c0;    % 频率窗口
opt.maxImag = 0.03 * p.g * p.c0;                % 最大允许虚部
opt.minM0Weight = 0.005;                        % 最小 m0 权重

% 从种子频率开始跟踪
omegaSeed = 0.344 * p.g * p.c0;
band = stpwe_track_band(sys, kGrid, omegaSeed, opt);

% 输出
figure;
subplot(2,1,1);
plot(kGrid/p.g, real(band.omega)/(p.g*p.c0), 'b-', 'LineWidth', 1.2);
xlabel('k\Lambda/(2\pi)');  ylabel('Re(\omega)\Lambda/(2\pi c)');
title('Tracked band');

subplot(2,1,2);
plot(kGrid/p.g, band.continuity, 'k-');
xlabel('k\Lambda/(2\pi)');  ylabel('link continuity');
title('Eigenvector continuity (near 1 = good tracking)');
```

> **跟踪原理：** 对第一个 k 点用频率最近邻选择种子模式，之后逐点选择使交叠 $|\langle R_{k_i}|R_{k_{i+1}}\rangle|$ 最大的模式，同时保持频率连续性。

---

## 6. 拓扑不变量计算

### 6.1 Zak 相位（一维 Wilson 环）

适用于一条孤立的能带，计算跨布里渊区的 Berry 相位：

```matlab
% === 完整流程 ===

% 步骤 1：构建 BZ 缝合矩阵（链接 k = -g/2 和 k = +g/2 的等价状态）
sewing = stpwe_bz_sewing_matrix(sys);

% 步骤 2：跟踪目标能带
Nk = 161;
kGrid = -p.g/2 + (0:Nk-1) * p.g/Nk;   % 不含重复右端点！
omegaSeed = 0.344 * p.g * p.c0;
band = stpwe_track_band(sys, kGrid, omegaSeed, opt);

% 步骤 3：计算双正交 Zak 相位
[zakPhase, info] = zak_phase_biorthogonal( ...
    band.R, band.L, sys.Bomega, sewing);

% 显示结果
fprintf('Zak phase = %.6f * pi\n', zakPhase / pi);
fprintf('Minimum Wilson link = %.3e\n', info.minimumLinkMagnitude);

% 重要诊断：如果 minimumLinkMagnitude < 0.1，说明能带接近简并，
% 单带 Zak 相位不可靠，需要使用多带 Wilson 环
if info.minimumLinkMagnitude < 0.1
    warning('Single-band Zak phase may be unreliable due to near-degeneracy.');
end
```

> **关键注意事项：**
> - `kGrid` 必须不包含重复的右端点（否则 BZ 边界链接重复计数）
> - Wilson 链接幅值 |link| 接近 1 时 Zak 相位可靠；接近 0 时说明能带在其他能带附近简并
> - Zak 相位依赖于晶胞原点选择——不同原点给出不同值（但差值应为 $\pi$ 的整数倍）

#### 双正交 Zak 相位的数学意义

对于非厄米系统（$\omega$ 可为复数），左右特征向量不同，必须使用双正交内积：

$$ \gamma_{\text{Zak}} = -\arg\prod_{j=1}^{N_k} \frac{L^\dagger(k_j) \cdot B \cdot R(k_{j+1})}{|L^\dagger(k_j) \cdot B \cdot R(k_{j+1})|} $$

其中 $R(k_{N_k+1}) = \text{sewing} \cdot R(k_1)$。

### 6.2 Chern 数（二维参数空间）

Chern 数描述二维参数空间中的拓扑量子化。这里用 Rice-Mele 模型做验证：

```matlab
% === Rice-Mele 电荷泵示例 ===
% 这是一个紧束缚模型，k 和调制相位构成二维拓扑空间
Nk = 61;  Np = 61;
kGrid = -pi + (0:Nk-1)*2*pi/Nk;
phaseGrid = (0:Np-1)*2*pi/Np;

t0 = 1;  delta0 = 0.6;  mass0 = 1.0;
states = complex(zeros(2, 1, Nk, Np));   % [Nbasis, Nocc, Nk, Np]

for ik = 1:Nk
    for ip = 1:Np
        phase = phaseGrid(ip);
        t1 = t0 + delta0*cos(phase);
        t2 = t0 - delta0*cos(phase);
        mass = mass0*sin(phase);

        H = [mass, t1 + t2*exp(-1i*kGrid(ik));
             conj(t1 + t2*exp(-1i*kGrid(ik))), -mass];
        [V, D] = eig(H, 'vector');
        [~, order] = sort(real(D));
        states(:, 1, ik, ip) = V(:, order(1));   % 占据态
    end
end

% 计算 Chern 数
chern = fhs_chern_number(states);
fprintf('Chern number = %.12f\n', chern);      % 应 ≈ 1
```

> **参数说明：**
>
> `fhs_chern_number(rightStates, leftStates, metric)`
>
> | 参数 | 说明 | 大小 |
> |------|------|------|
> | `rightStates` | 右本征态 | `[Nbasis, Noccupied, Nk, Np]` |
> | `leftStates` | 左本征态（厄米时可不填）| 同 `rightStates` |
> | `metric` | 双正交度规矩阵（厄米时可不填）| `[Nbasis, Nbasis]` |
>
> 返回值 | 含义
> |-------|------
> `chern` | 陈数（应为整数）
> `curvature` | 每个 plaquette 的 Berry 曲率（`[Nk, Np]`）
> `links.k` | k 方向上的非阿贝尔链接矩阵行列式

---

## 7. 场分布图绘制

### 7.1 选择特定模式

在能带图上选定一个 (k, ω) 位置，提取对应的模式：

```matlab
% 选择最接近目标 (k, ω) 的模式
kTarget = 0.175 * p.g;                        % 目标波矢
omegaTarget = 0.10 * p.g * p.c0;              % 目标频率
omegaWindow = [0, 0.9 * p.g * p.c0];          % 可选：频率窗口

mode = stpwe_select_mode(sys, kTarget, omegaTarget, omegaWindow);

% 输出模式信息
fprintf('Selected mode:\n');
fprintf('  k = %.4f\n', mode.k);
fprintf('  ω = %.6f %+.6ei\n', real(mode.omega), imag(mode.omega));
fprintf('  m=0 weight = %.3f\n', mode.m0Weight);
```

### 7.2 重构电场时空分布

```matlab
% 设置空间和时间网格
xField = linspace(0, 3*p.Lambda, 451);       % 空间范围：3 个晶胞
tField = linspace(0, 3*p.T, 361);             % 时间范围：3 个周期

% 重构电场（不包含增长因子 -> 显示稳定传播模式）
field = stpwe_reconstruct_field(mode, xField, tField, false);

% 或者包含增长因子（-> 显示动量禁带内的增益/衰减）
fieldWithGrowth = stpwe_reconstruct_field(mode, xField, tField, true);

% 绘制时空图
figure('Color','w');
imagesc(xField/p.Lambda, tField/p.T, field);
set(gca, 'YDir', 'normal');
caxis([-1 1]);
colormap(stm_redblue(256));
xlabel('x/\Lambda');
ylabel('t/T');
title('E(x,t) normalized');
colorbar;

% 叠加调制边界线
hold on;
for cellNo = 0:2
    xline(cellNo + p.xModStart/p.Lambda, 'k--', 'LineWidth', 0.45);
end
```

> **`stpwe_reconstruct_field` 参数：**
>
> | 参数 | 含义 |
> |------|------|
> | `mode` | `stpwe_select_mode` 或手动构建的模式结构体 |
> | `x` | 空间网格向量 |
> | `t` | 时间网格向量 |
> | `includeGrowth` | `false`= 用 Re(ω) 做载波（稳定传播）；`true`= 用 ω（含 Im(ω) 增长因子） |

---

## 8. TMM 交叉验证

时间转移矩阵方法（TMM）独立于 PWE，是验证频散关系的有力工具。

### 8.1 二元时间晶体的能带（TMM 方法）

```matlab
% 定义材料参数
epsA = 1.0;  epsB = 4.0;          % 介电常数值
muA = 1.0;   muB = 1.0;           % 磁导率
dutyA = 0.5;                       % εA 的占空比
T_period = 1;                      % 周期
durations = [dutyA, 1-dutyA] * T_period;
Omega = 2*pi / T_period;

% 扫描 k 空间
kNorm = linspace(0.02, 1.45, 181);  % k*c0/Ω
kValues = kNorm * Omega;

% TMM 能带计算
tmmBands = temporal_crystal_bands(kValues, [epsA epsB], [muA muB], durations);

% 绘制结果
figure;
subplot(1,2,1);
plot(kNorm, real(tmmBands.omegaF.')/Omega, 'k-', 'LineWidth', 1.5);
xlabel('kc/\Omega');  ylabel('Re(\omega_F)/\Omega');
title('TMM Floquet bands');  grid on;

subplot(1,2,2);
plot(kNorm, imag(tmmBands.omegaF.')/Omega, 'r-', 'LineWidth', 1.5);
xlabel('kc/\Omega');  ylabel('Im(\omega_F)/\Omega');
title('Momentum gaps');  grid on;  yline(0, 'k:');
```

### 8.2 PWE 与 TMM 交叉验证

```matlab
% 用 PWE 计算同一系统（注意：方波需要大 Mtime）
Mtime = 19;
epsCoeff = @(m,n) double(n==0) * temporal_binary_eps_coeff(m, epsA, epsB, dutyA);
sysPTC = stpwe_build_system(epsCoeff, [], 0, Mtime, 1, Omega);

% 对每个 k 匹配 PWE 和 TMM 结果
for ik = 1:numel(kValues)
    sol = stpwe_solve_omega(sysPTC, kValues(ik));
    pweFolded = stpwe_fold_frequency(sol.omega, Omega);
    % ... 然后通过频率最近邻匹配到 TMM 结果
end
```

### 8.3 多层时间结构 TMM

```matlab
% 设计一个四层时间结构
omega0 = 2*pi;                              % 归一化频率
epsInitial = 1;  epsFinal = 1;
epsSlabs = [6.25, 19.36, 4, 2.56];         % 各层的 ε
durations = [0.5, 1.0, 3.5, 0.25];          % 各层的时间长度

% 计算复合转移矩阵
[TM, details] = temporal_multilayer_tmm(omega0, ...
    epsInitial, 1, epsSlabs, ones(size(epsSlabs)), durations, epsFinal, 1);

% 透射/反射频谱扫描
fNorm = linspace(0.5, 1.5, 401);
spectrum = temporal_tmm_spectrum(2*pi*fNorm, ...
    epsInitial, 1, epsSlabs, ones(size(epsSlabs)), durations, epsFinal, 1);

figure('Color','w');
plot(fNorm, abs(spectrum.forward), 'b-', 'LineWidth', 1.5);
hold on;
plot(fNorm, abs(spectrum.backward), 'r--', 'LineWidth', 1.5);
xlabel('f/f_0');  ylabel('Electric-field amplitude');
legend('Forward |E^+_{out}/E^+_{in}|', 'Backward |E^-_{out}/E^+_{in}|');
grid on;  box on;
title('Temporal multilayer transfer function');
```

> **Temporal 多层结构原理（类比空间多层膜）：**
>
> | 空间多层膜 | 时间多层结构 |
> |-----------|-------------|
> | 空间界面（不同 ε 的边界） | 时间界面（不同 ε 切换的时刻） |
> | 传输矩阵法串联各层 | `temporal_multilayer_tmm` 按时间顺序串联 |
> | k 频率守恒，ω 波矢改变 | k 波矢守恒，ω 频率改变 |
> | Fresnel 公式 | Morgenthaler 矩阵 (`temporal_interface_matrix`) |

---

## 9. FDTD 时域仿真验证

FDTD 提供独立的时间域验证，尤其适合观察波包演化和参量过程。

### 9.1 快速启动 FDTD

```matlab
% 网格设置
lambda0 = 1;                         % 参考波长
dx = lambda0 / 32;                   % 空间步长
x = 0:dx:60;                         % 空间网格
dt = 0.80 * dx;                      % 时间步长（CFL 条件）
nSteps = ceil(30 / dt);              % 总时间步数
CFL = 0.8;                           % 保证稳定性

% 初始波包：高斯包络 × 平面波
k0 = 2*pi/lambda0;                   % 中心波矢
x0 = 17;  sigma = 3.5;
profile = @(xq) exp(-((xq-x0)/sigma).^2) .* exp(1i*k0*(xq-x0));
nBefore = 1.5;                       % 初始折射率

cfg.x = x;
cfg.dt = dt;
cfg.nSteps = nSteps;
cfg.epsFun = @(xq,tq) nBefore^2 * ones(size(xq));  % 均匀介质
cfg.muFun = @(xq,tq) ones(size(xq));
cfg.E0 = profile(x);
cfg.boundary = 'sponge';             % 海绵吸收边界
cfg.spongeCells = 100;
cfg.spongeStrength = 0.08;
cfg.recordEvery = 2;

out = fdtd1d_db(cfg);

% 绘制演化
figure('Color','w');
imagesc(out.x, out.t, abs(out.E));
set(gca, 'YDir', 'normal');
xlabel('x');  ylabel('Time');  title('|E(x,t)|');
colorbar;
```

### 9.2 时间界面上的波分裂

当介质在 $t = t_{\text{switch}}$ 时刻发生突变，入射波分裂为前向波和后向波：

```matlab
tSwitch = 12;  tEnd = 30;
nBefore = 1.5;  nAfter = 2.5;

cfg.epsFun = @(xq,tq) (tq < tSwitch) * nBefore^2 + (tq >= tSwitch) * nAfter^2;
% 或使用 temporal_eps 辅助函数

out = fdtd1d_db(cfg);

% 提取时间界面前后的场，做方向分解
[~, beforeId] = min(abs(out.t - (tSwitch - 2*dt)));
[~, afterId]  = min(abs(out.t - (tSwitch + 2*dt)));
Ebefore = out.E(beforeId, :);
Eafter  = out.E(afterId, :);
Hafter  = out.H(afterId, :);

% 方向分解：E± = 0.5*(E ± H/n)
nAfter = sqrt(nAfter^2);  % 注意这里假设 μ=1
Eplus  = 0.5 * (Eafter + Hafter / nAfter);
Eminus = 0.5 * (Eafter - Hafter / nAfter);

% 与解析理论的 Morgenthaler 系数对比
[~, tauExact, rhoExact] = temporal_interface_matrix(nBefore^2, 1, nAfter^2, 1);
tauNum = norm(Eplus) / norm(Ebefore);
rhoNum = norm(Eminus) / norm(Ebefore);
fprintf('|τ|: analytic=%.5f, FDTD=%.5f\n', abs(tauExact), tauNum);
fprintf('|ρ|: analytic=%.5f, FDTD=%.5f\n', abs(rhoExact), rhoNum);
```

> **FDTD 配置参数详解：**
>
> | 字段 | 含义 | 默认值 |
> |------|------|--------|
> | `cfg.x` | 均匀空间网格（至少 5 点）| 必填 |
> | `cfg.dt` | 时间步长 | 必填（建议 `0.8*dx/c0`） |
> | `cfg.nSteps` | 总步数 | 必填 |
> | `cfg.epsFun(x,t)` | 时变介电常数函数 | 必填 |
> | `cfg.muFun(x,t)` | 时变磁导率函数 | `@(x,t) ones(size(x))` |
> | `cfg.boundary` | `'sponge'` 或 `'periodic'` | `'sponge'` |
> | `cfg.spongeCells` | 吸收边界宽度 | `min(80, floor(Nx/8))` |
> | `cfg.spongeStrength` | 吸收强度 | `0.12` |
> | `cfg.E0` | 初始电场 | 全零 |
> | `cfg.Hhalf0` | 初始磁场（H 网格，t=-dt/2）| 全零 |
> | `cfg.sourceD(x,t,step)` | 可选的 D 源项 | 无 |
> | `cfg.recordEvery` | 记录间隔 | 1 |

### 9.3 STFT 频谱分析

分析波包通过时空介质后的频率成分变化：

```matlab
% 在某个探针点提取信号
[~, probeId] = min(abs(out.x - 14));
probe = out.E(:, probeId);
dtRecord = mean(diff(out.t));

% 短时傅里叶变换
[Sprobe, fProbe, tProbe] = stm_stft(probe, dtRecord, 128, 12, 512);

% 绘制频谱图
figure('Color','w');
imagesc(tProbe, fProbe, 20*log10(abs(Sprobe)/max(abs(Sprobe(:))) + 1e-8));
set(gca, 'YDir', 'normal');
xlabel('Time');  ylabel('Frequency');
title('Spectrogram (log scale)');
colormap(parula);  colorbar;
caxis([-60 0]);
```

---

## 10. 完整科研流程示例

下面是一个 "从零开始" 的完整流程：设计一个自定义时空晶体 → 计算能带 → 提取模式 → 计算拓扑 → 绘制场图。

```matlab
%% ========== 自定义时空晶体的完整分析流程 ==========
clear;  close all;
rootDir = stm_init();

%% ===== 第一步：设计参数 =====
% 定义：空间周期中 1/3 的区域以正弦调制，调制频率 0.15
p.Lambda = 1;
p.c0 = 1;
p.g = 2*pi/p.Lambda;

p.epsBack = 2;          % 静态背景
p.epsMod = 6;           % 调制区域基准
p.modDepth = 0.5;       % 调制深度
p.xModStart = 2*p.Lambda/3;
p.xModEnd = p.Lambda;

p.OmegaBar = 0.15;      % 归一化调制频率
p.Omega = 2*pi*p.c0/p.Lambda * p.OmegaBar;
p.T = 2*pi/p.Omega;

fprintf('Lambda=%.2f, T=%.2f, Omega=%.4f\n', p.Lambda, p.T, p.Omega);

%% ===== 第二步：定义材料函数和傅里叶系数 =====

% 实空间介电常数
epsilon_fun = @(x,t) custom_epsilon(x, t, p);

% 解析傅里叶系数（只需要支持 m=0,±1）
epsCoeff = @(m,n) custom_eps_coeff(m, n, p);

%% ===== 第三步：构建 ST-PWE 系统 =====
Nspace = 12;            % 空间截断
Mtime = 1;              % 正弦调制，1 阶足够
sys = stpwe_build_system(epsCoeff, [], Nspace, Mtime, p.g, p.Omega);
fprintf('System matrix size: %d x %d\n', 2*sys.S, 2*sys.S);

%% ===== 第四步：计算能带结构 =====
Nk = 301;
kBar = linspace(-0.5, 0.5, Nk);
kValues = p.g * kBar;
fMax = 0.6;
imagTol = 3e-3;
weightThresh = 0.02;

allK = [];  allF = [];  allW = [];  allFI = [];
for ik = 1:Nk
    sol = stpwe_solve_omega(sys, kValues(ik));
    fBar = sol.omega / (p.g*p.c0);
    keep = isfinite(fBar) & real(fBar)>=0 & real(fBar)<=fMax ...
        & abs(imag(fBar))<=imagTol & sol.m0Weight>=weightThresh;
    allK = [allK; repmat(kBar(ik), sum(keep), 1)];
    allF = [allF; real(fBar(keep))];
    allW = [allW; sol.m0Weight(keep)];
    allFI = [allFI; abs(imag(fBar(keep)))];
end

% 静态参考能带
epsCoeff0 = @(n) epsCoeff(0, n);
fStatic = stpwe_static_bands(kValues, epsCoeff0, Nspace, p.g, p.c0, 6);

% 绘制能带图
figure('Color','w', 'Position', [100 100 800 500]);
hold on;
for ib = 1:size(fStatic,1)
    plot(kBar, fStatic(ib,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
    for s = [-2 -1 1 2]
        plot(kBar, fStatic(ib,:) + s*p.OmegaBar, '--', ...
            'Color', [0.7 0.7 0.7], 'LineWidth', 0.4);
    end
end
scatter(allK, allF, 6, allW, 'filled');
xlim([-0.5 0.5]);  ylim([0 fMax]);
xlabel('k\Lambda/(2\pi)');  ylabel('\omega\Lambda/(2\pi c)');
colormap(flipud(gray));  colorbar;
title('Custom spacetime crystal: Floquet band structure');
grid on;  box on;

%% ===== 第五步：拓扑分析（Zak 相位）=====
NkTopo = 161;
kGrid = -p.g/2 + (0:NkTopo-1)*p.g/NkTopo;
omegaSeed = 0.25 * p.g * p.c0;

opt.omegaWindow = [0.15 0.50] * p.g * p.c0;
opt.maxImag = 0.02 * p.g * p.c0;
opt.minM0Weight = 0.01;

band = stpwe_track_band(sys, kGrid, omegaSeed, opt);
sewing = stpwe_bz_sewing_matrix(sys);
[zak, info] = zak_phase_biorthogonal(band.R, band.L, sys.Bomega, sewing);
fprintf('Zak phase = %.4f * pi\n', zak/pi);
fprintf('Min Wilson link = %.3e\n', info.minimumLinkMagnitude);

%% ===== 第六步：选择模式并绘制场图 =====
kTargetBar = 0.15;
fTargetBar = 0.20;
mode = stpwe_select_mode(sys, p.g*kTargetBar, p.g*p.c0*fTargetBar, [0 0.7*p.g*p.c0]);
fprintf('\nSelected mode:\n');
fprintf('  kBar=%.4f, fBar=%.5f%+.5ei, m0Weight=%.3f\n', ...
    mode.k/p.g, real(mode.omega/(p.g*p.c0)), imag(mode.omega/(p.g*p.c0)), mode.m0Weight);

xField = linspace(0, 3*p.Lambda, 301);
tField = linspace(0, 2*p.T, 201);
field = stpwe_reconstruct_field(mode, xField, tField, false);

figure('Color','w', 'Position', [200 200 600 400]);
imagesc(xField/p.Lambda, tField/p.T, field);
set(gca, 'YDir', 'normal');
caxis([-1 1]);  colormap(stm_redblue(256));
xlabel('x/\Lambda');  ylabel('t/T');
title(sprintf('E(x,t) at k=%.3f, f=%.3f', kTargetBar, fTargetBar));
colorbar;
hold on;
for c = 0:2
    xline(c + p.xModStart/p.Lambda, 'k--', 'LineWidth', 0.5);
end

%% ===== 辅助函数 =====
function epsr = custom_epsilon(x, t, p)
    xCell = mod(x, p.Lambda);
    mask = (xCell >= p.xModStart) & (xCell < p.xModEnd);
    epsr = p.epsBack * ones(size(x));
    epsr(mask) = p.epsMod * (1 + p.modDepth * sin(p.Omega * t));
end

function coeff = custom_eps_coeff(m, n, p)
    staticFrac = (p.xModStart)/p.Lambda;
    modFrac = (p.xModEnd - p.xModStart)/p.Lambda;
    switch m
        case 0
            coeff = (n==0)*(p.epsBack*staticFrac + p.epsMod*modFrac);
        case 1
            coeff = 1i*(p.epsMod*p.modDepth/2)*stm_interval_fourier(n, p.xModStart, p.xModEnd, p.Lambda);
        case -1
            coeff = -1i*(p.epsMod*p.modDepth/2)*stm_interval_fourier(n, p.xModStart, p.xModEnd, p.Lambda);
        otherwise
            coeff = 0;
    end
end
```

---

## 11. 函数索引

### `core/` — ST-PWE 核心计算

| 函数 | 输入要点 | 输出要点 | 依赖 |
|------|---------|---------|------|
| `stpwe_build_system` | `(epsCoeff, muCoeff, Nspace, Mtime, g, Omega)` | `sys` 结构体 | — |
| `stpwe_solve_omega` | `(sys, k)` | 复 `ω`，`R`/`L` 本征向量，`m0Weight` | `stpwe_build_system` |
| `stpwe_solve_k` | `(sys, omega)` | 复 `k`，`R`/`L` 本征向量 | `stpwe_build_system` |
| `stpwe_fold_frequency` | `(omega, Omega)` | 折叠到 $[-\Omega/2,\Omega/2)$ 的频率 | — |
| `stpwe_select_mode` | `(sys, kTarget, omegaTarget, [omegaWindow])` | `mode` 结构体 | `stpwe_solve_omega` |
| `stpwe_reconstruct_field` | `(mode, x, t, [includeGrowth])` | 归一化 $E(x,t)$ 矩阵 | — |
| `stpwe_static_bands` | `(kValues, epsCoeff0, Nspace, g, c0, nBands)` | `fBands` 矩阵 | — |
| `stpwe_sample_fourier_coefficients` | `(materialFun, mOrders, nOrders, Lambda, T, Nx, Nt)` | `table` 结构体 | — |
| `stpwe_lookup_coefficient` | `(table, m, n)` | 傅里叶系数值 | 需要 `table` |
| `stm_preset_modulated_slab` | — | Park-Min 模型参数结构体 `p` | — |
| `stm_permittivity_modulated_slab` | `(x, t, p)` | $\varepsilon(x,t)$ 矩阵 | — |
| `stm_fourier_modulated_slab` | `(m, n, p)` | $\varepsilon_{mn}$ 解析系数 | `stm_interval_fourier` |
| `stm_interval_fourier` | `(n, xa, xb, Lambda)` | 子区间傅里叶积分 | — |
| `stm_redblue` | `(n)` | $n\times 3$ 色图 | — |
| `stm_stft` | `(signal, dt, [window], [hop], [nFFT])` | 频谱，频率轴，时间轴 | — |

### `tmm/` — 时间转移矩阵

| 函数 | 输入要点 | 输出要点 | 依赖 |
|------|---------|---------|------|
| `temporal_interface_matrix` | `(ε_before, μ_before, ε_after, μ_after)` | `MM, τ, ρ` | — |
| `temporal_delay_matrix` | `(ω_init, n_init, n_slab, duration)` | 延迟矩阵 | — |
| `temporal_crystal_monodromy` | `(k, eps_seq, mu_seq, durations)` | 单值矩阵 $U$ | — |
| `temporal_crystal_bands` | `(kValues, eps_seq, mu_seq, durations)` | 能带结构体 | `temporal_crystal_monodromy` |
| `temporal_multilayer_tmm` | `(ω, ε_i, μ_i, ε_slabs, μ_slabs, durations, ε_f, μ_f)` | `TM, details` | `temporal_interface_matrix`, `temporal_delay_matrix` |
| `temporal_tmm_spectrum` | `(ω_range, ε_i, μ_i, ε_slabs, μ_slabs, durations, ε_f, μ_f)` | 频谱结构体 | `temporal_multilayer_tmm` |
| `temporal_binary_eps_coeff` | `(m, εA, εB, dutyA)` | 傅里叶系数值 | — |

### `fdtd/` — 时域有限差分

| 函数 | 输入要点 | 输出要点 | 依赖 |
|------|---------|---------|------|
| `fdtd1d_db` | `cfg` 结构体（见 9.1 节） | `out` 结构体（含 E/H/能量历史） | — |

### `topology/` — 拓扑不变量

| 函数 | 输入要点 | 输出要点 | 依赖 |
|------|---------|---------|------|
| `stpwe_bz_sewing_matrix` | `sys` | 置换矩阵 | `stpwe_build_system` |
| `stpwe_track_band` | `(sys, kGrid, omegaSeed, options)` | `band` 结构体 | `stpwe_solve_omega` |
| `zak_phase_biorthogonal` | `(rightStates, leftStates, metric, sewing)` | Zak 相位，`info` | — |
| `fhs_chern_number` | `(rightStates, [leftStates], [metric])` | Chern 数，曲率，链接 | — |

### 顶层脚本

| 脚本 | 功能 | 依赖 |
|------|------|------|
| `stm_init()` | 设置路径，创建输出目录 | — |
| `stm_run_examples(runHeavy)` | 按顺序运行所有示例 | 全部 |
| `stm_selftest()` | 自检测试 | `core/`, `tmm/`, `fdtd/`, `topology/` |

### 示例脚本 (`examples/`)

| 脚本 | 功能 | 依赖 |
|------|------|------|
| `example_floquet_bands_and_fields(quality)` | 重现 Park-Min 论文 Fig.2 | `core/` |
| `example_complex_gaps()` | 复频率动量带隙和复波数频率带隙 | `core/` |
| `example_ptc_pwe_vs_tmm()` | PWE 与 TMM 交叉验证 | `core/`, `tmm/` |
| `example_tmm_multilayer()` | 四类时间多层结构 | `tmm/` |
| `example_fdtd_interface()` | 单个时间界面的波分裂 | `fdtd/`, `tmm/` |
| `example_fdtd_wavepacket()` | 时空晶体中的波包演化 | `fdtd/`, `core/` |
| `example_zak_phase()` | 双正交 Zak 相位计算 | `core/`, `topology/` |
| `example_chern_pump()` | Rice-Mele 泵 Chern 数 | `topology/` |

---

> **引用和进一步阅读：**
>
> - Park & Min, "Space-time photonic crystals and topological phenomena," *Physical Review B*, 2024. — ST-PWE 方法的理论基础
> - Ramaccia et al., "Temporal multilayer structures," *Applied Physics Letters* 118, 101901 (2021). — TMM 方法参考
> - Fukui, Hatsugai & Suzuki, "Chern numbers in discretized Brillouin zones," *J. Phys. Soc. Jpn.* 74, 1674 (2005). — 离散陈数算法
> - Morgenthaler, "Velocity modulation of electromagnetic waves," *IRE Trans. Microwave Theory Tech.* 6, 167 (1958). — 时间界面边界条件
