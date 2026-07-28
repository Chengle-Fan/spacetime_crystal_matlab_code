# STM 工具包用户手册：时空光子晶体计算

> **STM (Space-Time Media) Toolbox** — 一维时空调制介质的光子能带结构、拓扑不变量和时域仿真
>
> 本文档合并原 `MANUAL.md` 和 `README_CN.md`（§3-10），去重后形成完整的中文技术手册。
> 返回 [README.md](../README.md) · 参见 [conventions.md](conventions.md)（物理约定）· [paper-map.md](paper-map.md)（文献覆盖）

---

## 目录

1. [统一场与 Fourier 约定](#1-统一场与-fourier-约定)
2. [ST-PWE 的两种本征问题](#2-st-pwe-的两种本征问题)
3. [参数设计：定义时空晶体](#3-参数设计定义你的时空晶体)
4. [构建系统矩阵](#4-构建系统矩阵)
5. [能带结构计算](#5-能带结构计算)
6. [拓扑不变量计算](#6-拓扑不变量计算)
7. [场分布图绘制](#7-场分布图绘制)
8. [TMM 交叉验证](#8-tmm-交叉验证)
9. [FDTD 时域仿真验证](#9-fdtd-时域仿真验证)
10. [完整科研流程示例](#10-完整科研流程示例)
11. [科研计算的收敛顺序](#11-科研计算的收敛顺序)
12. [适用范围与尚未包含的物理](#12-适用范围与尚未包含的物理)
13. [核心参考文献](#13-核心参考文献)
    附录：[函数索引](#附录-函数索引)

---

## 1. 统一场与 Fourier 约定

### 1.1 场展开

ST-PWE 使用

$$
\begin{pmatrix}E\\H\end{pmatrix}
=e^{i(kx-\omega t)}
\sum_{m,n}\widetilde{\Psi}_{m,n}
e^{i(ng x-m\Omega t)} .
$$

材料参数写为

$$
p(x,t)=\sum_{m,n}p_{m,n}e^{i(ng x-m\Omega t)} .
$$

因此 Fourier 系数是

$$
p_{m,n}=\frac{1}{\Lambda T}
\int_0^\Lambda\int_0^T
p(x,t)e^{-ingx+im\Omega t}\,dt\,dx .
$$

### 1.2 增长符号

在此约定下：

- $\operatorname{Im}\omega>0$ 表示 $e^{-i\omega t}$ 随时间**增长**；
- $\operatorname{Im}\omega<0$ 表示衰减；
- $\operatorname{Im}k>0$ 或 $<0$ 的物理解读取决于空间传播方向和所选边界条件。

> **注意：** Ramaccia 时间多层示例沿用其论文的 $e^{+i\omega t}$ 符号。比较幅度不受影响；比较复相位时必须先统一约定（见 §8）。

### 1.3 归一化

Park-Min 图 2 的默认归一化为

$$
\Lambda=c=1,\quad
\bar{k}=\frac{k\Lambda}{2\pi},\quad
\bar{f}=\frac{\omega\Lambda}{2\pi c},\quad
\bar{\Omega}=0.2 .
$$

### 1.4 三种数值方法概览

| 方法             | 目录      | 适用场景                                                  |
| ---------------- | --------- | --------------------------------------------------------- |
| **ST-PWE** | `core/` | Fourier 域求解，快速得到全频带结构，支持复数频率/波矢分析 |
| **TMM**    | `tmm/`  | 时间转移矩阵，精确求解多层时间结构，验证复数能带          |
| **FDTD**   | `fdtd/` | 时域仿真，直观展示波包演化，验证理论预测                  |

---

## 2. ST-PWE 的两种本征问题

### 2.1 固定 k → 复 ω（动量禁带分析）

`stpwe_solve_omega` 实现固定 $k$ 的广义本征问题

$$
A(k)\Phi=\omega B\Phi ,
$$

适合寻找动量带隙中的复频率和参量增长/衰减模。

### 2.2 固定 ω → 复 k（频率禁带分析）

`stpwe_solve_k` 实现固定 $\omega$ 的普通本征问题

$$
\mathcal{K}(\omega)\Phi=k\Phi ,
$$

适合寻找频率带隙中的复 Bloch 波数。

运行以下 demo 可同时看到两类复谱：

```matlab
demo02_complex_frequency_and_momentum_gaps
```

二者不能简单互换：一个描述时间不稳定性，另一个描述空间倏逝。

### 2.3 复现论文图 2

```matlab
% 快速截断
demo01_reproduce_fig2_stpwe('quick')

% 论文质量截断
demo01_reproduce_fig2_stpwe('paper')
```

`paper` 使用 20 阶正负空间谐波和一阶正负时间谐波，即 $41\times3=123$ 个时空 Fourier 状态、246 阶 Maxwell 广义本征矩阵。

原始 `time.m` 中的着色量实际是 $m=0$ Fourier 扇区占比，并不严格等于论文所述的"投影到所有静态本征模后的能量权重"。本代码把它明确标记为 `m=0 sector participation`，避免把近似量误称为精确静态带权重。能带本征值、复带隙和场重构不依赖这一着色近似。

### 2.4 改成任意时空材料

若可以解析写出系数，只需提供函数：

```matlab
epsCoeff = @(m,n) my_epsilon_coefficient(m,n);
muCoeff  = @(m,n) double(m==0 && n==0);
sys = stpwe_build_system(epsCoeff,muCoeff,12,3,g,Omega);
solution = stpwe_solve_omega(sys,k);
```

若材料分布较复杂，可直接数值取系数：

```matlab
Lambda = 1; T = 5;
epsFun = @(x,t) 3 + 0.4*cos(2*pi*x/Lambda-2*pi*t/T);

table = stpwe_sample_fourier_coefficients( ...
    epsFun,-6:6,-12:12,Lambda,T,512,512);
epsCoeff = @(m,n) stpwe_lookup_coefficient(table,m,n);

sys = stpwe_build_system(epsCoeff,[],12,6,2*pi/Lambda,2*pi/T);
```

采样区间不包含重复终点，避免 FFT/积分时重复一个周期边界。

---

## 3. 参数设计：定义你的时空晶体

这是整个流程的起点。你需要做三件事：

1. 定义介电常数 $\varepsilon(x,t)$ 的空间分布和时空调制方式
2. （可选）定义磁导率 $\mu(x,t)$
3. 设定截断阶数 $N_{\text{space}}$ 和 $M_{\text{time}}$

### 3.1 方式 A：提供解析 Fourier 系数（推荐，精度最高）

直接提供 $\varepsilon_{mn}$ 的解析表达式。以 Park-Min Fig.2 为例，介电常数在一个周期内的构型：

| $x$ 范围                        | $\varepsilon(x,t)$                                                  |
| --------------------------------- | --------------------------------------------------------------------- |
| $[0, \frac{3}{4}\Lambda)$       | $\varepsilon_1 = 2$（静态）                                         |
| $[\frac{3}{4}\Lambda, \Lambda)$ | $\varepsilon_c[1 + \text{modDepth}\cdot\sin(\Omega t)]$（时空调制） |

#### 步骤 1：加载预设参数

```matlab
p = stm_preset_modulated_slab();
% p 包含: Lambda, c0, g, eps1, epsc, modDepth, xModStart, xModEnd, OmegaBar, Omega, T
```

> **关键参数设置建议：**
>
> - `modDepth` ∈ (0, 1)：调制深度越大，动量禁带越宽
> - `OmegaBar`：调制频率与空间布里渊区边界的频率之比，决定禁带位置
> - `Nspace` 和 `Mtime`：截断阶数，决定计算精度。正弦调制 `Mtime=1` 已精确；方波调制需要更大的 `Mtime`

#### 步骤 2：编写 $\varepsilon(x,t)$ 的实空间函数

```matlab
function epsr = stm_permittivity_modulated_slab(x, t, p)
    xCell = mod(x, p.Lambda);
    mask = (xCell >= p.xModStart) & (xCell < p.xModEnd);
    epsr = p.eps1 * ones(size(x));
    epsr(mask) = p.epsc * (1 + p.modDepth * sin(p.Omega * t));
end
```

#### 步骤 3：提供解析 Fourier 系数

```matlab
function coeff = stm_fourier_modulated_slab(m, n, p)
    staticPart = stm_interval_fourier(n, 0, p.xModStart, p.Lambda);
    modPart = stm_interval_fourier(n, p.xModStart, p.xModEnd, p.Lambda);
    switch m
        case 0
            coeff = p.eps1*staticPart + p.epsc*modPart;
        case 1
            coeff = 1i*(p.epsc*p.modDepth/2)*modPart;
        case -1
            coeff = -1i*(p.epsc*p.modDepth/2)*modPart;
        otherwise
            coeff = 0;
    end
end
```

#### 辅助工具：空间区间 Fourier 积分

```matlab
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

对于没有解析 Fourier 系数的任意 $\varepsilon(x,t)$ 分布：

```matlab
materialFun = @(X, T) stm_permittivity_modulated_slab(X, T, p);
mOrders = -1:1;  nOrders = -3:3;
table = stpwe_sample_fourier_coefficients(materialFun, ...
    mOrders, nOrders, p.Lambda, p.T, 512, 256);
epsCoeff = @(m, n) stpwe_lookup_coefficient(table, m, n);
```

> | 参数   | 含义                     | 建议值                         |
> | ------ | ------------------------ | ------------------------------ |
> | `Nx` | 每个空间周期内的采样点数 | `max(128, 8*numel(nOrders))` |
> | `Nt` | 每个时间周期内的采样点数 | `max(128, 8*numel(mOrders))` |
>
> 采样使用均匀端点外网格（避免 Gibbs 效应）。

### 3.3 方式 C：二元（方波型）时间晶体

```matlab
epsA = 1.0;  epsB = 4.0;  dutyA = 0.5;
epsCoeff = @(m, n) double(n == 0) * temporal_binary_eps_coeff(m, epsA, epsB, dutyA);
```

> **三种方式对比：**
>
> | 方式        | 优点               | 缺点               | 适用场景                 |
> | ----------- | ------------------ | ------------------ | ------------------------ |
> | A. 解析系数 | 精度最高，计算最快 | 需要推导解析式     | 标准波形（正弦、方波等） |
> | B. 数值采样 | 通用性好           | 需要较高采样率     | 任意复杂波形             |
> | C. 二元解析 | 对间断波形精确     | 收敛慢，需大 Mtime | 方波时间晶体             |

---

## 4. 构建系统矩阵

将 Fourier 系数组装成特征值问题所需的大矩阵。

### 4.1 调用方式

```matlab
sys = stpwe_build_system(epsCoeff, muCoeff, Nspace, Mtime, g, Omega)

% 例 1：非磁性时空晶体
epsCoeff = @(m,n) stm_fourier_modulated_slab(m, n, p);
sys = stpwe_build_system(epsCoeff, [], 10, 1, p.g, p.Omega);

% 例 2：磁性材料
muCoeff = @(m,n) double(m==0 && n==0);
sys = stpwe_build_system(epsCoeff, muCoeff, 12, 2, p.g, p.Omega);
```

### 4.2 参数说明

| 参数              | 含义                                                        | 典型值                                    |
| ----------------- | ----------------------------------------------------------- | ----------------------------------------- |
| `epsCoeff(m,n)` | 返回介电常数$mn$ 阶 Fourier 系数的函数句柄                | —                                        |
| `muCoeff(m,n)`  | 返回磁导率$mn$ 阶 Fourier 系数的函数句柄                  | `[]` 或 `@(m,n) double(m==0 && n==0)` |
| `Nspace`        | 空间截断阶数$n \in [-N_{\text{space}}, N_{\text{space}}]$ | 静态区 5-15，时空区 8-20                  |
| `Mtime`         | 时间截断阶数$m \in [-M_{\text{time}}, M_{\text{time}}]$   | 正弦 1-2，方波 15-30                      |
| `g`             | 空间倒格矢$2\pi/\Lambda$                                  | —                                        |
| `Omega`         | 时间调制角频率$2\pi/T$                                    | —                                        |

> **截断阶数选择技巧：**
>
> - **正弦调制**：$M_{\text{time}}=1$ 已精确（仅含 ±1 阶 Fourier 系数）
> - **方波调制**：$M_{\text{time}}$ 需要 15-30 以收敛
> - $N_{\text{space}}$ 取 8-12 通常已足够，论文级精度取 20
> - 总矩阵大小 $2S \times 2S$，$S = (2N_{\text{space}}+1)(2M_{\text{time}}+1)$

### 4.3 sys 结构体字段

| 字段                          | 大小             | 含义                                                           |
| ----------------------------- | ---------------- | -------------------------------------------------------------- |
| `sys.Nspace`, `sys.Mtime` | 标量             | 截断阶数                                                       |
| `sys.nList`, `sys.mList`  | $S \times 1$   | 每个谐波对应的$(n,m)$ 索引                                   |
| `sys.S`                     | 标量             | 谐波总数                                                       |
| `sys.g`, `sys.Omega`      | 标量             | 倒格矢和调制频率                                               |
| `sys.Ceps`, `sys.Cmu`     | $S \times S$   | 介电/磁导率卷积矩阵                                            |
| `sys.G`                     | $S \times S$   | 空间波数矩阵$\text{diag}(n·g)$                              |
| `sys.W`                     | $S \times S$   | 时间频率矩阵$\text{diag}(m·\Omega)$                         |
| `sys.I`, `sys.Z`          | $S \times S$   | 单位阵和零矩阵                                                 |
| `sys.Bomega`                | $2S \times 2S$ | $\begin{bmatrix}0 & C_\mu \\ C_\varepsilon & 0\end{bmatrix}$ |

---

## 5. 能带结构计算

### 5.1 固定 k → 复 ω（动量禁带分析）

```matlab
kBar = linspace(-0.5, 0.5, 201);
kValues = p.g * kBar;

allOmega = [];  allK = [];  allWeight = [];
for ik = 1:numel(kValues)
    sol = stpwe_solve_omega(sys, kValues(ik));
    fBar = sol.omega / (p.g * p.c0);
    keep = isfinite(fBar) & real(fBar) >= 0 & real(fBar) <= fMax ...
        & abs(imag(fBar)) <= imagTolerance & sol.m0Weight >= weightThreshold;
    allK = [allK; repmat(kBar(ik), sum(keep), 1)];
    allOmega = [allOmega; real(fBar(keep))];
    allWeight = [allWeight; sol.m0Weight(keep)];
end

figure('Color','w');
scatter(allK, allOmega, 8, allWeight, 'filled');
xlabel('k\Lambda/(2\pi)');  ylabel('\omega\Lambda/(2\pi c)');
colormap(flipud(gray));  colorbar;
```

> **`sol` 结构体字段：**
>
> | 字段                 | 含义                                     |
> | -------------------- | ---------------------------------------- |
> | `sol.omega`        | 所有本征值（复频率）                     |
> | `sol.R`, `sol.L` | 右/左特征向量矩阵，$L^\dagger B R = I$ |
> | `sol.m0Weight`     | 各模式的$m=0$ 分量占比                 |

#### 筛选参数建议

| 参数                | 作用           | 默认值                     |
| ------------------- | -------------- | -------------------------- |
| `fMax`            | 显示的最大频率 | 取决于研究范围             |
| `imagTolerance`   | 允许的         | Im(ω)                     |
| `weightThreshold` | 最小 m0 权重   | `0.035` 用于筛选物理模式 |

### 5.2 固定 ω → 复 k（频率禁带分析）

```matlab
fBarSweep = linspace(0.38, 0.43, 121);
kStore = [];  fStore = [];  kiStore = [];
for iw = 1:numel(fBarSweep)
    solK = stpwe_solve_k(sys, p.g * p.c0 * fBarSweep(iw));
    kNorm = solK.k / p.g;
    keep = isfinite(kNorm) & real(kNorm) > -0.5 & real(kNorm) < 0.5 ...
        & abs(imag(kNorm)) < 0.12;
    kStore = [kStore; real(kNorm(keep))];
    fStore = [fStore; repmat(fBarSweep(iw), sum(keep), 1)];
    kiStore = [kiStore; imag(kNorm(keep))];
end

figure('Color','w');
scatter(kStore, fStore, 8, abs(kiStore), 'filled');
xlabel('Re(k)\Lambda/(2\pi)');  ylabel('\omega\Lambda/(2\pi c)');
colorbar;  title('Im(k) 指示频率禁带');
```

### 5.3 能带折叠

```matlab
sol = stpwe_solve_omega(sys, kValue);
foldedOmega = stpwe_fold_frequency(sol.omega, sys.Omega);
% real(foldedOmega) ∈ [-Omega/2, Omega/2)，imag 保持不变
```

### 5.4 静态参考能带

```matlab
epsCoeff0 = @(n) stm_fourier_modulated_slab(0, n, p);
fStatic = stpwe_static_bands(kValues, epsCoeff0, Nspace, p.g, p.c0, 8);

figure;  hold on;
for ib = 1:size(fStatic,1)
    plot(kBar, fStatic(ib,:), 'k-', 'LineWidth', 1.5);
    plot(kBar, fStatic(ib,:) + p.OmegaBar, 'k--', 'LineWidth', 0.5);
    plot(kBar, fStatic(ib,:) - p.OmegaBar, 'k--', 'LineWidth', 0.5);
end
xlabel('k\Lambda/(2\pi)');  ylabel('\omega\Lambda/(2\pi c)');
```

### 5.5 跟踪单条能带

```matlab
kGrid = -p.g/2 + (0:160) * p.g/161;

opt.omegaWindow = [0.15 0.70] * p.g * p.c0;
opt.maxImag = 0.03 * p.g * p.c0;
opt.minM0Weight = 0.005;

omegaSeed = 0.344 * p.g * p.c0;
band = stpwe_track_band(sys, kGrid, omegaSeed, opt);

figure;
subplot(2,1,1);
plot(kGrid/p.g, real(band.omega)/(p.g*p.c0), 'b-', 'LineWidth', 1.2);
xlabel('k\Lambda/(2\pi)');  ylabel('Re(\omega)\Lambda/(2\pi c)');

subplot(2,1,2);
plot(kGrid/p.g, band.continuity, 'k-');
xlabel('k\Lambda/(2\pi)');  ylabel('link continuity');
```

> **跟踪原理：** 对第一个 k 点用频率最近邻选择种子模式，之后逐点选择使交叠 $|\langle R_{k_i}|R_{k_{i+1}}\rangle|$ 最大的模式，同时保持频率连续性。v2 使用前后双正交重叠的几何平均，不再用普通右—右欧氏重叠。

---

## 6. 拓扑不变量计算

### 6.1 Zak 相位（一维 Wilson 环）

适用于一条孤立的能带，计算跨布里渊区的 Berry 相位：

```matlab
% 步骤 1：BZ 缝合矩阵
sewing = stpwe_bz_sewing_matrix(sys);

% 步骤 2：跟踪目标能带
Nk = 161;
kGrid = -p.g/2 + (0:Nk-1) * p.g/Nk;   % 不含重复右端点
omegaSeed = 0.344 * p.g * p.c0;
band = stpwe_track_band(sys, kGrid, omegaSeed, opt);

% 步骤 3：计算双正交 Zak 相位
[zakPhase, info] = zak_phase_biorthogonal(band.R, band.L, sys.Bomega, sewing);
fprintf('Zak phase = %.6f * pi\n', zakPhase / pi);
fprintf('Minimum Wilson link = %.3e\n', info.minimumLinkMagnitude);

if info.minimumLinkMagnitude < 0.1
    warning('Single-band Zak phase may be unreliable due to near-degeneracy.');
end
```

> **关键注意事项：**
>
> - `kGrid` 不包含重复右端点
> - Wilson 链接幅值接近 1 时可靠；接近 0 时说明能带简并
> - Zak 相位依赖于晶胞原点选择

#### 双正交 Zak 相位的数学意义

对于非 Hermitian 系统（$\omega$ 可为复数），左右特征向量不同，必须使用双正交内积：

$$
\gamma_{\text{Zak}} = -\arg\prod_{j=1}^{N_k} \frac{L^\dagger(k_j) \cdot B \cdot R(k_{j+1})}{|L^\dagger(k_j) \cdot B \cdot R(k_{j+1})|}
$$

其中 $R(k_{N_k+1}) = \text{sewing} \cdot R(k_1)$。

#### 空间 vs 时间 Zak 相

`demos/demo07` 的闭合路径是空间 Brillouin 区的 $k$，因此它是**空间 Zak 相**，沿空间 Bloch 波数积分，适用于空间 Bloch 带，依赖空间单胞原点和 BZ 端点 sewing。

Lustig 等人 PTC 文献中的**时间 Zak 相**以准频率为时间倒空间坐标，需固定准频率求 $k$ 并重构一周期时间 Bloch 模。当前例程不能互换使用。参见 [roadmap.md](roadmap.md) 里程碑 1。

### 6.2 Chern 数（二维参数空间）

Chern 数描述二维参数空间中的拓扑量子化。这里用 Rice-Mele 模型做算法基准：

```matlab
Nk = 61;  Np = 61;
kGrid = -pi + (0:Nk-1)*2*pi/Nk;
phaseGrid = (0:Np-1)*2*pi/Np;

t0 = 1;  delta0 = 0.6;  mass0 = 1.0;
states = complex(zeros(2, 1, Nk, Np));

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
        states(:, 1, ik, ip) = V(:, order(1));
    end
end

chern = fhs_chern_number(states);
fprintf('Chern number = %.12f\n', chern);      % 应 ≈ 1
```

> `fhs_chern_number(rightStates, leftStates, metric)`
>
> | 参数            | 说明                                 | 大小                            |
> | --------------- | ------------------------------------ | ------------------------------- |
> | `rightStates` | 右本征态                             | `[Nbasis, Noccupied, Nk, Np]` |
> | `leftStates`  | 左本征态（Hermitian 时可省略）       | 同上                            |
> | `metric`      | 双正交度规矩阵（Hermitian 时可省略） | `[Nbasis, Nbasis]`            |

`demos/demo08` 的 $C=1$ 来自 Rice–Mele 泵，是对通用 FHS 实现的独立测试。除非把目标 Maxwell/电路模型的本征态真正送入同一算法，否则这个整数不代表目标时空介质。

---

## 7. 场分布图绘制

### 7.1 选择特定模式

```matlab
kTarget = 0.175 * p.g;
omegaTarget = 0.10 * p.g * p.c0;
omegaWindow = [0, 0.9 * p.g * p.c0];

mode = stpwe_select_mode(sys, kTarget, omegaTarget, omegaWindow);
fprintf('k = %.4f, ω = %.6f%+.6ei, m0Weight = %.3f\n', ...
    mode.k, real(mode.omega), imag(mode.omega), mode.m0Weight);
```

### 7.2 重构电场时空分布

```matlab
xField = linspace(0, 3*p.Lambda, 451);
tField = linspace(0, 3*p.T, 361);

% 不包含增长因子
field = stpwe_reconstruct_field(mode, xField, tField, false);
% 包含增长因子
fieldWithGrowth = stpwe_reconstruct_field(mode, xField, tField, true);

figure('Color','w');
imagesc(xField/p.Lambda, tField/p.T, field);
set(gca, 'YDir', 'normal');
caxis([-1 1]);  colormap(stm_redblue(256));
xlabel('x/\Lambda');  ylabel('t/T');
colorbar;

% 叠加调制边界线
hold on;
for cellNo = 0:2
    xline(cellNo + p.xModStart/p.Lambda, 'k--', 'LineWidth', 0.45);
end
```

> **`stpwe_reconstruct_field(mode, x, t, includeGrowth)`:**
>
> - `includeGrowth=false`：用 $\operatorname{Re}(\omega)$ 做载波（稳定传播）
> - `includeGrowth=true`：用完整 $\omega$（含 $\operatorname{Im}(\omega)$ 增长因子）

---

## 8. TMM 交叉验证

时间转移矩阵方法（TMM）独立于 PWE，是验证频散关系的有力工具。Ramaccia 模型假定每个时间层在空间上均匀，并在整个空间同时切换。时间界面保持 $D$ 与 $B$ 连续。代码采用列向量，因此按实际时间顺序写成：

$$
\mathrm{TM}= \mathrm{MM}_{M+1}\mathrm{DM}_{M}\cdots \mathrm{MM}_{2}\mathrm{DM}_{1}\mathrm{MM}_{1}.
$$

### 8.1 二元时间晶体的能带（TMM 方法）

```matlab
epsA = 1.0;  epsB = 4.0;  muA = 1.0;  muB = 1.0;
dutyA = 0.5;  T_period = 1;
durations = [dutyA, 1-dutyA] * T_period;
Omega = 2*pi / T_period;

kNorm = linspace(0.02, 1.45, 181);
kValues = kNorm * Omega;

tmmBands = temporal_crystal_bands(kValues, [epsA epsB], [muA muB], durations);

figure;
subplot(1,2,1);
plot(kNorm, real(tmmBands.omegaF.')/Omega, 'k-', 'LineWidth', 1.5);
xlabel('kc/\Omega');  ylabel('Re(\omega_F)/\Omega');
grid on;

subplot(1,2,2);
plot(kNorm, imag(tmmBands.omegaF.')/Omega, 'r-', 'LineWidth', 1.5);
xlabel('kc/\Omega');  ylabel('Im(\omega_F)/\Omega');
grid on;  yline(0, 'k:');
```

对周期时间晶体，更稳健的方法是直接推进连续的状态 $[D,B]^T$：

$$
\frac{d}{dt}\begin{pmatrix}D\\B\end{pmatrix}
=-ik
\begin{pmatrix}
0&1/\mu(t)\\
1/\epsilon(t)&0
\end{pmatrix}
\begin{pmatrix}D\\B\end{pmatrix}.
$$

若单周期矩阵的本征值为 $\lambda$，则 $\lambda=e^{-i\omega_FT}$，$\omega_F=\frac{i}{T}\log\lambda$。

### 8.2 PWE 与 TMM 交叉验证

```matlab
Mtime = 19;
epsCoeff = @(m,n) double(n==0) * temporal_binary_eps_coeff(m, epsA, epsB, dutyA);
sysPTC = stpwe_build_system(epsCoeff, [], 0, Mtime, 1, Omega);

for ik = 1:numel(kValues)
    sol = stpwe_solve_omega(sysPTC, kValues(ik));
    pweFolded = stpwe_fold_frequency(sol.omega, Omega);
    % 通过频率最近邻匹配到 TMM 结果
end
```

`demos/demo03_ptc_bands_pwe_vs_tmm` 会把 2×2 精确单周期结果与时间 Fourier PWE 结果画在同一张图上，是修改代码后最重要的交叉验证之一。

### 8.3 多层时间结构 TMM

```matlab
omega0 = 2*pi;
epsInitial = 1;  epsFinal = 1;
epsSlabs = [6.25, 19.36, 4, 2.56];
durations = [0.5, 1.0, 3.5, 0.25];

[TM, details] = temporal_multilayer_tmm(omega0, ...
    epsInitial, 1, epsSlabs, ones(size(epsSlabs)), durations, epsFinal, 1);

fNorm = linspace(0.5, 1.5, 401);
spectrum = temporal_tmm_spectrum(2*pi*fNorm, ...
    epsInitial, 1, epsSlabs, ones(size(epsSlabs)), durations, epsFinal, 1);

figure('Color','w');
plot(fNorm, abs(spectrum.forward), 'b-', 'LineWidth', 1.5);
hold on;
plot(fNorm, abs(spectrum.backward), 'r--', 'LineWidth', 1.5);
xlabel('f/f_0');  ylabel('Electric-field amplitude');
legend('Forward', 'Backward');
grid on;
```

> **Temporal 多层结构原理（类比空间多层膜）：**
>
> | 空间多层膜                 | 时间多层结构                                      |
> | -------------------------- | ------------------------------------------------- |
> | 空间界面（不同 ε 的边界） | 时间界面（不同 ε 切换的时刻）                    |
> | 传输矩阵法串联各层         | `temporal_multilayer_tmm` 按时间顺序串联        |
> | k 频率守恒，ω 波矢改变    | k 波矢守恒，ω 频率改变                           |
> | Fresnel 公式               | Morgenthaler 矩阵 (`temporal_interface_matrix`) |

---

## 9. FDTD 时域仿真验证

FDTD 提供独立的时间域验证，尤其适合观察波包演化和参量过程。

### 9.1 基本原理

`fdtd1d_db` 更新的是

$$
B^{n+1/2}=B^{n-1/2}-\Delta t\,\partial_xE^n,
\qquad
D^{n+1}=D^n-\Delta t\,\partial_xH^{n+1/2},
$$

然后在对应时间用 $E=D/\epsilon(x,t)$，$H=B/\mu(x,t)$ 恢复场。这样在 $\epsilon$ 或 $\mu$ 突变时，$D$ 与 $B$ 不会被错误地重置。

### 9.2 快速启动

```matlab
lambda0 = 1;  dx = lambda0 / 32;
x = 0:dx:60;  dt = 0.80 * dx;
nSteps = ceil(30 / dt);

k0 = 2*pi/lambda0;  x0 = 17;  sigma = 3.5;
profile = @(xq) exp(-((xq-x0)/sigma).^2) .* exp(1i*k0*(xq-x0));

cfg.x = x;  cfg.dt = dt;  cfg.nSteps = nSteps;
cfg.epsFun = @(xq,tq) 2.25 * ones(size(xq));
cfg.muFun = @(xq,tq) ones(size(xq));
cfg.E0 = profile(x);
cfg.boundary = 'sponge';
cfg.spongeCells = 100;  cfg.spongeStrength = 0.08;
cfg.recordEvery = 2;

out = fdtd1d_db(cfg);

figure('Color','w');
imagesc(out.x, out.t, abs(out.E));
set(gca, 'YDir', 'normal');
xlabel('x');  ylabel('Time');  title('|E(x,t)|');
colorbar;
```

### 9.3 时间界面上的波分裂

```matlab
tSwitch = 12;  nBefore = 1.5;  nAfter = 2.5;
cfg.epsFun = @(xq,tq) (tq < tSwitch) * nBefore^2 + (tq >= tSwitch) * nAfter^2;
out = fdtd1d_db(cfg);

[~, beforeId] = min(abs(out.t - (tSwitch - 2*dt)));
[~, afterId]  = min(abs(out.t - (tSwitch + 2*dt)));
Ebefore = out.E(beforeId, :);
Eafter  = out.E(afterId, :);
Hafter  = out.H(afterId, :);

nAfter = sqrt(nAfter^2);
Eplus  = 0.5 * (Eafter + Hafter / nAfter);
Eminus = 0.5 * (Eafter - Hafter / nAfter);

[~, tauExact, rhoExact] = temporal_interface_matrix(nBefore^2, 1, nAfter^2, 1);
fprintf('|τ|: analytic=%.5f, FDTD=%.5f\n', abs(tauExact), norm(Eplus)/norm(Ebefore));
fprintf('|ρ|: analytic=%.5f, FDTD=%.5f\n', abs(rhoExact), norm(Eminus)/norm(Ebefore));
```

> **FDTD 配置参数详解：**
>
> | 字段                      | 含义                           | 默认值                    |
> | ------------------------- | ------------------------------ | ------------------------- |
> | `cfg.x`                 | 均匀空间网格（至少 5 点/波长） | 必填                      |
> | `cfg.dt`                | 时间步长                       | 必填（建议`0.8*dx/c0`） |
> | `cfg.nSteps`            | 总步数                         | 必填                      |
> | `cfg.epsFun(x,t)`       | 时变介电常数函数               | 必填                      |
> | `cfg.muFun(x,t)`        | 时变磁导率函数                 | `@(x,t) ones(size(x))`  |
> | `cfg.boundary`          | `'sponge'` 或 `'periodic'` | `'sponge'`              |
> | `cfg.spongeCells`       | 吸收边界宽度                   | `min(80, floor(Nx/8))`  |
> | `cfg.spongeStrength`    | 吸收强度                       | `0.12`                  |
> | `cfg.E0`                | 初始电场                       | 全零                      |
> | `cfg.Hhalf0`            | 初始磁场（H 网格，t=-dt/2）    | 全零                      |
> | `cfg.sourceD(x,t,step)` | 可选的 D 源项                  | 无                        |
> | `cfg.recordEvery`       | 记录间隔                       | 1                         |

### 9.4 FDTD 稳定性

时间步必须由整个仿真期间的最大相速度决定。归一化一维非色散介质中，保守条件是

$$
\Delta t < \frac{\Delta x}{\max_{x,t}\left[1/\sqrt{\epsilon_r(x,t)\mu_r(x,t)}\right]} .
$$

若允许 $\epsilon_r<1$、接近零、为负值、具有材料色散或有源增益，当前简单更新式不再充分，需要 ADE、卷积或专门的色散 FDTD。

v2 新增：正介质的采样 CFL 审计、`probeIndices` 探针、`storeFields=false` 的长时低内存模式。

### 9.5 STFT 频谱分析

```matlab
[~, probeId] = min(abs(out.x - 14));
probe = out.E(:, probeId);
dtRecord = mean(diff(out.t));

[Sprobe, fProbe, tProbe] = stm_stft(probe, dtRecord, 128, 12, 512);

figure('Color','w');
imagesc(tProbe, fProbe, 20*log10(abs(Sprobe)/max(abs(Sprobe(:))) + 1e-8));
set(gca, 'YDir', 'normal');  xlabel('Time');  ylabel('Frequency');
colormap(parula);  colorbar;  caxis([-60 0]);
```

---

## 10. 完整科研流程示例

下面是一个"从零开始"的完整流程：设计一个自定义时空晶体 → 计算能带 → 提取模式 → 计算拓扑 → 绘制场图。

```matlab
%% ========== 自定义时空晶体的完整分析流程 ==========
clear;  close all;
rootDir = stm_init();

%% ===== 第一步：设计参数 =====
% 空间周期中 1/3 的区域以正弦调制，调制频率 0.15
p.Lambda = 1;  p.c0 = 1;  p.g = 2*pi/p.Lambda;
p.epsBack = 2;
p.epsMod = 6;  p.modDepth = 0.5;
p.xModStart = 2*p.Lambda/3;  p.xModEnd = p.Lambda;
p.OmegaBar = 0.15;
p.Omega = 2*pi*p.c0/p.Lambda * p.OmegaBar;
p.T = 2*pi/p.Omega;

%% ===== 第二步：定义材料函数和 Fourier 系数 =====
epsilon_fun = @(x,t) custom_epsilon(x, t, p);
epsCoeff = @(m,n) custom_eps_coeff(m, n, p);

%% ===== 第三步：构建 ST-PWE 系统 =====
Nspace = 12;  Mtime = 1;
sys = stpwe_build_system(epsCoeff, [], Nspace, Mtime, p.g, p.Omega);
fprintf('System matrix size: %d x %d\n', 2*sys.S, 2*sys.S);

%% ===== 第四步：计算能带结构 =====
Nk = 301;  kBar = linspace(-0.5, 0.5, Nk);
kValues = p.g * kBar;  fMax = 0.6;
imagTol = 3e-3;  weightThresh = 0.02;

allK = [];  allF = [];  allW = [];
for ik = 1:Nk
    sol = stpwe_solve_omega(sys, kValues(ik));
    fBar = sol.omega / (p.g*p.c0);
    keep = isfinite(fBar) & real(fBar)>=0 & real(fBar)<=fMax ...
        & abs(imag(fBar))<=imagTol & sol.m0Weight>=weightThresh;
    allK = [allK; repmat(kBar(ik), sum(keep), 1)];
    allF = [allF; real(fBar(keep))];
    allW = [allW; sol.m0Weight(keep)];
end

% 静态参考能带
epsCoeff0 = @(n) epsCoeff(0, n);
fStatic = stpwe_static_bands(kValues, epsCoeff0, Nspace, p.g, p.c0, 6);

figure('Color','w', 'Position', [100 100 800 500]);  hold on;
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
colormap(flipud(gray));  colorbar;  grid on;

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
kTargetBar = 0.15;  fTargetBar = 0.20;
mode = stpwe_select_mode(sys, p.g*kTargetBar, p.g*p.c0*fTargetBar, [0 0.7*p.g*p.c0]);

xField = linspace(0, 3*p.Lambda, 301);
tField = linspace(0, 2*p.T, 201);
field = stpwe_reconstruct_field(mode, xField, tField, false);

figure('Color','w', 'Position', [200 200 600 400]);
imagesc(xField/p.Lambda, tField/p.T, field);
set(gca, 'YDir', 'normal');  caxis([-1 1]);
colormap(stm_redblue(256));
xlabel('x/\Lambda');  ylabel('t/T');
colorbar;  hold on;
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

## 11. 科研计算的收敛顺序

建议每次修改模型后依次检查：

1. `test_smoke` 是否通过。
2. 材料 Fourier 系数是否满足实函数的共轭关系 $p_{-m,-n}=p_{m,n}^*$。
3. 增大空间截断 `Nspace`，确认主要带和复带隙收敛。
4. 增大时间截断 `Mtime`。时间方波比正弦调制需要更多谐波。
5. 用 `demos/demo03` 一类的独立 TMM 结果校验时间 PWE。
6. FDTD 依次减半 $\Delta x$ 与 $\Delta t$，检查带边、增长率和散射幅度。
7. 拓扑量同时检查能隙、Wilson 链接、网格和 Fourier 截断。

只看一张"像论文"的彩图不足以证明数值正确；复谱、独立方法交叉验证和收敛测试应一起保留。`demos/demo12_ptc_convergence_audit` 是最小模板：它同时保存图片和误差数据。

---

## 12. 适用范围与尚未包含的物理

当前代码针对一维、各向同性、局域、非色散的时空介质。以下问题需要扩展：

- 二维/三维矢量 Maxwell ST-PWE；
- 各向异性、双各向异性或磁电耦合；
- Lorentz/Drude 色散和因果时变色散；
- 严格 CPML 与开放系统准正常模；
- Floquet 异常拓扑相的完整 Rudner winding；
- 简并复合带的非 Abelian 双正交 Wilson 环；
- 非线性调制与泵浦耗尽。

库中的传输线、时变超表面、合成光子网格、共形移动边界和量子紧束缚时空晶体还需要独立状态方程。它们不应通过替换一个 `epsFun` 就被称为已复现。逐篇判断见 [paper-map.md](paper-map.md)。

尤其要注意：时间 TMM 不能直接替代空间也周期变化的 ST-PWE；反之，有限 Fourier 截断的 PWE 也不应被当作突变时间层的"精确解"。

详见 [roadmap.md](roadmap.md) 中的课题扩展建议。

---

## 13. 核心参考文献

1. J. Park and B. Min, "Spatiotemporal plane wave expansion method for arbitrary space-time periodic photonic media," *Optics Letters* **46**, 484–487 (2021). DOI: [10.1364/OL.411622](https://doi.org/10.1364/OL.411622).
2. D. Ramaccia, A. Alù, A. Toscano, and F. Bilotti, "Temporal multilayer structures for designing higher-order transfer functions using time-varying metamaterials," *Applied Physics Letters* **118**, 101901 (2021). DOI: [10.1063/5.0042567](https://doi.org/10.1063/5.0042567).
3. G. R. Morgenthaler, "Velocity modulation of electromagnetic waves," *IRE Trans. Microwave Theory Tech.* **6**, 167 (1958).
4. T. Fukui, Y. Hatsugai & H. Suzuki, "Chern numbers in discretized Brillouin zones," *J. Phys. Soc. Jpn.* **74**, 1674 (2005).
5. Park & Min, "Space-time photonic crystals and topological phenomena," *Physical Review B*, 2024.
6. E. Lustig et al., "Topology of photonic time-crystals," arXiv:1803.08731.

其余 18 份库文件的覆盖状态见 [paper-map.md](paper-map.md)。

---

## 附录：函数索引

### `core/` — ST-PWE 核心计算

| 函数                                  | 输入要点                                               | 输出要点                                         | 依赖                     |
| ------------------------------------- | ------------------------------------------------------ | ------------------------------------------------ | ------------------------ |
| `stpwe_build_system`                | `(epsCoeff, muCoeff, Nspace, Mtime, g, Omega)`       | `sys` 结构体                                   | —                       |
| `stpwe_solve_omega`                 | `(sys, k)`                                           | 复$\omega$，`R`/`L` 本征向量，`m0Weight` | `stpwe_build_system`   |
| `stpwe_solve_k`                     | `(sys, omega)`                                       | 复$k$，`R`/`L` 本征向量                    | `stpwe_build_system`   |
| `stpwe_fold_frequency`              | `(omega, Omega)`                                     | 折叠到$[-\Omega/2,\Omega/2)$ 的频率            | —                       |
| `stpwe_select_mode`                 | `(sys, kTarget, omegaTarget, [omegaWindow])`         | `mode` 结构体                                  | `stpwe_solve_omega`    |
| `stpwe_reconstruct_field`           | `(mode, x, t, [includeGrowth])`                      | 归一化$E(x,t)$ 矩阵                            | —                       |
| `stpwe_static_bands`                | `(kValues, epsCoeff0, Nspace, g, c0, nBands)`        | `fBands` 矩阵                                  | —                       |
| `stpwe_sample_fourier_coefficients` | `(materialFun, mOrders, nOrders, Lambda, T, Nx, Nt)` | `table` 结构体                                 | —                       |
| `stpwe_lookup_coefficient`          | `(table, m, n)`                                      | Fourier 系数值                                   | 需要`table`            |
| `stm_preset_modulated_slab`         | —                                                     | Park-Min 模型参数结构体`p`                     | —                       |
| `stm_permittivity_modulated_slab`   | `(x, t, p)`                                          | $\varepsilon(x,t)$ 矩阵                        | —                       |
| `stm_fourier_modulated_slab`        | `(m, n, p)`                                          | $\varepsilon_{mn}$ 解析系数                    | `stm_interval_fourier` |
| `stm_interval_fourier`              | `(n, xa, xb, Lambda)`                                | 子区间 Fourier 积分                              | —                       |
| `stm_redblue`                       | `(n)`                                                | $n\times 3$ 色图                               | —                       |
| `stm_stft`                          | `(signal, dt, [window], [hop], [nFFT])`              | 频谱，频率轴，时间轴                             | —                       |

### `tmm/` — 时间转移矩阵

| 函数                                 | 输入要点                                                        | 输出要点        | 依赖                                                     |
| ------------------------------------ | --------------------------------------------------------------- | --------------- | -------------------------------------------------------- |
| `temporal_interface_matrix`        | `(ε_before, μ_before, ε_after, μ_after)`                  | `MM, τ, ρ`  | —                                                       |
| `temporal_delay_matrix`            | `(ω_init, n_init, n_slab, duration)`                         | 延迟矩阵        | —                                                       |
| `temporal_crystal_monodromy`       | `(k, eps_seq, mu_seq, durations)`                             | 单值矩阵$U$   | —                                                       |
| `temporal_crystal_bands`           | `(kValues, eps_seq, mu_seq, durations)`                       | 能带结构体      | `temporal_crystal_monodromy`                           |
| `temporal_interface_matrix_jump`   | 可选跃迁律                                                      | 界面跃迁矩阵    | —                                                       |
| `temporal_finite_crystal_response` | 有限周期参数                                                    | 复方向输出      | —                                                       |
| `temporal_domain_wall_mode`        | 两个 PTC 参数                                                   | 局域态匹配      | —                                                       |
| `temporal_multilayer_tmm`          | `(ω, ε_i, μ_i, ε_slabs, μ_slabs, durations, ε_f, μ_f)` | `TM, details` | `temporal_interface_matrix`, `temporal_delay_matrix` |
| `temporal_tmm_spectrum`            | `(ω_range, ...)`                                             | 频谱结构体      | `temporal_multilayer_tmm`                              |
| `temporal_binary_eps_coeff`        | `(m, εA, εB, dutyA)`                                        | Fourier 系数值  | —                                                       |

### `fdtd/` — 时域有限差分

| 函数          | 输入要点                 | 输出要点                          | 依赖 |
| ------------- | ------------------------ | --------------------------------- | ---- |
| `fdtd1d_db` | `cfg` 结构体（见 §9） | `out` 结构体（含 E/H/能量历史） | —   |

### `topology/` — 拓扑不变量

| 函数                       | 输入要点                                      | 输出要点             | 依赖                   |
| -------------------------- | --------------------------------------------- | -------------------- | ---------------------- |
| `stpwe_bz_sewing_matrix` | `sys`                                       | 置换矩阵             | `stpwe_build_system` |
| `stpwe_track_band`       | `(sys, kGrid, omegaSeed, options)`          | `band` 结构体      | `stpwe_solve_omega`  |
| `zak_phase_biorthogonal` | `(rightStates, leftStates, metric, sewing)` | Zak 相位，`info`   | —                     |
| `fhs_chern_number`       | `(rightStates, [leftStates], [metric])`     | Chern 数，曲率，链接 | —                     |

### 顶层脚本

| 函数                           | 功能                   | 依赖                                          |
| ------------------------------ | ---------------------- | --------------------------------------------- |
| `stm_init()`                 | 设置路径，创建输出目录 | —                                            |
| `stm_run_examples(runHeavy)` | 按顺序运行所有示例     | 全部                                          |
| `stm_selftest()`             | 自检测试               | `core/`, `tmm/`, `fdtd/`, `topology/` |

### 示例脚本 (`examples/`)

| 脚本                                          | 功能                     | 依赖                     |
| --------------------------------------------- | ------------------------ | ------------------------ |
| `example_floquet_bands_and_fields(quality)` | 重现 Park-Min 论文 Fig.2 | `core/`                |
| `example_complex_gaps()`                    | 复频率/复波数带隙        | `core/`                |
| `example_ptc_pwe_vs_tmm()`                  | PWE 与 TMM 交叉验证      | `core/`, `tmm/`      |
| `example_tmm_multilayer()`                  | 四类时间多层结构         | `tmm/`                 |
| `example_fdtd_interface()`                  | 单个时间界面的波分裂     | `fdtd/`, `tmm/`      |
| `example_fdtd_wavepacket()`                 | 时空晶体中的波包演化     | `fdtd/`, `core/`     |
| `example_zak_phase()`                       | 双正交 Zak 相位计算      | `core/`, `topology/` |
| `example_chern_pump()`                      | Rice-Mele 泵 Chern 数    | `topology/`            |

### 演示脚本 (`demos/`)

| 脚本                                           | 功能                          | 依赖                     |
| ---------------------------------------------- | ----------------------------- | ------------------------ |
| `demo01_reproduce_fig2_stpwe`                | 重构 Park-Min 论文图 2        | `core/`                |
| `demo02_complex_frequency_and_momentum_gaps` | 复频率/复波数带隙             | `core/`                |
| `demo03_ptc_bands_pwe_vs_tmm`                | PWE 与精确单周期 TMM 交叉验证 | `core/`, `tmm/`      |
| `demo04_temporal_multilayer_tmm`             | Ramaccia 四类时间多层结构     | `tmm/`                 |
| `demo05_fdtd_temporal_interface`             | 单时间界面波分裂              | `fdtd/`, `tmm/`      |
| `demo06_fdtd_spacetime_wavepacket`           | 时空晶体中的波包演化          | `core/`, `fdtd/`     |
| `demo07_zak_phase_stpwe`                     | 双正交 Zak 相位与诊断         | `core/`, `topology/` |
| `demo08_fhs_chern_thouless_pump`             | Rice-Mele 泵 FHS Chern 数     | `topology/`            |
| `demo09_finite_ptc_order_and_phase`          | AB/BA 单胞的有限周期响应      | `tmm/`                 |
| `demo10_temporal_domain_wall_mode`           | 时间畴壁增长/衰减模匹配       | `tmm/`                 |
| `demo11_coherent_time_interface`             | 双输入相干碰撞与跃迁律        | `tmm/`, `fdtd/`      |
| `demo12_ptc_convergence_audit`               | 二元时间晶体 PWE 收敛审计     | `core/`, `tmm/`      |
