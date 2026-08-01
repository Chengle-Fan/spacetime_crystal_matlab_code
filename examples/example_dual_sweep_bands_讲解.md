# `example_dual_sweep_bands.m` 完整讲解

## 一、物理背景：什么是光子时间晶体（PTC）？

### 1.1 从空间周期到时空周期

- **普通光子晶体**：介电常数只在空间上周期变化，$\varepsilon(x+\Lambda) = \varepsilon(x)$。Bloch 定理告诉我们 $k$ 是好量子数，求解 $\omega(k)$ 得到能带结构。在某些频率区间无解，形成**频率带隙（frequency bandgap）**。
- **光子时间晶体**：介电常数只在时间上周期变化，$\varepsilon(t+T) = \varepsilon(t)$。此时 $\omega$ 是好量子数，求解 $k(\omega)$。在某些动量区间无传播解，形成**动量带隙（momentum bandgap / k-gap）**。
- **时空光子晶体**（本代码研究对象）：$\varepsilon(x+\Lambda, t+T) = \varepsilon(x,t)$，同时具有空间和时间周期性。两种带隙可以共存。

### 1.2 为什么需要"双向扫描"？

正是因为时空周期介质里 $\omega$ 和 $k$ 不再有一一对应关系：

- **固定 k 求 ω**（常规思路）：在频率带隙处出现复数 ω，虚部代表衰减/放大。
- **固定 ω 求 k**（时间晶体思路）：在动量带隙处出现复数 k，虚部代表空间衰逝。

两种扫描互为补充。把它们画在同一张 $(k, \omega)$ 平面上，重合的地方是传播模，分歧的地方就是各种带隙。这就是本文代码的核心目的。

---

## 二、数学模型：ST-PWE（时空平面波展开）

### 2.1 一维 Maxwell 方程组

在 1D 中，设 $\mu=1$：

$$
\partial_x H = -\partial_t(\varepsilon E), \quad \partial_x E = -\partial_t H
$$

### 2.2 双 Fourier 展开

由于介质是时空周期的，将场和介电常数都做双 Fourier 展开：

$$
\varepsilon(x,t) = \sum_{m,n} \varepsilon_{m,n} \, e^{i n g x - i m \Omega t}
$$

$$
E(x,t) = e^{i(kx - \omega t)} \sum_{m,n} E_{m,n} \, e^{i n g x - i m \Omega t}
$$

$$
H(x,t) = e^{i(kx - \omega t)} \sum_{m,n} H_{m,n} \, e^{i n g x - i m \Omega t}
$$

其中 $g = 2\pi/\Lambda$ 是空间倒格矢，$\Omega = 2\pi/T$ 是调制频率。

### 2.3 两个本征值问题

代入 Maxwell 方程后得到：

**固定 k → 求 ω（Park-Min Eq.6，广义本征值问题）**：

$$
A(k) \Phi = \omega B \Phi
$$

$$
A = \begin{bmatrix} K & -W C_\mu \\ -W C_\varepsilon & K \end{bmatrix}, \quad B = \begin{bmatrix} 0 & C_\mu \\ C_\varepsilon & 0 \end{bmatrix}
$$

其中 $K = kI + G$（$G = \text{diag}(n g)$），$W = \text{diag}(m\Omega)$，$C_\varepsilon$ 和 $C_\mu$ 是 Fourier 系数的卷积矩阵。

**固定 ω → 求 k（Park-Min Eq.7，标准本征值问题）**：

$$
A'(\omega) \Phi = k \Phi
$$

$$
A' = \begin{bmatrix} -G & (\omega I + W) C_\mu \\ (\omega I + W) C_\varepsilon & -G \end{bmatrix}
$$

---

## 三、代码逐段详解

### 3.1 第 1-13 行：函数签名与帮助文档

```matlab
function example_dual_sweep_bands(quality)
%EXAMPLE_DUAL_SWEEP_BANDS 分别用"扫 k 求 ω"和"扫 ω 求 k"计算能带，并合并对比。
```

- 唯一输入参数 `quality`：`'quick'`（快速测试，默认）或 `'paper'`（论文级精度）。
- 产生三张图，分别对应固定 k 扫描、固定 ω 扫描、叠加对比。
- 图片自动保存到 `output/` 目录。

### 3.2 第 17-37 行：参数设置

```matlab
if nargin < 1
    quality = 'quick';
end
rootDir = stm_init();
```

`stm_init()`（[stm_init.m](../stm_init.m)）将所有子目录（`core/`、`tmm/`、`fdtd/`、`topology/`、`examples/`、`tests/`）加入 MATLAB 搜索路径，并创建 `output/` 目录。

```matlab
p = stm_preset_modulated_slab();
```

**这是最关键的结构体**，定义了一个具体的时空介质（[stm_preset_modulated_slab.m](../core/stm_preset_modulated_slab.m)）：

| 参数            | 值                                        | 含义                                                      |
| --------------- | ----------------------------------------- | --------------------------------------------------------- |
| `p.Lambda`    | 1                                         | 归一化空间周期                                            |
| `p.c0`        | 1                                         | 归一化光速                                                |
| `p.g`         | $2\pi/\Lambda = 2\pi$                   | 空间倒格矢                                                |
| `p.eps1`      | 2                                         | 静态区域的介电常数                                        |
| `p.epsc`      | 6                                         | 调制区域的背景介电常数                                    |
| `p.modDepth`  | 0.6                                       | 调制深度                                                  |
| `p.xModStart` | $3\Lambda/4 = 0.75$                     | 调制区起始位置                                            |
| `p.xModEnd`   | $\Lambda = 1$                           | 调制区结束位置                                            |
| `p.OmegaBar`  | 0.20                                      | 归一化调制频率$\bar{\Omega} = \Omega\Lambda/(2\pi c_0)$ |
| `p.Omega`     | $2\pi c_0/\Lambda \times 0.20 = 0.4\pi$ | 物理调制角频率                                            |
| `p.T`         | $2\pi/\Omega$                           | 调制周期                                                  |

**物理图像**：一个空间周期 $\Lambda=1$ 内：

- $[0,\; 0.75)$ — 静态介质，$\varepsilon = 2$
- $[0.75,\; 1)$ — 时间调制介质，$\varepsilon = 6 \times \bigl(1 + 0.6\sin(\Omega t)\bigr)$

```matlab
switch lower(quality)
    case 'paper'
        Nspace = 20;    % 空间 Fourier 截断：保留 n = -20..20（共 41 个谐波）
        Nk = 181;       % k 扫描点数
        Nw = 161;       % ω 扫描点数
    case 'quick'
        Nspace = 10;    % 空间 Fourier 截断：保留 n = -10..10（共 21 个谐波）
        Nk = 101;       % k 扫描点数
        Nw = 81;        % ω 扫描点数
end
Mtime = 1;              % 时间 Fourier 截断：保留 m = -1..1（共 3 个谐波）
```

`Mtime=1` 只需要 3 个时间谐波（$m=-1, 0, +1$），因为调制是纯正弦的 $\sin(\Omega t)$，介电常数在时间方向只有 DC 分量 + 基频分量（$m=0, \pm 1$）。

### 3.3 第 40-41 行：构造系统矩阵

```matlab
epsCoeff = @(m,n) stm_fourier_modulated_slab(m,n,p);
sys = stpwe_build_system(epsCoeff, [], Nspace, Mtime, p.g, p.Omega);
```

#### 3.3.1 `stm_fourier_modulated_slab(m,n,p)` — 解析计算 Fourier 系数

该函数将空间积分拆为静态区和调制区（[stm_fourier_modulated_slab.m](../core/stm_fourier_modulated_slab.m)）：

```matlab
staticPart = stm_interval_fourier(n, 0, p.xModStart, p.Lambda);          % [0, 0.75) 的 Fourier 积分
modPart    = stm_interval_fourier(n, p.xModStart, p.xModEnd, p.Lambda);  % [0.75, 1) 的积分
```

`stm_interval_fourier`（[stm_interval_fourier.m](../core/stm_interval_fourier.m)）计算：

$$
\frac{1}{\Lambda}\int_{x_a}^{x_b} e^{-i n g x} \, dx
$$

- 当 $n=0$：直接返回 $(x_b - x_a)/\Lambda$（即区间长度占比）
- 当 $n \neq 0$：返回封闭形式的解析积分值

$$
\text{value} = \frac{e^{-i n g x_b} - e^{-i n g x_a}}{-i n g \Lambda}
$$

然后按时间谐波 $m$ 分配系数：

| $m$ | 系数$\varepsilon_{m,n}$                                                        | 物理含义                |
| ----- | -------------------------------------------------------------------------------- | ----------------------- |
| 0     | $\varepsilon_1 \times \text{staticPart} + \varepsilon_c \times \text{modPart}$ | 直流（时间平均）分量    |
| +1    | $+i \cdot (\varepsilon_c \cdot \text{modDepth}/2) \times \text{modPart}$       | $e^{-i\Omega t}$ 分量 |
| -1    | $-i \cdot (\varepsilon_c \cdot \text{modDepth}/2) \times \text{modPart}$       | $e^{+i\Omega t}$ 分量 |
| 其他  | 0                                                                                | 纯正弦调制只有基频      |

> **注意** $m=\pm 1$ 系数的 $\pm i$ 相位差——这来自 $\sin(\Omega t) = (e^{i\Omega t} - e^{-i\Omega t})/(2i)$ 的 Fourier 展开。

#### 3.3.2 `stpwe_build_system` — 构造卷积矩阵

该函数做两件事（[stpwe_build_system.m](../core/stpwe_build_system.m)）：

**(1) 建立谐波索引网格**

```matlab
[NN, MM] = ndgrid(-Nspace:Nspace, -Mtime:Mtime);
nList = NN(:);   % 所有空间谐波索引
mList = MM(:);   % 所有时间谐波索引
S = numel(nList); % 总谐波数
```

以 `'quick'` 模式为例：$(2 \times 10 + 1) \times (2 \times 1 + 1) = 21 \times 3 = 63$ 个谐波。

**(2) 填充卷积矩阵 Ceps 和 Cmu**

```matlab
for row = 1:S
    for col = 1:S
        dm = mList(row) - mList(col);
        dn = nList(row) - nList(col);
        Ceps(row, col) = epsCoeff(dm, dn);
    end
end
```

这是一个**双重 Toeplitz 结构**：矩阵元只依赖于行列索引的**差值** $(m_{\text{row}}-m_{\text{col}},\; n_{\text{row}}-n_{\text{col}})$。这对应 Fourier 空间的卷积——在实空间中 $\varepsilon(x,t)E(x,t)$ 的逐点乘积，在 Fourier 空间中变成了 $\varepsilon_{m,n}$ 与 $E_{m,n}$ 之间的离散卷积。

**(3) 存储导数矩阵和 B 矩阵**

```matlab
sys.G = diag(nList * g);     % 空间导数对角矩阵
sys.W = diag(mList * Omega); % 时间导数对角矩阵
sys.Bomega = [Z, Cmu; Ceps, Z]; % 用于 ω 本征问题的 B 矩阵
```

### 3.4 第 44-51 行：公共参数

```matlab
fMax = 0.82;              % 频率显示上限（归一化 fBar）
kHalfRange = 0.5;         % kBar 显示半宽
weightThreshold = 0.035;  % m=0 参与度阈值（仅 ω 扫描使用）
imagTolKscan = 5e-2;      % 扫 k 时容许的 |Im(ω)| 上限
imagTolWscan = 0.12;      % 扫 ω 时容许的 |Im(k)| 上限
```

关键参数说明：

- **`fMax = 0.82`**：只显示归一化频率小于 0.82 的能带（前几个能带），高于此频率的模不显示。
- **`kHalfRange = 0.5`**：k 的范围 $[-0.5, 0.5]$（以 $g=2\pi/\Lambda$ 为单位），恰好覆盖整个第一 Brillouin 区。
- **`weightThreshold = 0.035`**：`m0Weight` 衡量 $m=0$（时间直流分量）的场能量占总能量的比例。低于此阈值的模被滤除，因为这些主要是高频 Floquet 边带而非我们关心的基频传播模。
- **`imagTolKscan = 0.05`**：在固定 k 扫描中，复数 ω（带隙内的衰逝模）的虚部上限。设 0.05 是为了保留近带隙的衰逝模用于着色，同时过滤纯数值伪解。
- **`imagTolWscan = 0.12`**：ω 扫描时 |Im(k)| 的阈值设得更大（0.12 vs 0.05），因为动量带隙中的 Im(k) 天然就比频率带隙中的 Im(ω) 更大。

### 3.5 第 54-77 行：扫 k 求 ω（Sweep 1 — 常规能带）

```matlab
kBarSweep = linspace(-kHalfRange, kHalfRange, Nk);
```

归一化波数 kBar $= k\Lambda/(2\pi)$ 的范围是 $[-0.5, 0.5]$，覆盖整个第一 Brillouin 区。

```matlab
for ik = 1:Nk
    sol = stpwe_solve_omega(sys, p.g * kBarSweep(ik));
```

每次调用 `stpwe_solve_omega` 执行以下步骤（[stpwe_solve_omega.m](../core/stpwe_solve_omega.m)）：

**1. 构造矩阵并求解广义本征值问题**

```matlab
K = k * sys.I + sys.G;   % 对角矩阵 diag(k + n*g)
A = [K,       -sys.W*sys.Cmu; ...
    -sys.W*sys.Ceps,  K         ];
B = sys.Bomega;           % [Z, Cmu; Ceps, Z]
[R, D, L] = eig(A, B);
omega = diag(D);
```

矩阵规模为 $2S \times 2S$。前 $S$ 行对应电场分量，后 $S$ 行对应磁场分量。

**2. 双正交归一化**

```matlab
for q = 1:numel(omega)
    overlap = L(:,q)' * B * R(:,q);
    if isfinite(overlap) && abs(overlap) > 1e-12
        R(:,q) = R(:,q) / overlap;   % 归一化使 L'*B*R = 1
    else
        % 自正交或退化情况：分别归一化
        R(:,q) = R(:,q) / norm(R(:,q));
        L(:,q) = L(:,q) / norm(L(:,q));
    end
end
```

对于非 Hermitian 系统，左右本征矢满足双正交关系。这个归一化确保后续的 m0Weight 计算有意义。

**3. 计算 m0Weight（时间直流分量占比）**

```matlab
mask0 = (sys.mList == 0);   % 选出 m=0 的时间谐波分量
E = R(1:S, :);              % 电场分量（前 S 行）
H = R(S+1:end, :);          % 磁场分量（后 S 行）

numerator   = sum(|E(m=0)|² + |H(m=0)|²);   % m=0 分量的场能量
denominator = sum(|E|² + |H|²);              % 总场能量
m0Weight    = numerator / denominator;
```

$$
\text{m0Weight} = \frac{\sum_{n} \bigl(|E_{n,m=0}|^2 + |H_{n,m=0}|^2\bigr)}{\sum_{n,m} \bigl(|E_{n,m}|^2 + |H_{n,m}|^2\bigr)}
$$

**物理含义**：m0Weight 越大，说明该模式越接近静态介质的常规 Bloch 模；m0Weight 很小意味着该模主要是高频 Floquet 边带。

**4. 筛选解**

```matlab
fBar = sol.omega / (p.g * p.c0);   % 转为归一化频率

keep = isfinite(fBar) ...
    & real(fBar) >= 0 & real(fBar) <= fMax ...     % 只取正频率，且在显示范围内
    & abs(imag(fBar)) <= imagTolKscan ...            % 虚部不太大（非严重衰减/放大模）
    & sol.m0Weight >= weightThreshold;               % 有足够的 DC 时间分量
```

- `real(fBar) >= 0`：只保留正频率解（时间反演对称性保证负频率是冗余的）
- `imagTolKscan` 滤除虚部过大的数值噪声，同时保留近带隙的微弱 Im(ω) 用于着色
- `weightThreshold` 滤除高次 Floquet 边带的伪解

**5. 汇总结果**

```matlab
kscan_k  = vertcat(solKscan_k{:});    % 所有 k 值
kscan_fr = vertcat(solKscan_fr{:});   % 所有 Re(ω) 值
kscan_fi = vertcat(solKscan_fi{:});   % 所有 Im(ω) 值
```

### 3.6 第 80-102 行：扫 ω 求 k（Sweep 2 — 复动量能带）

```matlab
fBarSweep = linspace(0, fMax, Nw);
```

归一化频率从 0 到 0.82。

```matlab
for iw = 1:Nw
    sol = stpwe_solve_k(sys, p.g * p.c0 * fBarSweep(iw));
```

`stpwe_solve_k`（[stpwe_solve_k.m](../core/stpwe_solve_k.m)）与 `stpwe_solve_omega` 的关键区别：

| 特性       | `stpwe_solve_omega`     | `stpwe_solve_k`     |
| ---------- | ------------------------- | --------------------- |
| 本征值类型 | 广义本征值`eig(A, B)`   | 标准本征值`eig(A)`  |
| 本征值     | $\omega$                | $k$                 |
| 归一化     | $L^H B R = 1$（双正交） | $L^H R = 1$（标准） |
| m0Weight   | 计算                      | 不计算                |
| 矩阵构造   | 对角块含$k$             | 对角块含$\omega$    |

```matlab
OW = omega * sys.I + sys.W;   % diag(omega + m*Omega)
A = [-sys.G,       OW*sys.Cmu; ...
      OW*sys.Ceps, -sys.G      ];
[R, D] = eig(A);
k = diag(D);
```

**筛选条件**：

```matlab
kNorm = sol.k / p.g;   % kBar = k/g

keep = isfinite(kNorm) ...
    & real(kNorm) >= -kHalfRange & real(kNorm) <= kHalfRange ...
    & abs(imag(kNorm)) <= imagTolWscan;
```

- 这里不筛选 `m0Weight`——因为固定 ω 时，时间频率已给定，"DC 分量"概念不直观
- Re(k) 取正也取负，覆盖 $[-0.5, 0.5]$ 整个 Brillouin 区
- `imagTolWscan = 0.12` 比 k 扫描的虚部阈值（0.05）更宽松

---

## 四、三张图详解

### 图 1（第 105-127 行）：固定 k 能带图（2D）

```
scatter(k, Re(ω), 颜色 = |Im(ω)|)
```

- **横轴**：$\bar{k} = k\Lambda/(2\pi)$（归一化 Bloch 波数）
- **纵轴**：$\bar{f} = \omega\Lambda/(2\pi c)$（归一化频率）
- **颜色**：
  - 深蓝：$|\text{Im}(\omega)| \approx 0$ → 纯传播模
  - 青 → 黄 → 红：$|\text{Im}(\omega)|$ 逐渐增大 → 近带隙区域，模开始衰减

颜色映射是自定义的双色渐变：

```matlab
colormap([linspace(0,1,256)', linspace(0.447,0,256)', linspace(0.741,0,256)'])
%        R: 0→1           G: 0.447→0              B: 0.741→0
%        深蓝 → 青 → 红
```

这是固体物理中最常见的"能带图"——纵轴频率、横轴波数。颜色编码的 |Im(ω)| 额外标示了频率带隙的位置和强弱。

### 图 2（第 130-183 行）：固定 ω 复动量能带图（3D）

```
scatter3(Re(k), Im(k), ω)
```

**三个轴**：

- **X 轴**：$\text{Re}(k)\Lambda/(2\pi)$
- **Y 轴**：$\text{Im}(k)\Lambda/(2\pi)$
- **Z 轴**：$\omega\Lambda/(2\pi c)$

**四种视觉元素**：

| 元素         | 外观                   | 含义                                   |
| ------------ | ---------------------- | -------------------------------------- |
| 绿色散点     | Im(k)=0 平面内的绿色点 | 传播模（纯实数 k）                     |
| 橙色散点     | Im(k)≠0 处的橙色点    | 复动量模（k 有非零虚部），指示动量带隙 |
| 灰色半透明面 | Im(k)=0 参考面         | 传播面基底，复模偏离此面               |
| 灰色投影点   | 所有点投影到 Im(k)=0   | 辅助确认能带在传播面的位置             |
| 黑色竖线     | k=0 处的参考线         | Γ 点标记                              |

**物理含义**：这是光子时间晶体最标志性的图像。在动量带隙（k-gap）处，给定 $\omega$ 找不到纯实数的 k 解，而是出现复数 k。复数 k 的虚部 $\text{Im}(k)$ 代表空间上的指数衰减/增长——这就是"动量禁带"的直观表现。能带在 $\text{Im}(k) \neq 0$ 处形成偏离传播面的环路结构，是 PTC 区别于普通光子晶体的核心特征。

**视角**：`view(-38, 24)` 从左前下方观察，可以同时看到传播面上的能带结构和 Im(k) 方向的偏离。

**动态 Y 轴范围**：

```matlab
imSpanK = max([abs(min(wscan_ki)), abs(max(wscan_ki)), 0.01]) * 1.3;
ylim(ax2, [-imSpanK, imSpanK]);
```

自适应调整 Im(k) 轴的显示范围，确保复模的偏离清晰可见。

### 图 3（第 188-210 行）：叠加对比图（2D）

```
蓝色 scatter: k 扫描的 (k, Re(ω))
绿色 scatter: ω 扫描的 (Re(k), ω)
```

两种互补方法的结果画在同一张 $(k,\omega)$ 平面上：

- **重合区域**（蓝绿交叠）：两种方法都找到了实数解 → **传播模**，结果互相印证
- **蓝色独有区域**（有 k 但该频率附近 ω 扫描找不到传播的 k）：可能对应**频率带隙**的边缘——k 扫描找到的解 Im(ω) 较小仍被保留，但 ω 扫描在该频率附近找不到 Im(k)≈0 的解
- **绿色独有区域**（有 ω 但该 k 附近 k 扫描找不到传播的 ω）：可能对应**动量带隙**的边缘
- **两者都不在的区域**：完全禁带（full gap），两种方法都找不到传播解

> **核心思想**：由于滤波阈值（`imagTolKscan`、`imagTolWscan`、`weightThreshold`）的存在，两种方法得到的"传播解"集合不会完全相同。重合程度是两种数值方法一致性的直观验证。

---

## 五、输出与保存（第 213-244 行）

### 5.1 控制台统计

```matlab
fprintf('\n=== Dual-sweep summary ===\n');
fprintf('Figure 1: Fixed-k scan — %d k-points, %d modes  [2D]\n', Nk, numel(kscan_fr));
fprintf('Figure 2: Fixed-w scan — %d w-points, %d modes (%d prop + %d complex)  [3D]\n', ...
    Nw, numel(wscan_f), sum(isPropK), sum(isGapK));
fprintf('Figure 3: Overlay  [2D, (k,w) plane]\n');
```

打印每个扫描找到的模数量、传播模 vs 复模的比例。

### 5.2 图片保存

```matlab
tags = {'fix_k_scan', 'fix_w_scan_3d', 'overlay'};
for i = 1:3
    outputFile = fullfile(rootDir, 'output', ...
        ['example_dual_sweep_bands_' tags{i} '_' lower(quality) '.png']);
    save_example_figure(figs(i), outputFile);
end
```

生成三张 PNG 文件（以 `'quick'` 为例）：

- `output/example_dual_sweep_bands_fix_k_scan_quick.png`
- `output/example_dual_sweep_bands_fix_w_scan_3d_quick.png`
- `output/example_dual_sweep_bands_overlay_quick.png`

### 5.3 辅助函数 `save_example_figure`

```matlab
function save_example_figure(fig, outputFile)
    try
        exportgraphics(fig, outputFile, 'Resolution', 220);
    catch
        print(fig, outputFile, '-dpng', '-r220');
    end
end
```

先尝试 `exportgraphics`（MATLAB R2020a+ 的现代导出函数，支持更好的抗锯齿），失败则回退到传统的 `print` 函数。分辨率为 220 DPI。

---

## 六、关键概念总结

| 概念              | 符号/公式                                                                                                         | 物理含义                              |
| ----------------- | ----------------------------------------------------------------------------------------------------------------- | ------------------------------------- |
| 空间倒格矢        | $g = 2\pi/\Lambda$                                                                                              | 空间周期性的 Fourier 基本频率         |
| 调制频率          | $\Omega$                                                                                                        | 时间周期性的角频率                    |
| 归一化波数        | $\bar{k} = k\Lambda/(2\pi)$                                                                                     | 以 Brillouin 区宽度为单位的无量纲波数 |
| 归一化频率        | $\bar{f} = \omega\Lambda/(2\pi c)$                                                                              | 无量纲频率                            |
| 第一 Brillouin 区 | $\bar{k} \in [-0.5, 0.5]$                                                                                       | Bloch 波数的不等价范围                |
| m0Weight          | $\frac{\sum_n \vert E_{n,0}\vert^2+\vert H_{n,0}\vert^2}{\sum_{n,m} \vert E_{n,m}\vert^2+\vert H_{n,m}\vert^2}$ | $m=0$ 时间谐波的能量占比            |
| 频率带隙          | $\text{Im}(\omega) \neq 0$ at fixed $k$                                                                       | 空间周期性的传统带隙                  |
| 动量带隙          | $\text{Im}(k) \neq 0$ at fixed $\omega$                                                                       | 时间周期性的特有带隙                  |

---

## 七、数据流图

```
stm_preset_modulated_slab()
        │
        ▼
  p (参数结构体: ε分布, Λ, Ω, …)
        │
        ▼
  epsCoeff = @(m,n) stm_fourier_modulated_slab(m,n,p)
        │
        ▼
  stpwe_build_system(epsCoeff, [], Nspace, Mtime, g, Omega)
        │
        ▼
  sys (系统矩阵: Ceps, Cmu, G, W, Bomega, …)
        │
        ├──────────────────────────────────┐
        ▼                                  ▼
  Sweep 1: 扫 k 求 ω                    Sweep 2: 扫 ω 求 k
  stpwe_solve_omega(sys, k_i)           stpwe_solve_k(sys, ω_j)
        │                                  │
        ▼                                  ▼
  (k, Re(ω), Im(ω), m0Weight)         (Re(k), Im(k), ω)
        │                                  │
        ├──────┬──────┐                    ├──────┬──────┐
        ▼      ▼      ▼                    ▼      ▼      ▼
      图1    图3     --                  图2    图3     --
   (k-Re(ω)  叠加图                    (3D复k)  叠加图
    着色Im)
```

---

## 八、调用方式

```matlab
% 快速测试（默认）
example_dual_sweep_bands('quick')

% 论文级精度
example_dual_sweep_bands('paper')

% 无参数等同于 'quick'
example_dual_sweep_bands
```

---

## 九、依赖文件一览

| 文件                                  | 作用                                                     |
| ------------------------------------- | -------------------------------------------------------- |
| `stm_init.m`                        | 工具箱初始化，添加路径                                   |
| `core/stm_preset_modulated_slab.m`  | 定义时空介质的物理参数                                   |
| `core/stm_fourier_modulated_slab.m` | 解析计算介电常数的时空 Fourier 系数$\varepsilon_{m,n}$ |
| `core/stm_interval_fourier.m`       | 计算空间子区间上的 Fourier 积分                          |
| `core/stpwe_build_system.m`         | 构造 ST-PWE 的系统矩阵（卷积矩阵 + 导数矩阵）            |
| `core/stpwe_solve_omega.m`          | 固定 k 求解广义本征值问题$A\Phi = \omega B\Phi$        |
| `core/stpwe_solve_k.m`              | 固定 ω 求解标准本征值问题$A'\Phi = k\Phi$             |
