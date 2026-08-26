# 05 · 时间传输矩阵（temporal TMM）精讲

> 目标读者：**会麦克斯韦方程、能带、Fourier 级数，但还没写过几行 MATLAB** 的物理学生。
> 这一章我们做一件"把空间换成时间"的事：空间传输矩阵把一层层**空间**串起来，时间传输矩阵把一段段**时间**串起来。物理直觉、2×2 矩阵、能带判据，全都平移过去。
>
> 读完这一章，你应该能做到：
> 1. 说出"空间 TMM 串层"与"时间 TMM 串时间段"的一一对应关系；
> 2. 读懂时间界面匹配矩阵 `temporal_interface_matrix` 的 τ、ρ 从哪来，为什么矩阵是对称的；
> 3. 理解单值矩阵 `temporal_crystal_monodromy` 为什么用 `expm`、为什么 `det(U)=1`；
> 4. 自己跑一个二元时间晶体的能带 + `halfTrace` 判据，并指出动量带隙在哪、多宽。

所有代码都在 `tmm/` 目录，运行前先 `startup_stm`（若忘了，翻 [01-matlab-primer.md](01-matlab-primer.md) 复习）。本章对应的 demo 是 [demo03](../demos/demo03_ptc_bands_pwe_vs_tmm.m)、[demo04](../demos/demo04_temporal_multilayer_tmm.m)、[demo09](../demos/demo09_finite_ptc_order_and_phase.m)、[demo10](../demos/demo10_temporal_domain_wall_mode.m)，建议跑完一章读一个。

---

## 一、从空间传输矩阵到时间传输矩阵

### 1.1 复习：空间 TMM 在做什么

空间多层膜（多层介质堆叠）的传输矩阵法你一定见过：每一层里 $\varepsilon$ 常数，波在里面自由传播，攒相位 $e^{\pm i k_j d_j}$；每两个层之间的界面上用 Fresnel 公式做匹配。整个堆叠就是一堆 2×2 矩阵的乘积：

```
入射 ──► [界面1][层1][界面2][层2]...[层N][界面N+1] ──► 出射
```

- 在**每一层内**：波矢 $k_j$ 变（由层内折射率决定），频率 $\omega$ 守恒；
- 在**每一个界面**：$E$、$H$ 连续（切向分量连续）。

关键点：**空间界面上匹配、层内传播**，这两种"动作"交替出现，各自用一个 2×2 矩阵表示，然后串起来。

### 1.2 时间 TMM：把"分层"换成"分时间段"

时间 TMM 处理的介质是：**空间均匀，但介电常数随时间台阶式跳变**。比如二元时间晶体（PTC）：

$$
\varepsilon(t)=
\begin{cases}
\varepsilon_A, & 0 \le t < d_A T,\\[2pt]
\varepsilon_B, & d_A T \le t < T,
\end{cases}
\qquad \varepsilon(t+T)=\varepsilon(t).
$$

把时间轴分成一段段，每一段里 $\varepsilon(t)$ 是常数。这和空间多层膜的分层完全对偶：

| 空间多层膜 | 时间多层结构 |
|---|---|
| 每一"层"在空间上常数 $\varepsilon_j$ | 每一"时间段"在时间上常数 $\varepsilon_j$ |
| 层内波矢 $k_j$ 变，$\omega$ 守恒 | 段内频率 $\omega_j$ 变，$k$ 守恒 |
| 空间界面：$E,H$ 连续 | 时间界面：$D,B$ 连续（默认模型） |
| 界面匹配：Fresnel 公式 | 界面匹配：Morgenthaler 矩阵 `temporal_interface_matrix` |
| 传输矩阵串起来：`TM = ... MM_j DM_j ...` | `temporal_multilayer_tmm` 按时间顺序串起来 |

这张对照表也写在 [本文章节 2.3](#23-通用版本-temporal_interface_matrix_jump跳变律是第-5-个参数) 里，值得背下来。**关键的对偶**是：

- 空间上介质均匀 ⇒ **$k$ 是好量子数、跨时间界面守恒**；
- 时间上介质周期 ⇒ 波在该段的频率 $\omega_j$ 由色散决定：$\omega_j = k / n_j$，其中 $n_j=\sqrt{\varepsilon_j \mu_j}$（归一化单位下 $c_0=1$）。

所以在一个时间分段常数介质里，每一段内的解就是**平面波**：

$$
E(x,t)=\bigl(E^{+}\,e^{i k x}+E^{-}\,e^{-i k x}\bigr)\,e^{-i\omega_j t},
\qquad \omega_j=\frac{k}{n_j}.
$$

两个传播方向（$+x$ 与 $-x$）都要保留——因为时间界面会把波"折"出反向分量（下一节细说）。这就是整个 2×2 方法的基本状态：**一个 2 分量列向量，装下 $+x$ 和 $-x$ 两个分支**。

### 1.3 状态怎么装：方向态 vs D/B 态

工具箱里同一套物理有两个等价的 2 分量表示，看清它们是读后续所有函数的关键：

| 名称 | 列向量 | 物理含义 | 在哪个链条里用 |
|---|---|---|---|
| **方向态**（directional） | `[E+; E-]` | $+x$ 前向振幅、$-x$ 后向振幅 | 时间界面、时间多层（`temporal_multilayer_tmm` 链） |
| **D/B 态** | `[D; B]` | $D=\varepsilon E$、$B=\mu H$ 两个场量 | 单值矩阵（`temporal_crystal_monodromy` 链） |

两个表示在任意介质 $(\varepsilon_r,\mu_r)$ 下可以互转，代码就是这两个函数（[`temporal_directional_to_db.m`](../tmm/temporal_directional_to_db.m) 与 [`temporal_db_to_directional.m`](../tmm/temporal_db_to_directional.m)，互为逆）：

```matlab
n = sqrt(epsr*mur);
% [E+;E-] -> [D;B]      D = epsr*(E+ + E-),   B = n*(E+ - E-)
state      = temporal_directional_to_db(amplitudes, epsr, mur);
% [D;B] -> [E+;E-]      E+ = 0.5*(D/epsr + B/n),  E- = 0.5*(D/epsr - B/n)
amplitudes = temporal_db_to_directional(state,      epsr, mur);
```

为什么方向态里要用 $B=n(E^{+}-E^{-})$？因为 $B=\mu H$，而 $+x$ 波 $H=E/\eta$、$-x$ 波 $H=-E/\eta$（$\eta=\sqrt{\mu/\varepsilon}$ 是波阻抗），于是 $H=(E^{+}-E^{-})/\eta$，$B=\mu H=n(E^{+}-E^{-})$。这正是"前向、后向两个分支"的物理编码。

> 💡 **两条独立链条**（tmm 目录的整体骨架，本系列 05 章的阅读主线）：
> - **D/B 链（精确、无截断）**：`temporal_crystal_monodromy` → `temporal_crystal_bands` → `temporal_finite_crystal_response` / `temporal_domain_wall_mode`。因为 $D,B$ 在时间界面上自动连续，这条链**不需要任何界面矩阵**。
> - **方向态链（频率域多层）**：`temporal_interface_matrix` → `temporal_delay_matrix` → `temporal_multilayer_tmm` → `temporal_tmm_spectrum`。这条链模拟"单输入频率下，频率在每个时间段内被压缩/拉伸"的 Ramaccia 式传输矩阵。

---

## 二、时间界面匹配矩阵

### 2.1 时间界面：$\varepsilon$ 在某一时刻突变

时间界面（temporal interface）指：**空间均匀的介质，在某一时刻 $t_0$ 整个空间的 $\varepsilon$ 从 $\varepsilon_b$ 突跳到 $\varepsilon_a$**。这是所有 PTC 物理的最小单元（物理背景见 [02-physics-background.md](02-physics-background.md)）。

波跨过这个界面的行为：

- **$k$ 守恒**：空间平移对称性在调制前后都成立，所以 $k$ 不变；
- **$\omega$ 改变**：$\omega_b=k/n_b$，$\omega_a=k/n_a$，于是 $\omega_a=(n_b/n_a)\,\omega_b$；
- **$D,B$ 连续**：Maxwell 方程里 $\partial_t D$、$\partial_t B$ 有限，跳变瞬间 $D$、$B$ 不能跳。而 $E=D/\varepsilon$、$H=B/\mu$ 会跳。

> ⚠️ 注意：$D,B$ 连续是**本构切换模型**（constitutive-switching model），不是放之四海皆准。真实开关电路可能保 $E$、可能注入电荷，跳变律得重新指定——这就是 `temporal_interface_matrix_jump` 存在的意义（见 2.3）。

### 2.2 `temporal_interface_matrix`：标准 Morgenthaler 矩阵

最常用的入口是 [`temporal_interface_matrix.m`](../tmm/temporal_interface_matrix.m)：

```matlab
[MM, tau, rho] = temporal_interface_matrix(epsBefore, muBefore, epsAfter, muAfter)
% 调用约定：after = MM * before，列向量是 [E+; E-]
```

它做的事情只有一行：**转发到 `temporal_interface_matrix_jump(..., 'DB')`**，即"$D,B$ 连续"模型。`'DB'` 之外还有别的跳变律，都收敛在 [`temporal_interface_matrix_jump.m`](../tmm/temporal_interface_matrix_jump.m) 里。

### 2.3 通用版本 `temporal_interface_matrix_jump`：跳变律是第 5 个参数

```matlab
[MM, tau, rho, details] = temporal_interface_matrix_jump( ...
    epsBefore, muBefore, epsAfter, muAfter, model)
```

前四个参数是界面两侧的 $\varepsilon,\mu$，**第 5 个参数 `model` 指定跳变律**。它定义两个跳变因子：

$$
D_{\text{after}} = \texttt{jumpD}\,D_{\text{before}},
\qquad
B_{\text{after}} = \texttt{jumpB}\,B_{\text{before}}.
$$

| `model` 取值 | 连续量 | `jumpD` | `jumpB` |
|---|---|---|---|
| `'DB'`（默认） | $D,B$ 连续 | 1 | 1 |
| `'EB'` | $E,B$ 连续 | $\varepsilon_a/\varepsilon_b$ | 1 |
| `'DH'` | $D,H$ 连续 | 1 | $\mu_a/\mu_b$ |
| `'EH'` | $E,H$ 连续 | $\varepsilon_a/\varepsilon_b$ | $\mu_a/\mu_b$ |
| `struct('jumpD',qD,'jumpB',qB)` | 自定义微观模型 | `qD` | `qB` |

代码核心就三行（抄自源码，`ratioElectric`/`ratioMagnetic` 是代码里的变量名）：

```matlab
ratioElectric  = jumpD*epsBefore/epsAfter;
ratioMagnetic  = jumpB*sqrt((epsBefore*muBefore)/(epsAfter*muAfter));
tau = 0.5*(ratioElectric + ratioMagnetic);
rho = 0.5*(ratioElectric - ratioMagnetic);
MM  = [tau, rho; rho, tau];   % 作用于 [E+; E-]
```

对默认 `'DB'` 模型（`jumpD=jumpB=1`），两个比值退化为

$$
r_E=\frac{\varepsilon_b}{\varepsilon_a},\qquad
r_M=\sqrt{\frac{\varepsilon_b\mu_b}{\varepsilon_a\mu_a}}=\frac{n_b}{n_a},
\qquad
\tau=\frac{1}{2}\bigl(r_E+r_M\bigr),\quad
\rho=\frac{1}{2}\bigl(r_E-r_M\bigr).
$$

这就是 Morgenthaler 1958 年的结果（对应 Ramaccia 等 APL 118, 101901 (2021) 的式 (5)，见 [`temporal_interface_matrix.m`](../tmm/temporal_interface_matrix.m) 注释）。$\tau$ 叫**时间透射**（同向分支、频率改变），$\rho$ 叫**时间反射**（反向分支、"逆向频率"）。

**动手验证一下**（$\varepsilon$ 从 1 跳到 9，$\mu=1$）：

```matlab
[MM, tau, rho] = temporal_interface_matrix(1, 1, 9, 1)
% ratioElectric = 1/9,  ratioMagnetic = 1/3
% tau = 0.2222,  rho = -0.1111
```

注意到 $\rho$ 是**负的**——时间反射系数可以取负值，别被吓到。

### 2.4 直觉：为什么是 2×2、为什么对称

**为什么正、反两列都出现？** 一个时间界面把入射波**劈成两支**：一支继续沿原方向走（时间透射，频率 $\omega_b\to\omega_a$），一支折返（时间反射，产生"逆向频率"分量）。所以"入射的两个分量 $[E^+;E^-]$"都会同时贡献给"出射的两个分量"，只能是一个 2×2 矩阵——矩阵的每一列对应一个入射方向，每一行对应一个出射方向。

**为什么矩阵对称（反对角元都是 $\rho$）？** 时间界面的介质**空间均匀**，把 $x$ 换成 $-x$（即交换 $+$ 与 $-$）系统完全不变。$+$ 方向入射产生反向分量的比例（$\rho$），必然等于 $-$ 方向入射产生正向分量的比例（$\rho$）。这就是空间各向同性强制出来的对称性。对比一下空间界面：Fresnel 匹配矩阵一般不对称，因为空间界面有"左右介质不同"这种不对称。

**还有个重要性质**：$\det(MM)=\tau^2-\rho^2=r_E r_M\neq 1$。时间调制会对波做功、交换能量，所以时间界面的矩阵**不保模**——这与空间界面散射矩阵的幺正性（保能量）是本质区别。上面例子 $\varepsilon:1\to9$ 时 $\det=1/27<1$；反过来 $\varepsilon:9\to1$ 时 $\det=27>1$，场被放大。这正是"时间开关放大波"的种子。

---

## 三、单值矩阵（monodromy）

### 3.1 对周期时间晶体，直接推进 $[D;B]$

对**周期**时间晶体，与其一个界面一个界面地匹配，不如用更省事的一条路：把状态选成 $[D;B]$，因为在时间界面上 $D,B$ **自动连续**，根本不用界面矩阵，只需要"每一段内推进、段与段直接接上"。

一维 Maxwell（$\partial_x H=-\partial_t D$、$\partial_x E=-\partial_t B$，代入 $E=D/\varepsilon$、$H=B/\mu$，对平面波 $e^{ikx}$ 令 $\partial_x\to ik$）：

$$
\frac{d}{dt}\begin{pmatrix}D\\B\end{pmatrix}
=-ik\begin{pmatrix}
0 & 1/\mu(t)\\
1/\varepsilon(t) & 0
\end{pmatrix}
\begin{pmatrix}D\\B\end{pmatrix}.
$$

验证第一行：$\partial_t D=-\partial_x H=-\partial_x(B/\mu)=-(1/\mu)(ikB)$。✓ 这个演化方程正是 [`temporal_crystal_monodromy.m`](../tmm/temporal_crystal_monodromy.m) 注释里写的那个，手册 §8 也有。

### 3.2 `temporal_crystal_monodromy`：一个周期 $T$ 的整体演化

```matlab
U = temporal_crystal_monodromy(k, epsSequence, muSequence, durations)
```

输入四个东西：固定波矢 `k`、时间上依次经历的介电常数 `epsSequence`、磁导率 `muSequence`（可为标量，自动扩展）、每段时长 `durations`。输出**单值矩阵** $U$：在固定 $k$ 下，$[D;B]$ 状态经过整整一个周期 $T=\sum_m \tau_m$ 后变成 $U[D;B]$。

每一段内 $\varepsilon,\mu$ 常数，生成元 $G_m$ 是常数矩阵：

$$
G_m=-ik\begin{pmatrix}
0 & 1/\mu_m\\
1/\varepsilon_m & 0
\end{pmatrix},
$$

整周期演化就是矩阵指数的乘积，**时间上最早的段在最右边**：

$$
U=\exp\!\bigl(G_M\tau_M\bigr)\,\cdots\,\exp\!\bigl(G_2\tau_2\bigr)\,\exp\!\bigl(G_1\tau_1\bigr).
$$

代码的核心循环（照抄自源码）正是"后乘式累积"：

```matlab
U = eye(2);
for m = 1:M
    generator = -1i*k*[0, 1/muSequence(m); 1/epsSequence(m), 0];
    U = expm(generator*durations(m))*U;   % 新的 U = expm(G_m*tau_m) * 旧的 U
end
```

### 3.3 为什么用 `expm` 而不是 `exp`？

这是新手最常问的。看 [01-matlab-primer.md](01-matlab-primer.md) §2.3 里强调过：`exp(A)` 是**逐元素**的 $e^{A_{ij}}$，`expm(A)` 才是**矩阵指数** $\sum_n A^n/n!$。这里状态是 2 分量向量 $[D;B]$，生成元 $G$ 是 2×2 矩阵，$D$ 和 $B$ 是**耦合**在一起演化的——不能用标量指数把每个分量单独推进。矩阵指数 $e^{G\tau}$ 才是方程 $d\psi/dt=G\psi$ 的真正传播子。这和"二维转动要用转动生成元的矩阵指数"是同一个道理。

> 💡 单段内的物理直觉：$G_m^2=-k^2/(\varepsilon_m\mu_m)\,I$，所以 $G_m$ 的本征值是 $\pm i k/n_m$（纯虚），单段的 $\exp(G_m\tau_m)$ 本征值为 $e^{\pm i k\tau_m/n_m}=e^{\pm i\omega_m\tau_m}$——**只在单位圆上转，不放大**。但注意单段矩阵**不是幺正矩阵**（$\varepsilon\ne\mu$ 时 $G_m$ 不是反 Hermitian），多个这样的矩阵乘起来，乘积的本征值就能离开单位圆——这就是下一节动量带隙的来处。

### 3.4 `det(U)=1`：为什么、怎么验证

因为每个 $\exp(G_m\tau_m)$ 的行列式是 $\exp(\mathrm{tr}(G_m\tau_m))=\exp(0)=1$（$G_m$ 无迹），所以

$$
\det U=\prod_m \det \exp(G_m\tau_m)=1.
$$

`test_smoke.m`（[tests/test_smoke.m](../tests/test_smoke.m)）里就有一句断言：

```matlab
U = temporal_crystal_monodromy(1.7, [1 4], [1 1], [0.4 0.6]);
assert(abs(det(U)-1) < 1e-12, ...
    'Temporal-crystal monodromy determinant is not unity.');
```

自己也能验：

```matlab
U = temporal_crystal_monodromy(1.7, [1 4], [1 1], [0.4 0.6]);
det(U)   % 1.0000
```

$\det U=1$ 意味着两个本征值互为倒数：$\lambda$ 与 $1/\lambda$。这一条是下面能带判据的根基。

---

## 四、能带：从 $U$ 到准频率 $\omega_F$

### 4.1 Floquet 乘法子与准频率

周期系统在时间方向的类比是 Floquet 定理：一整周期演化 $U$ 的本征值 $\lambda$（叫 **Floquet 乘法子**）给出"每周期乘一个 $\lambda$"的增长因子。约定（见 [本文章节 4.1](#41-floquet-乘法子与准频率)）：

$$
\lambda=e^{-i\,\omega_F\,T},
\qquad
\omega_F=\frac{i}{T}\log\lambda.
$$

取**主值对数**，让 $\mathrm{Re}\,\omega_F$ 落在**第一时域 Brillouin 区** $(-\Omega/2,\ \Omega/2]$（$\Omega=2\pi/T$）。这就是 [`temporal_crystal_bands.m`](../tmm/temporal_crystal_bands.m) 做的事：

```matlab
bands = temporal_crystal_bands(kValues, epsSequence, muSequence, durations)
```

对每个 `k`：`U = temporal_crystal_monodromy(...)` → `lam = eig(U)` → `omegaF = 1i*log(lam)/T`，再按 $(\mathrm{Re},\mathrm{Im})$ 排序两条分支（`sortrows`），保证曲线连续可复现。返回的结构体字段：

| 字段 | 形状 | 含义 |
|---|---|---|
| `bands.k` | 1×Nk | 扫描的波矢 |
| `bands.omegaF` | 2×Nk | 两条准频率带（第 1、2 行） |
| `bands.lambda` | 2×Nk | 两个 Floquet 乘法子 |
| `bands.halfTrace` | 1×Nk | $\mathrm{tr}(U)/2$ |
| `bands.T`、`bands.Omega` | 标量 | 周期与调制角频率 |

### 4.2 动量带隙判据：看 `halfTrace`

因为 $\det U=1$，$\lambda$ 满足 $\lambda^2-(\mathrm{tr}\,U)\lambda+1=0$，所以

$$
\mathrm{tr}\,U=2\cos(\omega_F T),
\qquad
\boxed{\ \texttt{halfTrace}=\frac{\mathrm{tr}(U)}{2}=\cos(\omega_F T)\ }.
$$

于是有了整套判据（`halfTrace` 是代码里现成的量，连 $U$ 都不用对角化）：

| 条件 | $\lambda$ | $\omega_F$ | 物理 |
|---|---|---|---|
| $|\texttt{halfTrace}| \le 1$ | $|\lambda|=1$（单位圆上） | 实数 | **允许带**（传播带），场只振荡不放大 |
| $|\texttt{halfTrace}| > 1$ | $|\lambda|\ne 1$（一对 $\lambda,1/\lambda$） | 复数 | **动量带隙**，场随时间指数增长/衰减 |

动量带隙里两个乘法子一个 $|\lambda|>1$（增长）一个 $|\lambda|<1$（衰减），成对出现；对应 $\mathrm{Im}\,\omega_F$ 一正一负、关于 0 对称。这是**时间晶体区别于普通光子晶体的标志性特征**——普通光子晶体的带隙是空间衰逝，这里是**时间上的指数增长**（也叫参数放大）。想确认这些符号约定，翻 [02 章 §4.3](02-physics-background.md#43-动量带隙场随时间指数增长衰减)。

### 4.3 为什么 `halfTrace` 是实数、曲线长什么样

$\texttt{halfTrace}$ 在代码里存成复数，但**算出来虚部恒为零**（这也是 [demo09](../demos/demo09_finite_ptc_order_and_phase.m) 直接画 `real(bands.halfTrace)` 的原因）。原因藏在单段矩阵的结构里：生成元形如 $-ik\begin{pmatrix}0&a\\b&0\end{pmatrix}$（这里 $a=1/\mu$、$b=1/\varepsilon$），其矩阵指数是

$$
\exp\!\Bigl(-ik\tau\begin{pmatrix}0&a\\b&0\end{pmatrix}\Bigr)
=\cos\theta\,I
-i\frac{\sin\theta}{\sqrt{ab}}\begin{pmatrix}0&a\\b&0\end{pmatrix},
\qquad \theta=k\tau\sqrt{ab},
$$

也就是**对角线是相等的实数、反对角线是纯虚数**。这一类矩阵相乘封闭（两两相乘后对角线仍是相等的实数、反对角线仍是纯虚数——你可以拿两个乘一乘验证），所以乘积 $U$ 也长这样：对角线实数，$U$ 的迹是实数。于是 $\texttt{halfTrace}=\cos(\omega_F T)$ 是实量。

曲线长什么样（带内）：

- **允许带**：$\lambda,\lambda^*$ 是单位圆上的共轭对，$\texttt{halfTrace}=\cos(\omega_F T)\in[-1,1]$，随 $k$ 像余弦一样在 $\pm 1$ 之间来回摆动；
- **动量带隙**：曲线**穿破 $\pm 1$**，$\lambda$ 变成实数且互为倒数（一个 $>1$ 一个 $<1$），$\omega_F$ 变复，实部贴在第一时域 BZ 的边界（$\pm\Omega/2$）或中心（0）上。

**画一下你就能看见**：`[1 4]` 二元晶体、占空比 0.5、$\Omega=2\pi$，第一个动量带隙大致落在 $kc/\Omega\approx0.55\text{–}0.75$，带内 `halfTrace` 是好看的余弦型曲线，进带隙处"鼓"出 $\pm 1$。第六节的完整例子就是画这个。

---

## 五、其余 TMM 函数一览

tmm 目录共 12 个函数，除了上面细讲的两个"重头戏"，其余逐个看一眼（都核对过源码，签名见 [docs/tool-reference.md](../docs/tool-reference.md)）：

| 函数 | 一句话作用 | 在哪用 / 依赖 |
|---|---|---|
| `temporal_multilayer_tmm` | 把 $M$ 段"时间界面 + 片内延迟"按时间顺序串成总矩阵 `TM`：`TM = MM_final*DM_M*...*DM_1*MM_first`（列向量约定，最早的矩阵在最右） | [demo04](../demos/demo04_temporal_multilayer_tmm.m) 复现 Ramaccia 四类设计 |
| `temporal_delay_matrix` | 单个时间片内的相位累积 `DM = diag([exp(1i*phi), exp(-1i*phi)])`，其中 `omegaSlab=(nInitial/nSlab)*omegaInitial` | 被 `temporal_multilayer_tmm` 内部调用 |
| `temporal_tmm_spectrum` | 批量扫一串输入频率，逐个调 `temporal_multilayer_tmm`，返回前向/后向电场系数 `forward`/`backward` 和 `detTM` | [demo04](../demos/demo04_temporal_multilayer_tmm.m) 画传递函数 |
| `temporal_finite_crystal_response` | 有限周期时间晶体（`nPeriods` 个周期）夹在相同静态背景之间的**精确**响应：`state = U^nPeriods * initialState`，输入输出都是背景里的方向态 `[E+;E-]` | [demo09](../demos/demo09_finite_ptc_order_and_phase.m) |
| `temporal_domain_wall_mode` | 时间畴壁局域态：在两侧晶体的**共同动量隙**里，把左侧的增长本征态与右侧的衰减本征态匹配，最小化 `abs(det([vGrowth,vDecay]))`，得到局域在时间界面处的态 | [demo10](../demos/demo10_temporal_domain_wall_mode.m) |
| `temporal_db_to_directional` / `temporal_directional_to_db` | `[D;B]` ↔ `[E+;E-]` 互转（见 1.3） | 两条链条之间的"翻译" |
| `temporal_binary_eps_coeff` | 二元时间晶体 $\varepsilon(t)$ 的第 `m` 个 Fourier 系数：$m=0$ 时是时间平均，$m\ne0$ 时是 `(epsA-epsB)*(exp(1i*2*pi*m*dutyA)-1)/(1i*2*pi*m)` | [demo03](../demos/demo03_ptc_bands_pwe_vs_tmm.m)、[demo12](../demos/demo12_ptc_convergence_audit.m) 喂给 PWE 侧 |

几点补充说明：

- **`temporal_delay_matrix` 沿用了 Ramaccia 的 $e^{+i\omega t}$ 时谐约定**（源码注释写明），其余链条用 $e^{ikx-i\omega t}$。比**幅度**不受影响；比**相位**时先统一约定。
- **`temporal_multilayer_tmm` 是"单输入频率"框架**：$k$ 守恒、片内 $\omega_{\text{slab}}=(n_{\text{initial}}/n_{\text{slab}})\,\omega_{\text{initial}}$ 逐段变，`details` 里能取到每段的 `tau`/`rho`/`omegaSlabs`。它和"$k$ 不守恒、$\omega$ 守恒"的空间 TMM 完全对偶。
- **`temporal_domain_wall_mode` 不证明量化的时间 Zak 相**（源码注释明确免责）：它建立的是"体–界面模"，要做拓扑还需要对称性与单胞分析。别把相位差误当拓扑不变量（[demo09](../demos/demo09_finite_ptc_order_and_phase.m) 专门演示了"能带相同、相位不同"）。

---

## 六、完整可运行例子：二元时间晶体的能带与 halfTrace

把前面全串起来：用 `demo03` 的参数（$\varepsilon_A=1.0$、$\varepsilon_B=4.0$、占空比 0.5、$T=1$）算一条能带，再画 `halfTrace` 判据。直接复制到 MATLAB 就能跑。

```matlab
% 05 章示例：二元时间晶体的能带 + halfTrace 动量带隙判据
startup_stm;                     % 把 core/tmm/... 加进路径（v2 系列）

% ---- 1) 二元时间晶体参数（同 demo03）----
epsA = 1.0;  epsB = 4.0;         % 介电常数两个取值
dutyA = 0.5;                     % epsA 占空比
T  = 1;                          % 调制周期
Omega = 2*pi/T;
durations = [dutyA, 1-dutyA]*T;  % 每周期两段时间片：[0.5 0.5]
epsSeq = [epsA epsB];  muSeq = [1 1];

% ---- 2) k 网格：用 kc/Omega（时间晶体的"归一化波矢"，见下面小字）----
kNorm   = linspace(0.02, 1.45, 181);
kValues = kNorm*Omega;

% ---- 3) 能带（精确单值矩阵）----
bands = temporal_crystal_bands(kValues, epsSeq, muSeq, durations);
% bands.omegaF    2×181：两条准频率带
% bands.halfTrace 1×181：trace(U)/2

% ---- 4) 画图：左=能带，右=halfTrace ----
fig = figure('Color','w','Position',[80 80 1100 460]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact');

ax1 = nexttile(tl);
plot(ax1, kNorm, real(bands.omegaF.')./Omega, 'k-', 'LineWidth', 1.6);
hold(ax1,'on');
plot(ax1, kNorm, imag(bands.omegaF.')./Omega, 'r--', 'LineWidth', 1.2);
yline(ax1, 0, 'k:');
xlabel(ax1, 'kc/\Omega');
ylabel(ax1, '\omega_F/\Omega');
title(ax1, 'Floquet 能带（实部黑线，虚部红线）');
legend(ax1, {'Re \omega_F','Im \omega_F'}, 'Location','best');
grid(ax1,'on'); box(ax1,'on');

ax2 = nexttile(tl);
plot(ax2, kNorm, real(bands.halfTrace), 'b-', 'LineWidth', 1.6);
hold(ax2,'on');
yline(ax2, 1, 'r:');  yline(ax2, -1, 'r:');   % ±1 是带隙边界
xlabel(ax2, 'kc/\Omega');
ylabel(ax2, 'Tr(U)/2');
title(ax2, '|Tr(U)/2|>1 即动量带隙');
grid(ax2,'on'); box(ax2,'on');
```

**怎么读这张图**：

- 左图：黑线是 $\mathrm{Re}\,\omega_F/\Omega$ 的两条分支，在所扫的 $kc/\Omega\in[0.02,1.45]$ 范围内像普通能带一样散开又收拢；红虚线是 $\mathrm{Im}\,\omega_F/\Omega$，**只在 $kc/\Omega\approx0.55\text{–}0.75$ 一带非零**——那里就是动量带隙，场随时间指数增长/衰减。
- 右图：`halfTrace` 曲线在带内于 $[-1,1]$ 间振荡，进带隙处**穿破 $\pm 1$**（红虚线）。`halfTrace` 穿过 $\pm 1$ 的位置，正是左图红虚线长出、黑线"钉住"的位置。两条判据完全一致。
- 想确认自己没算错：带隙内两个 $|\lambda|$ 一定是一个 $>1$、一个 $<1$（互为倒数），$\mathrm{Im}\,\omega_F$ 关于 0 对称。

> 💡 为什么横轴用 $kc/\Omega$ 而不是 $\bar{k}=k\Lambda/(2\pi)$？时间晶体没有空间周期，没有 $\Lambda$；它的"归一化波矢"天然是 $kc/\Omega$（$c_0=1$ 时即 $k/\Omega=kT/(2\pi)$），是时间方向的 Brillouin 区类比。ST-PWE 那边的 $\bar{k}=k\Lambda/(2\pi)$（见 [04-st-pwe.md](04-st-pwe.md)）是用空间周期归一化的——两套轴对应不同的物理设置，别混。

---

## 七、适用边界与坑

最后把新手最容易踩的坑集中讲一遍：

### 7.1 只对"时间分段常数"严格精确

单值矩阵用 $\exp(G_m\tau_m)$ 精确演化每一段，**前提是每段内 $\varepsilon,\mu$ 是常数**。对连续时间调制（比如 $\varepsilon(t)=\varepsilon_c(1+\delta\sin\Omega t)$），只能把它近似成很多个台阶：台阶越密、近似越好、但段数 $M$ 越大。所以对连续调制，"时间切片数"就是收敛参数——算完要像 01 章教的那样做**收敛性检查**。

### 7.2 空间均匀是根本前提

整条 TMM 链靠的是**空间均匀 ⇒ $k$ 守恒**。一旦 $\varepsilon$ 也随空间变化（时空晶体，$\varepsilon(x,t)$），$k$ 不再是好量子数，这整套 2×2 机制**直接失效**——那是 ST-PWE 的管辖范围（[04-st-pwe.md](04-st-pwe.md)）。反过来，PWE 用有限 Fourier 截断逼近突变时间层也不是"精确解"。三套引擎的适用范围见 [03-methods-and-map.md](03-methods-and-map.md)。

### 7.3 $D,B$ 连续只是"本构切换"模型

`temporal_interface_matrix` 默认 $D,B$ 连续。真实开关电路若保 $E$、注入电荷或空间局域开关，就得改用 `temporal_interface_matrix_jump` 的 `'EB'`/`'DH'`/`'EH'` 或自定义 `jumpD`/`jumpB`。这一条 [demo11](../demos/demo11_coherent_time_interface.m) 里 DB 与 EB 的相干相消点不同，是现成的反面教材。

### 7.4 和 ST-PWE 互相验证（demo03）

二元方波不连续，Fourier 系数只按 $1/|m|$ 衰减、且 50% 占空比时偶次谐波精确为 0，所以时间 PWE 需要很大 `Mtime`（demo03 用 19）才收敛。**TMM 是这一段上的精确参照**（无截断），[demo03](../demos/demo03_ptc_bands_pwe_vs_tmm.m) 把 PWE 散点画在 TMM 线上验证收敛，[demo12](../demos/demo12_ptc_convergence_audit.m) 定量扫 `Mtime`。修改 PWE 代码后，跑 demo03 是最重要的交叉验证。

### 7.5 其他约定与符号坑

- **`temporal_delay_matrix` 用 $e^{+i\omega t}$ 约定**（Ramaccia 论文），TMM 的 D/B 链和 FDTD/ST-PWE 用 $e^{ikx-i\omega t}$。比相位时先统一。
- **`halfTrace` 实数的结论依赖实、正、无损的 $\varepsilon,\mu$**。有损/增益材料（复 $\varepsilon$）时单段矩阵不再属于"实对角 + 纯虚反对角"类，`halfTrace` 会变复，判据要重新讨论。
- **时间调制会与场交换能量**，总场能量一般不守恒：界面矩阵 $\det(MM)\ne1$，带隙内 $|\lambda|\ne1$。这与"静态介质的能量守恒"完全不同，别拿错守恒律当检查（FDTD 章节 [06-fdtd.md](06-fdtd.md) 会再讲）。
- `tmm` 只给"本征乘法子/准频率"，不给场形；想看动量隙内场的**实际时空演化**和主导谐波，要回到 ST-PWE 的 `stpwe_reconstruct_field`（[04-st-pwe.md](04-st-pwe.md)）或 FDTD（[06-fdtd.md](06-fdtd.md)）。

---

## 下一步

- 想从 Maxwell 方程另一条路（Fourier 展开 → 本征值问题）交叉验证这里的能带 → [04-st-pwe.md](04-st-pwe.md)
- 想看时间界面/时间晶体/时空晶体的物理图景 → [02-physics-background.md](02-physics-background.md)
- 想用 FDTD 在时域"实况转播"动量带隙里的增长 → [06-fdtd.md](06-fdtd.md)
- 想逐个跑并读懂 12 个 demo → [08-demos-walkthrough.md](08-demos-walkthrough.md)
- 遇到报错、参数对不上 → [10-faq-troubleshooting.md](10-faq-troubleshooting.md)
