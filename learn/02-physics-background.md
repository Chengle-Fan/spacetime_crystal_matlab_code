# 02 · 时变介质物理背景：时间界面、光子时间晶体与时空晶体

> 目标读者：**物理基础 OK、还没开始看代码**的读者。
> 这一章一行数值方法都不深入，只做一件事——给你装上一套**物理直觉**：当"周期性"从空间搬到时间，电磁波会多出哪些不寻常的行为。这些直觉是后面 [03-methods-and-map.md](03-methods-and-map.md) 里三套数值引擎（ST-PWE、时间 TMM、FDTD）各自要算的东西。

读完这一章，你应该能回答三个问题：
1. 什么是"时间界面"？为什么一个波经过它会被拆成**前向 + 后向**两列？
2. 什么是"动量带隙"？为什么它里面的场会**随时间指数增长**，而不是像普通光子晶体那样在空间里衰减？
3. "光子时间晶体"和"时空晶体"差在哪？后面代码里的双 Fourier 展开 $\varepsilon_{m,n}$ 是从哪冒出来的？

全程不需要你运行任何代码；本章只出现几段"示意性质"的 MATLAB 小片段，帮你把物理和代码对上号。

---

## 一、为什么研究时变介质：把"空间周期"换成"时间周期"

### 1.1 一分钟回顾：空间光子晶体

你已经熟悉普通光子晶体：介电常数在空间上周期变化

$$
\varepsilon(x+\Lambda) = \varepsilon(x),
$$

$\Lambda$ 是晶格常数。Bloch 定理告诉我们**波矢 $k$ 是好量子数**（守恒量），于是可以定义一个函数 $\omega(k)$——能带结构。能带之间夹着**频率带隙（frequency bandgap）**：

```
   ω
   │          ╱╲
   │         ╱  ╲         ← 第二条能带
   │        ╱    ╲
   │ ──────╱──────╲─────  ← 频率带隙：这个 ω 区间没有传播解
   │       ╱      ╲
   │      ╱        ╲
   │   ╱╱╱╱         ╲╲    ← 第一条能带
   └────────────────────▶  k
    -g/2             g/2   （第一 Brillouin 区）
```

带隙的物理来源是 **Bragg 散射**：空间周期结构把 $k$ 和 $k-g$ 两个平面波分量"粘"在一起，在 Brillouin 区边界处原本相交的两条线劈开，中间空出来的频率区间就是带隙。带隙里没有传播模——你想在这个频率区间里传波，$k$ 会变成复数（渐逝波），场在**空间**上指数衰减。

### 1.2 直觉实验：把"周期"从空间搬到时间

现在做一件看起来很小、后果很不同的事：**把下标 $x$ 换成 $t$**。

$$
\varepsilon(t+T) = \varepsilon(t),
$$

介质在空间上完全均匀，但介电常数在时间上以周期 $T$ 起伏。直觉上会发生什么？

- **空间平移不变性还在**：介质在每个瞬间都是空间均匀的，所以 $k$ 依然是严格的好量子数。
- **时间平移不变性没了**：介质参数本身随时间变，一个"时间平移"前后的系统不是同一个系统，所以**频率 $\omega$ 不再是守恒量**。

于是能带结构"横纵轴对调"了：对静态晶体，你习惯"固定 $k$ 求 $\omega(k)$"；对时间晶体，更自然的是"固定 $\omega$ 求 $k(\omega)$"，而且 $\omega$ 作为准频率（quasifrequency）可以折叠进"时间 Brillouin 区"。这个"$\omega$ 和 $k$ 谁是好量子数"的对调，是本篇所有现象的根源。建议把这句话记下来：

> **空间周期 → 破坏空间平移 → 作用在 $k$ 上；时间周期 → 破坏时间平移 → 作用在 $\omega$ 上。**

### 1.3 三列对照表

| | **空间光子晶体** | **光子时间晶体（PTC）** | **时空晶体（spacetime crystal）** |
|---|---|---|---|
| 调制对象 | $\varepsilon(x+\Lambda)=\varepsilon(x)$ | $\varepsilon(t+T)=\varepsilon(t)$ | $\varepsilon(x+\Lambda,\,t+T)=\varepsilon(x,t)$ |
| 倒空间 | $k$ 是好量子数，第一 Brillouin 区 $k\in[-g/2,g/2]$ | $\omega$ 折叠进"时间 Brillouin 区" $\omega\in(-\Omega/2,\Omega/2]$，$k$ 守恒 | 二维倒空间 $(k,\omega)$，两者各有周期性 |
| 带隙类型 | **频率带隙**：固定 $\omega$ 无实 $k$，$k$ 变复 | **动量带隙**：固定 $k$ 无实 $\omega$，$\omega$ 变复 | 两类带隙**可以共存** |
| 场行为 | 带隙内场**空间**衰逝 | 带隙内场**时间**指数增长/衰减 | 取决于落在哪类带隙 |

> 关于"固定 $k$ 求 $\omega$"和"固定 $\omega$ 求 $k$"这两个本征值问题、以及它们怎么互相对照，工具箱里有一个专门演示（`demos/demo02_complex_frequency_and_momentum_gaps.m`），我们在 [04-st-pwe.md](04-st-pwe.md) 会细讲。这里先建立印象：**时间周期介质里，$\omega$ 和 $k$ 不再是一一对应的。**

---

## 二、最小单元：时间界面（time interface）

### 2.1 场景：$\varepsilon$ 在 $t=0$ 跳变

空间界面你已经很熟：在某个位置 $x=x_0$，折射率突变，入射波被劈成反射波 + 透射波。现在做一个完全对称的操作，把"位置"换成"时刻"：

**整个空间的 $\varepsilon$ 在 $t=0$ 那一瞬间，从 $\varepsilon_1$ 跳到 $\varepsilon_2$，其他时候都是均匀的。**

比如 `demos/demo05_fdtd_temporal_interface.m` 里就是 $\varepsilon$ 从 $1.5^2$ 跳变到 $2.5^2$。这就是**时间界面（time interface / temporal boundary）**。它是对波场施加的一个"全局、瞬时的踢"。

### 2.2 两件事：$k$ 守恒，$\omega$ 改变

时间界面发生的物理可以浓缩成两条：

1. **波矢 $k$ 守恒。** 跳变是"对整片空间同时"发生的，空间平移不变性没有被破坏，所以每个平面波分量 $e^{ikx}$ 继续是同一个 $e^{ikx}$——$k$ 在界面前后取同一个值。
2. **频率 $\omega$ 改变。** 跳变瞬间破坏了时间平移对称性，频率不再守恒。新频率由色散关系决定：归一化单位下 $c_0=1$，$\omega_i = k/n_i$（$n_i=\sqrt{\varepsilon_i\mu_i}$ 是折射率），所以

$$
\omega_2 = \frac{n_1}{n_2}\,\omega_1 .
$$

直觉：$\varepsilon$ 变大（折射率变大）→ 波速变慢 → **频率降低**；反过来 $\varepsilon$ 变小 → 频率升高。这个"$\omega$ 变了、$k$ 没变"的效果，正是 demo05 时空图里波包过界面后**波长不变、但波动快慢改变**的原因。

### 2.3 与空间界面的对照

| | **空间界面（$x=x_0$）** | **时间界面（$t=t_0$）** |
|---|---|---|
| 守恒量 | $\omega$（时间平移不变） | $k$（空间平移不变） |
| 改变量 | $k$（两侧折射率不同） | $\omega$（由色散 $\omega=k/n$ 决定） |
| 几何 | 空间上分成两个区域 | 时间上分成两个阶段 |
| 新出现的波 | 反射波（空间反向）+ 透射波 | **前向波 + 后向波**（见下） |

这张表就是 1.2 节那句话的量化版：空间界面守恒的是时间量，时间界面守恒的是空间量。**物理是对称的，只是把 $x$ 和 $t$ 的角色对调。**

### 2.4 时间反射 / 时间折射

关键问题：界面前只有一列向右传的波（$E_+$），界面后场必须满足"$k$ 相同但 $\omega$ 不同"。一个频率 $\omega_2$ 的平面波要同时满足"向右、向左都能以这个频率传播"，于是——**界面后同时出现两列波**：

```
    t
    │
    │    ↖ 时间反射（后向波，频率 ω₂）
    │    ↗ 时间透射（前向波，频率 ω₂）
    │
    │═══════════════════════════  t = 0  时间界面：ε₁ → ε₂
    │
    │    →  入射波（频率 ω₁）
    │
    └──────────────────────────────────▶  x
```

- **时间透射（temporal refraction）**：继续向前传的波，频率变成 $\omega_2$；
- **时间反射（temporal reflection）**：新冒出来的向后传的波，频率也是 $\omega_2$。

两列波**频率相同、空间方向相反**——这就是"波被时间劈开"的意思。`demo05_fdtd_temporal_interface.m` 的头部注释把这件事说得很直白：*"Wave splitting at one time boundary"*（在单个时间界面上的波劈裂），它用一个复高斯波包，把界面后的场干净地分解成前向/后向两列，再和下面要讲的解析系数对比。

> ⚠️ 一个常见的困惑：时间反射听起来像"波在时间上反弹"。更准确的图像是——**入射波在 $t=0$ 之后被迫"换挡"到新频率**，而麦克斯韦方程要求这个新频率的波同时存在两个传播方向的分量。你不必为"凭什么冒出后向波"纠结，下一节用界面条件一推就清楚了。

---

## 三、Morgenthaler 系数：时间界面的"透射 / 反射"

### 3.1 $D$、$B$ 连续（DB 模型）

现在把 2.4 节变成公式。记界面两侧的介电常数为 $\varepsilon_1,\varepsilon_2$、磁导率为 $\mu_1,\mu_2$，折射率 $n_i=\sqrt{\varepsilon_i\mu_i}$。因为 $k$ 守恒，把场写成两列反向波（下标 1、2 表示界面前、后）：

$$
E(x,t) = E_{1+}\,e^{i(kx-\omega_1 t)} + E_{1-}\,e^{i(kx+\omega_1 t)} \quad (t<0),
$$

$$
E(x,t) = E_{2+}\,e^{i(kx-\omega_2 t)} + E_{2-}\,e^{i(kx+\omega_2 t)} \quad (t>0).
$$

在**空间均匀**介质里，$D=\varepsilon E$，$B=\mu H$；再利用前向波 $H=E/\eta$、后向波 $H=-E/\eta$（$\eta=\sqrt{\mu/\varepsilon}$ 是波阻抗），可得方向态与 $(D,B)$ 的关系：

$$
D = \varepsilon(E_+ + E_-), \qquad B = n\,(E_+ - E_-), \qquad n=\sqrt{\varepsilon\mu}.
$$

在时间界面上，麦克斯韦方程里出现的是 $\partial D/\partial t$ 和 $\partial B/\partial t$，它们在跳变瞬间保持有限，所以 **$D$ 与 $B$ 连续**（这称为 **DB 模型**，本构切换模型 constitutive-switching model）：

$$
D_1 = D_2,\qquad B_1 = B_2.
$$

代入方向态：$\varepsilon_1(E_{1+}+E_{1-}) = \varepsilon_2(E_{2+}+E_{2-})$，$n_1(E_{1+}-E_{1-}) = n_2(E_{2+}-E_{2-})$。把后一组写成 $E_{2\pm} = \tfrac12(D_2/\varepsilon_2 \pm B_2/n_2)$，用 $D_2=D_1$、$B_2=B_1$ 代回去，整理成矩阵形式：

$$
\begin{pmatrix} E_{2+} \\ E_{2-} \end{pmatrix}
= \begin{pmatrix} \tau & \rho \\ \rho & \tau \end{pmatrix}
\begin{pmatrix} E_{1+} \\ E_{1-} \end{pmatrix},
\qquad
\tau = \tfrac12\Big(\tfrac{\varepsilon_1}{\varepsilon_2} + \tfrac{n_1}{n_2}\Big), \quad
\rho = \tfrac12\Big(\tfrac{\varepsilon_1}{\varepsilon_2} - \tfrac{n_1}{n_2}\Big).
$$

- $\tau$ 叫**时间透射系数**：入射的前向波有多少继续向前（$\rho=0$ 时全部透射，若 $\varepsilon_1=\varepsilon_2$ 且 $n_1=n_2$ 则退化为恒等映射）；
- $\rho$ 叫**时间反射系数**：有多少变成后向波。注意这里"反射"不是空间的镜面反射，而是"时间界面生成的一列反向波"。

### 3.2 代码里的 $\tau$、$\rho$：照代码写

工具箱把上面两步拆开：先规定"哪些量连续"，再组装矩阵。核心实现在 `tmm/temporal_interface_matrix_jump.m` 里，这几行就是**全文的骨架**：

```matlab
% D_after = jumpD*D_before,  B_after = jumpB*B_before
ratioElectric  = jumpD * epsBefore / epsAfter;
ratioMagnetic  = jumpB * sqrt((epsBefore*muBefore)/(epsAfter*muAfter));
tau = 0.5*(ratioElectric + ratioMagnetic);
rho = 0.5*(ratioElectric - ratioMagnetic);
MM  = [tau, rho; rho, tau];      % after = MM*before，作用在 [E+; E-] 列向量上
```

对照公式：`ratioElectric` 就是 $\varepsilon_1/\varepsilon_2$ 乘上"电跳变因子"`jumpD`，`ratioMagnetic` 就是 $n_1/n_2 = \sqrt{\varepsilon_1\mu_1/(\varepsilon_2\mu_2)}$ 乘上"磁跳变因子"`jumpB`。DB 模型下 `jumpD = 1`、`jumpB = 1`，于是 `ratioElectric = epsBefore/epsAfter`、`ratioMagnetic = sqrt((epsBefore*muBefore)/(epsAfter*muAfter))`，恰好回到 3.1 节的 $\tau=\tfrac12(\varepsilon_1/\varepsilon_2 + n_1/n_2)$。

平时直接调用便捷封装 `tmm/temporal_interface_matrix.m`（它内部就是 DB 模型）：

```matlab
[MM, tau, rho] = temporal_interface_matrix(epsBefore, 1, epsAfter, 1);
% demo05 就是这样调用的：temporal_interface_matrix(nBefore^2, 1, nAfter^2, 1)
```

> 注意到 `det(MM) = tau^2 - rho^2 = ratioElectric * ratioMagnetic`，一般不等于 1——**时间界面不像空间界面的散射矩阵那样保模（幺正）**。这个"不保模"就是能量可以注入/抽出的起点，第五节再展开。

### 3.3 "哪些量连续"不是唯一选择

3.1 节的 DB 模型是一个**物理假设**，不是普适口号。真实的开关电路可能保持的是电压、断开的是电荷、或由控制电路注入冲量——每种机制对应不同的"跳变定律"。所以工具箱在 `tmm/temporal_interface_matrix_jump.m` 里提供了四种命名模型和一种自定义方式，文档串和 [05 章 §2.3](05-temporal-tmm.md)（temporal_interface_matrix_jump 的跳变律）都写明了这一点：

| 模型 | 连续的物理量 | 物理图像 |
|---|---|---|
| `'DB'`（默认） | $D$、$B$ | 标准"本构切换"体模型 |
| `'EB'` | $E$、$B$ | 保持电场，$\varepsilon$ 变了 $D$ 跟着跳 |
| `'DH'` | $D$、$H$ | 保持电位移和磁场强度 |
| `'EH'` | $E$、$H$ | 场量连续，本构量跳变 |
| `struct('jumpD',qD,'jumpB',qB)` | 自定义 | 你的微观开关模型 |

（`jumpD`/`jumpB` 的具体取值属于代码层面，[05-temporal-tmm.md](05-temporal-tmm.md) 有完整表格，这里只讲物理含义。）

`demo11_coherent_time_interface.m` 就用 `'DB'` 和 `'EB'` 各算一遍，展示同一材料跳变、不同跳变定律会给出不同的 $\tau,\rho$，连相干相消点都不一样。

> **这是物理问题，不是数值问题。** 选哪种模型取决于你的实验装置（你的开关到底连续了什么），代码只负责把你选的定律忠实地算出来。研究前必读本章第 3.3 节。要提醒的是：若装置有**色散、内部电路状态、有限开关时间或空间局域开关**，两种标量跳变因子就不够用了，需要扩大状态空间重新积分——工具箱文档把这层边界写得清清楚楚。

---

## 四、时间多层与光子时间晶体：动量带隙

### 4.1 从"一次跳变"到"一串跳变"

单个时间界面是最小单元。现在让 $\varepsilon$ 周期性地反复跳变，比如一个周期 $T$ 内先 $\varepsilon_A$ 持续 $T/2$、再 $\varepsilon_B$ 持续 $T/2$，然后重复。这就像把一堆时间界面排成一串，形成**时间多层（temporal multilayer）**；当层数无限、模式严格周期 $\varepsilon(t+T)=\varepsilon(t)$ 时，就是**光子时间晶体（photonic time crystal, PTC）**——1.2 节那张对照表的中间一列。

TMM 侧处理它的方式是"时间传输矩阵"：每个时间片内 $(\varepsilon,\mu)$ 恒定，用一个 $2\times2$ 矩阵演化 $(D,B)$ 态，整周期乘起来得到**单值矩阵（monodromy）** $U$。`tmm/temporal_crystal_monodromy.m` 的头部注释把这个方程的雏形写了出来：

$$
\frac{d}{dt}\begin{pmatrix} D \\ B \end{pmatrix}
= -i k \begin{pmatrix} 0 & 1/\mu \\ 1/\varepsilon & 0 \end{pmatrix}
\begin{pmatrix} D \\ B \end{pmatrix},
$$

每个时间界面上 $D,B$ 连续，所以在这个表示里**不需要额外的界面矩阵**——连续量自己保证了演化无缝衔接。

### 4.2 周期调制与 Floquet：一段直觉

对时间周期系统，解不是简单的 $e^{-i\omega t}$，而是 Floquet 理论说的"准周期基函数"的组合：调制频率 $\Omega=2\pi/T$ 会把频率 $\omega$ 与 $\omega\pm\Omega,\ \omega\pm2\Omega,\dots$ 的边带耦合在一起。数学上，**$\omega$ 变成了准频率（quasifrequency）**，只定义到模 $\Omega$（折进"时间 Brillouin 区" $\omega\in(-\Omega/2,\Omega/2]$）。

数值上怎么做？`tmm/temporal_crystal_bands.m` 对每个 $k$：算单值矩阵 $U$ → 求它的特征值 $\lambda$（Floquet 乘子）→ 由

$$
\lambda = e^{-i\omega_F T} \quad\Longrightarrow\quad \omega_F = \frac{i\log\lambda}{T}
$$

反解出准频率 $\omega_F$（主值对数把 $\mathrm{Re}\,\omega_F$ 放进第一时间 Brillouin 区）。它返回的 `bands` 结构体里有 `.omegaF`（2×Nk 的两条带）、`.lambda`、`.halfTrace` 等字段。

### 4.3 动量带隙：场随时间指数增长/衰减

$U$ 是 $2\times2$ 且 $\det U=1$（每个时间片生成元无迹），所以特征值成对 $\lambda,1/\lambda$，满足 $\lambda + \lambda^{-1} = \mathrm{Tr}\,U$。于是能不能传播，由一个标量判别式决定：

$$
h(k) = \frac{\mathrm{Tr}\,U(k)}{2}.
$$

- $|h(k)|\le 1$：$\lambda$ 在单位圆上，$\omega_F$ 为实，是**传播带**；
- $|h(k)|>1$：$\lambda$ 不在单位圆上，$\omega_F$ 取**复值**，这个 $k$ 区间就是**动量带隙（momentum gap）**。

```
   h(k) = Tr(U)/2           ← 半迹判别式（代码返回 bands.halfTrace）
  1.5 ┤          ┌────┐
      ┤       ┌──┘    └──┐         动量带隙：|h| > 1
  1.0 ┤═══════╪══════════╪═══════  ← h = +1 判据线
      ┤     传播带 │       │  传播带
  0.0 ┤─────────────────────────
      ┤        │         │
 -1.0 ┤═══════╪══════════╪═══════  ← h = −1 判据线
  1.5 ┤      └──┐    ┌──┘
      ┤         └────┘
      └─────────────────────────▶  k
       ↑                 ↑
     动量带隙          动量带隙
```

动量带隙里 $\lambda$ 与 $1/\lambda$ 一个 $|\lambda|>1$、一个 $|\lambda|<1$，对应场在时间上**指数增长**和**指数衰减**两支（互为共轭的复准频率 $\pm i\,\mathrm{Im}\,\omega_F$）。这段能跑的最小代码（`tmm/temporal_crystal_bands.m` 的用法，演示脚本 `demos/demo03_ptc_bands_pwe_vs_tmm.m` 与 `demos/demo09_finite_ptc_order_and_phase.m` 都这么调）：

```matlab
% 二元时间晶体：ε 在 1 与 4 之间按占空比 0.5 跳变
epsA = 1;  epsB = 4;
T = 1;
epsSequence = [epsA, epsB];    % 时间上先 εA 再 εB
durations   = [0.5, 0.5]*T;    % 每个时间片持续半个周期

kNorm   = linspace(0.05, 1.4, 200);
kValues = kNorm * 2*pi/T;      % 物理波矢（TMM 用的是物理 k）

U = temporal_crystal_monodromy(kValues(50), epsSequence, [1 1], durations);
halfTrace = trace(U)/2;        % |halfTrace| > 1 ⇒ 这个 k 在动量带隙里

bands = temporal_crystal_bands(kValues, epsSequence, [1 1], durations);
% bands.omegaF   : 2×Nk 准频率（两条带）
% bands.lambda   : 2×Nk Floquet 乘子 λ
% bands.halfTrace: 1×Nk 带隙判别式 h(k)
```

### 4.4 一张最重要的对照

把两种带隙并排看（这是 `demos/demo02_complex_frequency_and_momentum_gaps.m` 的主题，也是本章反复强调的"$\omega$ 与 $k$ 角色对调"）：

| | **频率带隙（空间晶体）** | **动量带隙（PTC）** |
|---|---|---|
| 好量子数 | $\omega$（时间平移不变） | $k$（空间平移不变） |
| 在带隙里变复的量 | 固定实 $\omega$ → $k$ 变复 | 固定实 $k$ → $\omega_F$ 变复 |
| 场的表现 | **空间**衰逝（渐逝波） | **时间**指数增长 / 衰减 |
| 能量 | 反射回去，不积累 | 调制向场注入能量（放大） |

一句话记住：**静态晶体的带隙是"频率带隙"，时间晶体的带隙是"动量带隙"**——前者场在空间里出不去，后者场在时间里越演越大。这正是 PTC 最标志性的现象：**参数放大**。

---

## 五、参量放大与能量：调制就是泵浦

### 5.1 能量不守恒（但物理没被打破）

时间界面和 PTC 里，**电磁场单独的能量不守恒**。原因很简单：介质参数由外部驱动按 $\varepsilon(t)$ 变化，外部驱动就是对场做功的"泵浦"。`demo05_fdtd_temporal_interface.m` 的第三张子图就是瞬时能量曲线——在 $t=t_{\mathrm{switch}}$ 处能量有一个跳变，那就是调制瞬间注入/抽出的能量。`fdtd1d_db` 的文档串和本章第 5.1 节都直接写明：*"时间调制会与场交换能量，所以总场能量一般不守恒"*。

守恒的到底是什么？是"场 + 泵浦"这个更大的系统。我们数值上只看场，所以能量这条曲线会跳、会涨——**这不是数值错误，是物理**。

### 5.2 直觉：两个时间谐波被调制"粘"在一起

放大是怎么来的？给一个不引公式的直觉：

1. 时间调制以频率 $\Omega$ 起伏，它扮演的角色是一个**媒介**，能把频率为 $\omega$ 的波和一个"伙伴"——频率为 $\omega\pm\Omega$ 的波——**耦合**起来（这就是 4.2 节说的 Floquet 边带）。
2. 当两个波与调制满足某种"时间相位匹配"（波之间 + 调制的频率和相位对得齐）时，耦合最强，一个波**持续被放大**。
3. 放大对的正反两个分量就是 4.3 节的 $|\lambda|>1$ 增长支和 $|\lambda|<1$ 衰减支：**$|\lambda|>1$ 的支是"被调制泵上能量"的那一支**，它们是复准频率 $\pm i\,\mathrm{Im}\,\omega_F$ 的共轭对。
4. 光子图像里，这相当于"泵浦光子 $\Omega$ 劈成一对场光子"，两个场光子的频率加起来等于泵浦频率——**光子对产生**；反过来，两列波也可以合起来把能量还给泵浦（减幅）。

> 细节（耦合强度、相位匹配的具体条件、放大率随调制深度的变化）交给后面 [05-temporal-tmm.md](05-temporal-tmm.md) 和 [04-st-pwe.md](04-st-pwe.md) 里的公式。这里只要记住一句：**时间调制是主动介质，它既可能给波充电，也可能把波的能量抽走——这就是 PTC 与静态光子晶体最大的不同。**

---

## 六、时空晶体：空间与时间同时周期

### 6.1 定义与直觉

把第 1 节的两张表合并，就得到本工具箱的核心研究对象——**时空晶体（spacetime crystal）**：

$$
\varepsilon(x+\Lambda,\,t+T) = \varepsilon(x,t).
$$

介质同时具有空间周期 $\Lambda$ 和时间周期 $T$。它既"在空间上像光子晶体"（Bragg 散射），又"在时间上像 PTC"（Floquet 耦合、泵浦放大）。一个具体的例子（`demos/demo01_reproduce_fig2_stpwe.m` 复现的 Park–Min 论文介质）：单个元胞 $\Lambda$ 内前 $3/4$ 是静态介质，后 $1/4$ 被 $\varepsilon=\varepsilon_c(1+\delta\sin\Omega t)$ 时间调制。

### 6.2 二维能带结构

因为 $k$ 和 $\omega$ 都是"部分好的"量子数（各有一重周期性），能带不再是一条一维曲线 $\omega(k)$，而是 $(k,\omega)$ 平面上的**二维结构**：

- **固定 $k$ 求 $\omega$**：得到 Floquet 能带，频率带隙处出现复 $\omega$；
- **固定 $\omega$ 求 $k$**：得到复动量能带，动量带隙处出现复 $k$。

两种带隙可以在 $(k,\omega)$ 平面**共存**：落在某个 $k$ 区间的复 $\omega$ 是动量带隙（时间放大），落在某个 $\omega$ 区间的复 $k$ 是频率带隙（空间渐逝）。把两种扫描画在同一张图上，重合处是传播模，分歧处就是各种带隙——这正是 `examples/example_dual_sweep_bands.m`（[example_dual_sweep_bands_讲解.md](../examples/example_dual_sweep_bands_讲解.md)）做的事。

### 6.3 双 Fourier 展开：ST-PWE 的种子

时空晶体是双周期的，最自然的数学工具就是**双 Fourier 展开**：把 $\varepsilon$ 同时按空间倒格矢 $g=2\pi/\Lambda$ 和时间频率 $\Omega=2\pi/T$ 展开：

$$
\varepsilon(x,t) = \sum_{m,n} \varepsilon_{m,n}\; e^{\,i\,n\,g\,x - i\,m\,\Omega\,t}.
$$

- $n$ 是**空间谐波**指标（$n g$ 是倒格矢的整数倍，来自空间周期性）；
- $m$ 是**时间谐波**指标（$m\Omega$ 是调制频率的整数倍，来自时间周期性）。

系数 $\varepsilon_{m,n}$ 把整个时空介质浓缩成一串数。比如纯正弦调制 $\varepsilon(t)=\varepsilon_c(1+\delta\sin\Omega t)$ 只有 $m=0,\pm1$ 三个非零时间谐波；而方波跳变的系数按 $1/|m|$ 衰减（收敛慢，需要更多谐波）。

**这一行展开式就是后面 ST-PWE（时空平面波展开法）的全部数学起点**：`core/stpwe_build_system.m` 用 $\varepsilon_{m,n}$ 组装卷积矩阵，`core/stpwe_solve_omega.m` 把麦克斯韦方程变成"固定 $k$ 的广义本征值问题 $A\Phi=\omega B\Phi$"。你在本章先种下"$\varepsilon_{m,n}$ 是介质的身份证"这个念头，等读到 [04-st-pwe.md](04-st-pwe.md) 时会发现一切从它长出来。

---

## 七、收尾预告：三套引擎各管哪一段

物理图像齐了，该看数值工具了。工具箱的三套"引擎"恰好从三个角度研究上面这些现象，互相验证：

| 引擎 | 目录 | 研究什么 | 对应本章哪一节 |
|---|---|---|---|
| **ST-PWE**（时空平面波展开） | `core/` | 把场和介质双 Fourier 展开，化成本征值问题，算能带、复频率/复波数、场型 | 6.3 的展开式 → 时空晶体的二维能带与两类带隙 |
| **时间 TMM**（时间传输矩阵 / 单值矩阵） | `tmm/` | 精确演化和匹配时间界面与时间多层，由 `temporal_crystal_bands` 出精确能带与动量带隙 | 3 的时间界面系数、4 的动量带隙、4.1 的单值矩阵 |
| **FDTD**（时域有限差分） | `fdtd/` | 直接在时空网格上推进麦克斯韦方程，看波包如何演化 | 2 的时间界面波劈裂、5 的能量交换 |

它们怎么配合、代码放在哪、`v1`/`v2` 两套入口怎么分辨，是 [03-methods-and-map.md](03-methods-and-map.md) 的内容。想动手感受"波被时间劈开"，可以先去跑 `demos/demo05_fdtd_temporal_interface.m`；想马上看到动量带隙的指数增长，跑 `demos/demo03_ptc_bands_pwe_vs_tmm.m` 或 `demos/demo09_finite_ptc_order_and_phase.m`；想见识"时间畴壁"里的局域态（两列顺序相反的 PTC 拼接，在共同动量带隙里匹配增长支与衰减支），跑 `demos/demo10_temporal_domain_wall_mode.m`——不过更系统的路线还是按 [03-methods-and-map.md](03-methods-and-map.md) 一章章来。

---

## 下一步

- 先建立工具箱整体地图，知道每个函数住在哪个目录 → [03-methods-and-map.md](03-methods-and-map.md)
- 想复习 MATLAB 语法 → 回翻 [01-matlab-primer.md](01-matlab-primer.md)
- 想核对"时间界面连续量"等物理约定的完整清单 → [05 章 §2](05-temporal-tmm.md)（时间界面匹配矩阵 / 界面连续量），归一化约定见 [04 章 §1.1](04-st-pwe.md)
