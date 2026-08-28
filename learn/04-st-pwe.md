# 04 · ST-PWE 时空平面波展开精讲

> 目标读者：**会麦克斯韦方程、能带、Fourier 级数，但还没写过几行 MATLAB** 的物理学生。
> 这一章我们做一件"看起来很神奇、其实很机械"的事：把"解偏微分方程"变成"解一个矩阵本征值问题"，然后让 MATLAB 替我们把能带算出来。
>
> 读完这一章，你应该能做到：
>
> 1. 从一维 Maxwell 方程出发，亲手推出 ST-PWE 的广义本征值问题 $A(k)\Phi=\omega B\Phi$；
> 2. 说清"固定 $k$ 求 $\omega$"和"固定 $\omega$ 求 $k$"两个本征问题分别揭示哪类带隙；
> 3. 逐行读懂 `stpwe_build_system` → `stpwe_solve_omega` → `stpwe_select_mode` → `stpwe_reconstruct_field` 这条完整管线；
> 4. 自己写出一个扫 $k$ 画 Floquet 能带的脚本，并知道怎么判断"算得够不够准"。

所有代码都在 `core/` 目录，运行前先 `startup_stm`（若忘了，翻 [01-matlab-primer.md](01-matlab-primer.md) 复习）。

---

## 一、从物理到数学：把 Maxwell 方程变成一个本征值问题

### 1.1 我们研究的介质：时空周期

ST-PWE 处理的对象是**介电常数既随空间又随时间周期变化**的介质：

$$
\varepsilon(x+\Lambda,t)=\varepsilon(x,t),\qquad \varepsilon(x,t+T)=\varepsilon(x,t).
$$

$\Lambda$ 是空间周期（晶格常数），$T$ 是时间调制周期。定义两个基本频率：

$$
g=\frac{2\pi}{\Lambda}\quad\text{（空间倒格矢）},\qquad
\Omega=\frac{2\pi}{T}\quad\text{（时间调制角频率）}.
$$

本工具箱按惯例做归一化：$\varepsilon_0=\mu_0=c_0=1$，$\Lambda=1$，于是 $g=2\pi$。归一化波矢与频率（全系列统一约定）：

$$
\bar{k}=\frac{k\Lambda}{2\pi},\qquad \bar{f}=\frac{\omega}{g c_0}.
$$

> 💡 **为什么要归一化？** 物理问题里的"单位"在数值计算里只是常数因子。设 $c_0=\Lambda=1$ 后，一切量都以周期为单位来量，$g=2\pi$、$\Omega=2\pi\bar{\Omega}$，代码干净，图和论文的归一化轴也对得上。你只要记住口诀：**"见 $\omega/(gc_0)$ 就是 $\bar{f}$，见 $k\Lambda/(2\pi)$ 就是 $\bar{k}$"**。

### 1.2 一维横向电磁波：写出 Maxwell 方程

考虑沿 $x$ 传播、单一偏振（比如 $E\parallel \hat{y}$、$H\parallel \hat{z}$）的横向波。一维情形下 $E$、$H$ 都是标量，Maxwell 旋度方程（取归一化单位）化简为两个一阶方程：

$$
\frac{\partial E}{\partial x}=-\mu(x,t)\,\frac{\partial H}{\partial t},
\qquad
\frac{\partial H}{\partial x}=-\varepsilon(x,t)\,\frac{\partial E}{\partial t}.
$$

注意右边的系数是 $\varepsilon(x,t)$、$\mu(x,t)$ **本身**（介质参数随时空变化），不是常数——这正是时空晶体的"性格"所在。代码里 $E$、$H$ 之前的导数算子、以及 $\varepsilon,\mu$ 的卷积，都会从下面这套推导里逐项长出来。

### 1.3 双 Fourier 展开：场和介质的"乐高积木"

既然介质是时空双周期的，我们就把它展开成**双 Fourier 级数**（约定：空间指数取 $+$、时间指数取 $-$）：

$$
\varepsilon(x,t)=\sum_{m,n}\varepsilon_{m,n}\,e^{i n g x - i m \Omega t},
\qquad
\mu(x,t)=\sum_{m,n}\mu_{m,n}\,e^{i n g x - i m \Omega t}.
$$

场的 Floquet-Bloch 形式是"载波 × 双谐波组合"（时间因子统一为 $e^{ikx-i\omega t}$）：

$$
E(x,t)=e^{ikx-i\omega t}\sum_{m,n}E_{m,n}\,e^{i n g x - i m \Omega t},
\qquad
H(x,t)=e^{ikx-i\omega t}\sum_{m,n}H_{m,n}\,e^{i n g x - i m \Omega t}.
$$

于是每个 $(m,n)$ 分量的**实际波矢**是 $k+ng$，**实际频率**是 $\omega+m\Omega$——空间谐波 $n$ 把波矢平移 $n g$，时间谐波 $m$ 把频率平移 $m\Omega$。$k$ 与 $\omega$ 是待定的复参量（后面要让它们成为本征值）。

### 1.4 逐项代入：乘积变卷积（关键一步，别跳）

把场的形式代入 Maxwell。以第一个方程为例，注意到 $E=e^{ikx-i\omega t}u(x,t)$，其中 $u$ 是双谐波组合，所以

$$
\frac{\partial E}{\partial x}=e^{ikx-i\omega t}(ik+\partial_x)u.
$$

同理 $\partial_t H=e^{ikx-i\omega t}(-i\omega+\partial_t)v$。代入 $\partial_x E=-\mu\,\partial_t H$ 并约去公共的 $e^{ikx-i\omega t}$：

$$
(ik+\partial_x)u=-\mu\,(-i\omega+\partial_t)v.
$$

现在对每一项取 $(m,n)$ 谐波分量。三个"动作"分别是：

| 操作                                     | 在 Fourier 空间变成                                               | 说明                           |
| ---------------------------------------- | ----------------------------------------------------------------- | ------------------------------ |
| $\partial_x$                           | 乘以$i n g$                                                     | 每个谐波的对角因子             |
| $\partial_t$                           | 乘以$-i m\Omega$                                                | 每个谐波的对角因子             |
| $\mu(\cdot)$ 或 $\varepsilon(\cdot)$ | **离散卷积** $\sum_{m',n'}\mu_{m-m',n-n'}(\cdot)_{m',n'}$ | 实空间乘积 → Fourier 空间卷积 |

逐分量对齐，得到第 $(m,n)$ 个谐波的**代数方程**：

$$
(k+ng)E_{m,n}=\sum_{m',n'}(\omega+m'\Omega)\,\mu_{m-m',n-n'}\,H_{m',n'},
$$

$$
(k+ng)H_{m,n}=\sum_{m',n'}(\omega+m'\Omega)\,\varepsilon_{m-m',n-n'}\,E_{m',n'}.
$$

> 💡 **这一行是整个方法的核心。** 你看到右侧和号里 $(\omega+m'\Omega)$ 了吗？它把"时间调制在频率上平移了 $m'\Omega$"这件事，变成了一个**只差一个系数 $\omega$** 的矩阵问题。偏微分方程 → 代数本征值问题，靠的就是这一步。

### 1.5 拼成矩阵：$A(k)\Phi=\omega B\Phi$

把 $E_{m,n}$、$H_{m,n}$ 按某种顺序排成两个 $S$ 维列向量（$S=$ 保留的谐波总数），再竖着拼成 $2S$ 维向量 $\Phi=[E;\,H]$。定义矩阵：

| 符号                      | 定义                                                              | 名字                     |
| ------------------------- | ----------------------------------------------------------------- | ------------------------ |
| $K$                     | $kI+G=\mathrm{diag}(k+ng)$                                      | 空间波矢（含$k$ 平移） |
| $G$                     | $\mathrm{diag}(ng)$                                             | 空间导数                 |
| $W$                     | $\mathrm{diag}(m\Omega)$                                        | 时间导数                 |
| $C_\varepsilon,\ C_\mu$ | $C_\varepsilon(\text{row,col})=\varepsilon_{m_r-m_c,\,n_r-n_c}$ | 介电/磁导率卷积矩阵      |
| $B$                     | $\begin{bmatrix}0 & C_\mu \\ C_\varepsilon & 0\end{bmatrix}$    | 右端块矩阵               |

把上面两个代数方程逐行写进矩阵（第二项移到左边、$\omega$ 留在右边），就得到

$$
\underbrace{\begin{bmatrix}
K & -W C_\mu\\
-W C_\varepsilon & K
\end{bmatrix}}_{A(k)}
\begin{bmatrix}E\\H\end{bmatrix}
=
\omega
\underbrace{\begin{bmatrix}
0 & C_\mu\\
C_\varepsilon & 0
\end{bmatrix}}_{B}
\begin{bmatrix}E\\H\end{bmatrix},
\qquad
A(k)\Phi=\omega B\Phi .
$$

推导核对（请你按行验证一遍，这是最值得自己动手的地方）：

- 第一块行 $K E - W C_\mu H=\omega C_\mu H$，即 $K E=(\omega I+W)C_\mu H$。它的第 $(m,n)$ 行分量正是第一行代数方程 $(k+ng)E_{m,n}=(\omega+m\Omega)\sum\mu H$。✓
- 第二块行同理给出第二行方程。✓

好——**这就是 `stpwe_solve_omega` 里实际组装的矩阵**，我们稍后在代码里会逐行对上。

### 1.6 一分钟推导小抄

```
Maxwell (1D, 横向, 归一化)
   ∂E/∂x = -μ ∂H/∂t ,   ∂H/∂x = -ε ∂E/∂t
        │  代入 Floquet-Bloch 形式
        ▼
每谐波代数方程
   (k+ng)E_mn = Σ (ω+m'Ω) μ_{m-m',n-n'} H_m'n'
   (k+ng)H_mn = Σ (ω+m'Ω) ε_{m-m',n-n'} E_m'n'
        │  排成向量 Φ=[E;H]
        ▼
A(k) Φ = ω B Φ        （广义本征值问题，本征值=复频率 ω）
```

---

## 二、两个本征问题，用途不同

时空周期介质里，$\omega$ 和 $k$ 不再是"给定一个求另一个"的一一对应：两者都可能是复数。于是有两条互补的路线，各自揭示一种带隙。

### 2.1 固定波矢 $k$：广义本征值问题（动量带隙）

**`stpwe_solve_omega`** 解 $A(k)\Phi=\omega B\Phi$，输出 $2S$ 个复频率 $\omega$：

- **实部** $\mathrm{Re}\,\omega$ 画出能带结构；
- **虚部** $\mathrm{Im}\,\omega\neq 0$ 意味着这个 $k$ 下没有稳定的传播模——场随时间指数增长/衰减，这就是**动量带隙（momentum gap / k-gap）**，是时间周期介质（PTC）特有的"参数放大"带隙。

> ⚠️ **为什么是"广义"本征值问题？** 因为 $\omega$ 不是像 $A\Phi=\omega\Phi$ 那样直接乘在 $\Phi$ 上，而是乘在一个块对角为 0 的矩阵 $B=[0,C_\mu;C_\varepsilon,0]$ 上——它来自"$\omega C_\mu H$"这个耦合项。$B$ 不是单位阵，所以要用 MATLAB 的广义本征求解 `eig(A, B)`（内部走 QZ 算法）。

### 2.2 固定频率 $\omega$：普通本征值问题（频率带隙）

把代数方程重新整理，让 $k$ 当本征值：

$$
\underbrace{\begin{bmatrix}
-G & (\omega I+W)C_\mu\\
(\omega I+W)C_\varepsilon & -G
\end{bmatrix}}_{A'(\omega)}
\begin{bmatrix}E\\H\end{bmatrix}
=
k\begin{bmatrix}E\\H\end{bmatrix},
\qquad
A'(\omega)\Phi=k\Phi .
$$

因为 $k$ 只出现在对角线上，这是**普通本征值问题**，直接 `eig(A)` 即可。输出的复 $k$：$\mathrm{Im}\,k\neq 0$ 表示该频率在空间上只能以**倏逝波（evanescent）**存在——这就是**频率带隙（frequency gap）**，对应普通空间光子晶体里的禁带，只不过现在是 $k$ 变复。

### 2.3 一句话总结

| 本征问题                   | 固定量     | 本征值               | $\mathrm{Im}(\cdot)\neq 0$ 表示 | 揭示的带隙                         | 求解器                |
| -------------------------- | ---------- | -------------------- | --------------------------------- | ---------------------------------- | --------------------- |
| 广义$A\Phi=\omega B\Phi$ | $k$      | $\omega$（复频率） | 时间指数增长/衰减                 | **动量带隙**（时间不稳定性） | `stpwe_solve_omega` |
| 普通$A'\Phi=k\Phi$       | $\omega$ | $k$（复波矢）      | 空间倏逝                          | **频率带隙**（空间禁带）     | `stpwe_solve_k`     |

> 💡 一句话记忆：**固定 $k$ 看到"时间的带隙"，固定 $\omega$ 看到"空间的带隙"。** 二者不能互换：前者是时间不稳定性，后者是空间倏逝，物理完全不同。`demos/demo02_complex_frequency_and_momentum_gaps.m` 把两条路线画在同一套介质上，值得跑一遍。

---

## 三、代码导读：逐函数走一遍

> 行号以磁盘上 `core/` 的实际文件为准（写本文时如此；若你本地版本不同，以 `help` 和文件为准）。每个函数都是一小段独立讲，建议边读边在编辑器打开对应 `.m`。

### 3.1 `stpwe_build_system.m` —— 组装一切

这是所有求解器的**唯一数据源**：把 Fourier 系数、卷积矩阵、导数矩阵全部打包进结构体 `sys`。签名（第 1 行）：

```matlab
function sys = stpwe_build_system(epsCoeff, muCoeff, Nspace, Mtime, g, Omega)
```

**输入里最关键的是 `epsCoeff`——它必须是一个函数句柄** `@(m,n)`，返回系数 $\varepsilon_{m,n}$。也就是"给一对 $(m,n)$，给一个数"。第 9–11 行给 `muCoeff` 一个默认值 $\mu\equiv 1$：

```matlab
if nargin < 2 || isempty(muCoeff)
    muCoeff = @(m,n) double(m == 0 && n == 0);   % 只有 (m,n)=(0,0) 时返回 1，其余返回 0
end
```

**第 13–16 行：建谐波索引网格。** 这是理解一切排序的关键：

```matlab
[NN, MM] = ndgrid(-Nspace:Nspace, -Mtime:Mtime);
nList = NN(:);     % 空间谐波索引 n 的列向量
mList = MM(:);     % 时间谐波索引 m 的列向量
S = numel(nList);  % 谐波总数 = (2*Nspace+1)*(2*Mtime+1)
```

`ndgrid` 生成 $(2N_{\text{space}}+1)\times(2M_{\text{time}}+1)$ 的网格，再按列主序展平。**记住这个顺序：`mList` 变化得慢、`nList` 变化得快**——后面很多"第几行对应哪个 $(m,n)$"都靠它。

**第 18–27 行：填卷积矩阵（双重 Toeplitz 结构）。**

```matlab
for row = 1:S
    for col = 1:S
        dm = mList(row) - mList(col);   % 行、列的 m 差
        dn = nList(row) - nList(col);   % 行、列的 n 差
        Ceps(row,col) = epsCoeff(dm,dn);
        Cmu(row,col) = muCoeff(dm,dn);
    end
end
```

矩阵元只依赖 $(m_r-m_c,\,n_r-n_c)$——这正是实空间乘积 $\varepsilon(x,t)E(x,t)$ 在 Fourier 空间变成卷积的体现。

**第 29–42 行：打包 `sys` 结构体。** 字段一览：

| 字段                          | 大小            | 含义                                         |
| ----------------------------- | --------------- | -------------------------------------------- |
| `sys.Nspace`, `sys.Mtime` | 标量            | 截断阶数                                     |
| `sys.nList`, `sys.mList`  | $S\times 1$   | 每个谐波的$(n,m)$ 索引（列主序）           |
| `sys.S`                     | 标量            | 谐波总数                                     |
| `sys.g`, `sys.Omega`      | 标量            | 倒格矢与调制频率                             |
| `sys.Ceps`, `sys.Cmu`     | $S\times S$   | 介电/磁导率卷积矩阵                          |
| `sys.G`                     | $S\times S$   | $\mathrm{diag}(n\, g)$，空间导数           |
| `sys.W`                     | $S\times S$   | $\mathrm{diag}(m\,\Omega)$，时间导数       |
| `sys.I`, `sys.Z`          | $S\times S$   | 单位阵、零阵                                 |
| `sys.Bomega`                | $2S\times 2S$ | $[Z, C_\mu;\, C_\varepsilon, Z]$，即 $B$ |

第 38–42 行把它们一个个写进 `sys`，其中：

```matlab
sys.G = diag(nList*g);       % 对应推导里的 G
sys.W = diag(mList*Omega);   % 对应推导里的 W
sys.Bomega = [sys.Z, Cmu; Ceps, sys.Z];   % 对应推导里的 B
```

> 💡 结构体就是"集装箱"：`stpwe_build_system` 只负责**组装一次**，之后扫几百个 $k$ 都复用同一个 `sys`，不用重算卷积矩阵。这就是把"组装"和"求解"拆成两个函数的原因。

### 3.2 `stpwe_solve_omega.m` —— 固定 $k$，求复 $\omega$

**第 8–11 行：组装 $A$ 和 $B$，和推导完全对应。**

```matlab
K = k*sys.I + sys.G;                 % diag(k + n*g)
A = [K, -sys.W*sys.Cmu; ...
    -sys.W*sys.Ceps, K];             % A(k)
B = sys.Bomega;                      % B = [Z, Cmu; Ceps, Z]
```

**第 13–14 行：调用广义本征求解。**

```matlab
[R, D, L] = eig(A, B);   % 解 A*R = B*R*D，D 对角元即本征值
omega = diag(D);
```

`eig(A,B)` 返回三样东西：`R` 的**列**是右本征向量（就是 $\Phi=[E;H]$），`L` 的列是左本征向量，`D` 是对角阵。矩阵规模是 $2S\times 2S$：**前 $S$ 行对应 $E$ 系数，后 $S$ 行对应 $H$ 系数**。

**第 16–30 行：双正交归一化。** 因为系统非厄米，左右本征向量满足的是双正交关系 $L_q'B R_q=\delta_{qq}$。代码逐列算 `overlap = L(:,q)'*B*R(:,q)`，用 $1/\text{overlap}$ 缩放 `R(:,q)`，使 $L_q'BR_q=1$；若自正交（overlap 太小或非有限）则退化为分别按范数归一化。这一步保证后面算占比时分子分母量纲一致。

**第 32–38 行：计算 `m0Weight`（m=0 扇区参与度）。**

```matlab
mask0 = (sys.mList == 0);          % 选出时间直流（m=0）的谐波
E = R(1:S,:);  H = R(S+1:end,:);   % 拆出 E、H 系数块
numerator   = sum(abs(E(mask0,:)).^2 + abs(H(mask0,:)).^2, 1);
denominator = sum(abs(E).^2 + abs(H).^2, 1);
m0Weight    = real(numerator./max(denominator, eps));
```

$$
\text{m0Weight}=\frac{\sum_{n}\bigl(|E_{0,n}|^2+|H_{0,n}|^2\bigr)}{\sum_{n,m}\bigl(|E_{m,n}|^2+|H_{m,n}|^2\bigr)}.
$$

物理含义：这个模的能量有多大比例落在 $m=0$（时间直流）扇区。**越接近 1，越像普通静态介质的 Bloch 模；越小，越是高次 Floquet 边带的混合模。** 它是论文"投影到静态带的权重"的一个透明可复现代理。

**第 40–46 行：输出 `sol` 结构体。**

| 字段                 | 含义                                    |
| -------------------- | --------------------------------------- |
| `sol.k`            | 输入的波矢                              |
| `sol.omega`        | $2S\times 1$，全部本征值（复频率）    |
| `sol.R`, `sol.L` | 左右本征向量矩阵（已双正交归一化）      |
| `sol.A`, `sol.B` | 组装的矩阵（排错时有用）                |
| `sol.m0Weight`     | $2S\times 1$，每个本征值的 m=0 参与度 |

> ⚠️ **本征值排序不稳定。** `eig` 不保证本征值的顺序在相邻 $k$ 之间连续——所以不能"第一个本征值就是第一能带"。这正是 `stpwe_select_mode` 存在的理由。

### 3.3 `stpwe_solve_k.m` —— 固定 $\omega$，求复 $k$

**第 4–6 行：组装 $A'(\omega)$。**

```matlab
OW = omega*sys.I + sys.W;            % diag(omega + m*Omega)
A = [-sys.G, OW*sys.Cmu; ...
     OW*sys.Ceps, -sys.G];
```

**第 8–9 行：普通本征问题。**

```matlab
[R, D, L] = eig(A);   % A*R = R*D，没有第二个矩阵！
k = diag(D);
```

**行/列结构与 `solve_omega` 的对应关系：**

| 位置             | `solve_omega`（固定 $k$）   | `solve_k`（固定 $\omega$）     |
| ---------------- | ------------------------------- | ---------------------------------- |
| 矩阵左上块       | $K=kI+G$                      | $-G$                             |
| 矩阵对角上的变量 | $k$（在 $K$ 里）            | $\omega$（在 `OW` 里）         |
| 本征值           | $\omega$                      | $k$                              |
| 右端矩阵         | $B=[0,C_\mu;C_\varepsilon,0]$ | 无（单位阵）                       |
| 归一化           | $L'BR=I$                      | $L'R=I$                          |
| 返回字段         | 有`B`、`m0Weight`           | **没有** `B`、`m0Weight` |

> ⚠️ 用 `solve_k` 的结果画图时，`sol` 里没有 `nList`/`mList`——需要自己用 `sys.nList`、`sys.mList` 补。这是两个求解器返回结构体不同最容易踩的坑。

### 3.4 `stpwe_select_mode.m` —— 挑一个"干净"的模式

能带图里每个 $k$ 有 $2S$ 个本征值，绝大多数是冗余的 Floquet 副本或数值噪声。`stpwe_select_mode` 帮你在目标点 $(k_{\text{target}},\omega_{\text{target}})$ 附近挑出**物理上最相关的那一个**。

**第 4–7 行：`omegaWindow` 可省略**，默认 `[-inf, inf]`；demo 里传 `[0, 0.9*p.g*p.c0]` 把高频/负频副本滤掉。

**第 8–15 行：求解并过滤**——`stpwe_solve_omega(sys, kTarget)`，只保留有限值且在窗口内的本征值；找不到就 `error`。

**第 17–22 行：打分公式（本函数的核心）。**

```matlab
scale = max(abs(omegaTarget), sys.Omega);
score = abs(real(sol.omega(ids)) - omegaTarget)/scale ...   % 离目标频率近
      + 0.15*abs(imag(sol.omega(ids)))/scale ...            % 虚部尽量小
      + 0.02*(1 - sol.m0Weight(ids));                       % m=0 参与度尽量高
[~, localId] = min(score);
```

三项惩罚分别是：频率偏差、非零虚部（增长/衰减）、m=0 参与度过低。**取分数最小者。** 权重 0.15 和 0.02 是设计好的：在"够近"的前提下，优先要虚部小、再要像静态带。

> 💡 **为什么要挑模式？** 因为本征值排序不稳定，直接画"第 $j$ 个本征值 vs $k$"会得到一堆乱跳的折线。按物理标准（频率最近、虚部最小、m=0 主导）挑出来的才是能带图上那条干净的带。

**第 24–26 行：相位锁定。** 本征向量整体乘以任意 $e^{i\varphi}$ 仍是本征向量——这个自由度会导致不同 $k$ 的场图颜色随机翻转。代码找到 $E$ 系数中模最大的分量，整体旋转使它的相位归零：

```matlab
vec = sol.R(:,id);
[~, phaseId] = max(abs(vec(1:sys.S)));   % 找最强 E 分量
vec = vec*exp(-1i*angle(vec(phaseId)));  % 旋到正实轴
```

**第 28–37 行：输出 `mode` 结构体**——`mode.k`、`mode.omega`、`mode.Ecoef`（$E$ 系数，已锁相）、`mode.Hcoef`、`mode.nList`/`mode.mList`/`mode.g`/`mode.Omega`（从 `sys` 拷来，重构场时用）、`mode.m0Weight`。

### 3.5 `stpwe_reconstruct_field.m` —— 把模式变回场

有了一个模式的系数，就可以把实电场在任意 $(x,t)$ 网格上重建出来。它实现的就是 1.3 节那个公式：

$$
E(x,t)=\operatorname{Re}\left[\,e^{ikx-i\omega t}\sum_{m,n}E_{m,n}\,e^{ingx-im\Omega t}\right].
$$

**第 7–9 行：`includeGrowth` 默认为 `false`。**

**第 11–16 行：先叠周期部分**（对每个谐波求和）：

```matlab
[XX, TT] = meshgrid(x, t);
periodicPart = complex(zeros(size(XX)));
for q = 1:numel(mode.Ecoef)
    periodicPart = periodicPart + mode.Ecoef(q).*exp(1i*( ...
        mode.nList(q)*mode.g*XX - mode.mList(q)*mode.Omega*TT));
end
```

**第 18–24 行：再乘载波。`includeGrowth` 的作用在这里：**

```matlab
if includeGrowth
    omegaForPlot = mode.omega;        % 用完整的 ω（含虚部 → 指数增长/衰减）
else
    omegaForPlot = real(mode.omega);  % 丢掉虚部 → 干净的传播波形
end
carrier = exp(1i*(mode.k*XX - omegaForPlot*TT));
field = real(carrier.*periodicPart);
```

**第 26–29 行：归一化到 $[-1,1]$**——`field = field/max(abs(field(:)))`，让各子图用统一的 `caxis([-1 1])` 对比。

> 💡 输出矩阵形状是 `numel(t) × numel(x)`（行是时间、列是空间），配合 `imagesc(x, t, field)` 后要 `set(gca,'YDir','normal')`，否则纵轴颠倒。

### 3.6 辅助件一览

| 函数                                  | 一句话用途                                                                                 |
| ------------------------------------- | ------------------------------------------------------------------------------------------ |
| `stpwe_static_bands`                | 关掉时间调制（$m=0$ 系数）后算**静态光子能带**，当背景参照线                       |
| `stpwe_fold_frequency`              | 把$\mathrm{Re}\,\omega$ 折进 $[-\Omega/2,\Omega/2)$，虚部保留，用于可视化 Floquet 副本 |
| `stpwe_sample_fourier_coefficients` | 对任意$\varepsilon(x,t)$ 数值采样 Fourier 系数，得到系数表                               |
| `stpwe_lookup_coefficient`          | 从系数表查$(m,n)$ 系数，越界补零                                                         |

查表用的标准搭配（若你的介质写不出解析系数）：

```matlab
table = stpwe_sample_fourier_coefficients(epsFun, -1:1, -10:10, Lambda, T, 512, 256);
epsCoeff = @(m,n) stpwe_lookup_coefficient(table, m, n);
```

详见 [本文章节 3.6](#36-辅助件一览)。

---

## 四、完整可运行例子：扫 $k$，画 Floquet 能带

下面这个脚本仿照 `demos/demo01_reproduce_fig2_stpwe.m` 的 quick 档（`Nspace=10, Mtime=1`，介质取 `stm_fig2_parameters()`），把推导和代码导读串成一条能直接跑的线。**假设你已经运行过 `startup_stm`**（多跑一次也无妨）。

```matlab
%% 04-st-pwe 示例：ST-PWE 扫 k，画 Floquet 能带（仿 demo01 quick 档）
startup_stm;                       % 加入路径（可重复执行）

p = stm_fig2_parameters();         % 归一化参数：Lambda=1, c0=1, g=2*pi, Omega=0.4*pi
Nspace = 10;                       % 空间谐波 n = -10..10，共 21 个
Mtime  = 1;                        % 时间谐波 m = -1..1，共 3 个（正弦调制，1 阶已精确）
epsCoeff = @(m,n) stm_fig2_eps_coeff(m,n,p);   % 解析 Fourier 系数
sys = stpwe_build_system(epsCoeff, [], Nspace, Mtime, p.g, p.Omega);
fprintf('S = %d，广义本征矩阵 %d x %d\n', sys.S, 2*sys.S, 2*sys.S);

Nk = 101;                          % k 采样点数（奇数，让 k=0 落在采样点上）
kBar = linspace(-0.5, 0.5, Nk);    % 归一化波矢，覆盖第一布里渊区
kValues = p.g * kBar;              % 物理波矢 k = kBar * g

fMax = 0.82;                       % 频率上限（归一化 fBar）
allK1 = cell(Nk,1); allF1 = cell(Nk,1); allW1 = cell(Nk,1);   % 严格筛选 → 能带图
allK2 = cell(Nk,1); allF2 = cell(Nk,1); allIm2 = cell(Nk,1);  % 宽松筛选 → Im 图

for ik = 1:Nk
    sol = stpwe_solve_omega(sys, kValues(ik));   % 解 A*Phi = omega*B*Phi
    fBar = sol.omega / (p.g * p.c0);             % 归一化频率

    % 严格筛选：只留"干净"传播模（能带图用）
    keep1 = isfinite(fBar) & real(fBar)>=0 & real(fBar)<=fMax ...
          & abs(imag(fBar))<=3e-3 & sol.m0Weight>=0.035;
    allK1{ik} = repmat(kBar(ik), sum(keep1), 1);
    allF1{ik} = real(fBar(keep1));
    allW1{ik} = sol.m0Weight(keep1);

    % 宽松筛选：保留近带隙的复模，好让 Im(omega) 现形（图 2 用）
    keep2 = isfinite(fBar) & real(fBar)>=0 & real(fBar)<=fMax ...
          & abs(imag(fBar))<=0.03 & sol.m0Weight>=0.01;
    allK2{ik} = repmat(kBar(ik), sum(keep2), 1);
    allF2{ik} = real(fBar(keep2));
    allIm2{ik} = imag(fBar(keep2));
end
allK1 = vertcat(allK1{:}); allF1 = vertcat(allF1{:}); allW1 = vertcat(allW1{:});
allK2 = vertcat(allK2{:}); allF2 = vertcat(allF2{:}); allIm2 = vertcat(allIm2{:});

% 图 1：Re(omega) 能带，颜色 = m=0 参与度（红 = 干净传播模）
figure('Color','w');
scatter(allK1, allF1, 8, allW1, 'filled');
colormap([linspace(1,0,256)', zeros(256,1), ones(256,1)]);   % 品红(低权重) → 蓝(高权重)
caxis([0 1]); colorbar;
xlabel('$k\Lambda/(2\pi)$', 'Interpreter','latex');
ylabel('$\omega/(gc_0)$', 'Interpreter','latex');
title('Floquet band structure (color = m=0 participation)');

% 图 2：Im(omega) vs k —— 动量带隙一目了然
figure('Color','w');
scatter(allK2, allIm2, 8, allF2, 'filled');
yline(0, 'k:');                      % 零线
xlabel('$k\Lambda/(2\pi)$', 'Interpreter','latex');
ylabel('$\mathrm{Im}(\omega)/(gc_0)$', 'Interpreter','latex');
title('Imaginary part: momentum gap = Im(omega) ~= 0');

% 顺手演示 stpwe_select_mode + stpwe_reconstruct_field：
% 在 kBar=0.10 处挑一个 m=0 参与度高的传播模，重构它的电场时空图
mode = stpwe_select_mode(sys, p.g*0.10, p.g*p.c0*0.45, [0, 0.9*p.g*p.c0]);
fprintf('kBar=%.3f  fBar=%.4f%+.2ei  m0Weight=%.3f\n', ...
    mode.k/p.g, real(mode.omega/(p.g*p.c0)), imag(mode.omega/(p.g*p.c0)), mode.m0Weight);

xField = linspace(0, 3*p.Lambda, 451);
tauField = linspace(0, 3, 361);
field = stpwe_reconstruct_field(mode, xField, tauField*p.T, false);
figure('Color','w');
imagesc(xField/p.Lambda, tauField, field);
set(gca, 'YDir','normal');
caxis([-1 1]); colormap(stm_redblue(256));
xlabel('x/\Lambda'); ylabel('t/T');
title('Selected mode field (includeGrowth=false)');
```

### 怎么解释输出

1. **图 1** 里每一条"带"是一个 Floquet 能带，横轴 $-0.5$ 到 $0.5$ 就是整个第一布里渊区；颜色越偏蓝（m=0 参与度接近 1）的模式越像普通传播模，越偏品红越是高次时间扇区混合。
2. **带隙处带子会"断开"或变稀**——因为严格筛选把 $\mathrm{Im}\,\omega$ 大的模滤掉了，没有点落在真正的动量带隙内部。
3. **图 2** 是"现形"的地方：若某个 $k$ 区间出现**关于 $y=0$ 成对对称**的散点（$+\mathrm{Im}\,\omega$ 与 $-\mathrm{Im}\,\omega$），那就是动量带隙，场在这个区间随时间指数增长/衰减（参数放大对，互为共轭）。
4. 对照物理背景：时间调制频率 $\bar{\Omega}=0.2$ 决定了带隙的位置；把 `stm_fig2_parameters` 里 `modDepth` 调大，图 2 的带隙会变宽，调小会变窄甚至消失——这是检验你是否理解"调制深度 → 动量带隙宽度"的好实验。
5. 最后那行打印的 `fBar` 应落在你给定的目标频率（0.45）附近、虚部很小、`m0Weight` 接近 1——说明 `stpwe_select_mode` 挑对了"干净"模。

---

## 五、收敛与截断：怎么知道算得够不够

### 5.1 矩阵大小怎么涨

谐波总数和矩阵规模：

$$
S=(2N_{\text{space}}+1)(2M_{\text{time}}+1),
\qquad \text{本征矩阵 } 2S\times 2S.
$$

| 截断                            | $S$              | 本征矩阵          |
| ------------------------------- | ------------------ | ----------------- |
| `Nspace=10, Mtime=1`（quick） | $21\times 3=63$  | $126\times 126$ |
| `Nspace=20, Mtime=1`（paper） | $41\times 3=123$ | $246\times 246$ |
| `Nspace=20, Mtime=2`          | $41\times 5=205$ | $410\times 410$ |

`eig` 对稠密矩阵的开销约是矩阵阶数的三次方量级，所以**矩阵维数按乘积增长，代价上升非常快**：把 `Nspace` 从 10 加到 20，矩阵从 126 到 246，单次 `eig` 大约贵 $ (246/126)^3\approx 7$ 倍。这不是写错代码，是方法本身的性质——所以"够用就好"。

### 5.2 怎么判断够不够：收敛审计的思想

`demos/demo12_ptc_convergence_audit.m` 提供了一个可照搬的思路：**拿一个独立参照，扫描截断阶数，看误差怎么降**。

- 它研究的是方波时间晶体，参照是精确的 2×2 单周期矩阵（monodromy）；
- 对 `Mtime = [3 5 9 13 19]` 各算一遍能带，按分支匹配后记录误差 $|\Delta\omega|/\Omega$；
- 打印每个截断的 `median`/`max` 误差和耗时——**误差应随截断增大而下降，且你要知道下降到多少才够**。

**几个实用判据，够你的例子用：**

1. **变一个量，结果基本不变**：把 `Nspace` 从 10 加到 20（或 30），如果图 1 的带、图 2 的带隙位置/宽度都不变，说明空间方向已收敛。demo01 的 quick/paper 两档正是干这个。
2. **正弦调制下 `Mtime=1` 是精确的**：因为 $\varepsilon_{m,n}$ 只有 $m=0,\pm1$ 非零（看 `stm_fig2_eps_coeff` 的 `switch m`），时间方向不引入截断误差，要收敛测试的只有 `Nspace`。**方波调制则相反**：系数只按 $1/|m|$ 衰减，需要 `Mtime` 15–30，且收敛慢得多。
3. **对比独立方法**：拿 `stpwe_static_bands` 的静态带 + Floquet 副本（±$\bar{\Omega}$ 平移）当参照，看动态带是否落在预期位置；或像 demo03 那样和精确 TMM 并排画。

### 5.3 关于 `Nspace` 的奇偶选择

常见的一个问题是"介质有对称性，能不能只取半边谐波省一半计算"。**答案：本工具箱没这么做。** 我核对了 `stpwe_build_system` 和 `demo01`：代码总是保留完整区间 $n=-N_{\text{space}}\ldots N_{\text{space}}$，不利用空间对称性做非负/非正半边截断。原因有二：

- $2N_{\text{space}}+1$ 对任意整数 `Nspace` 都是奇数，**正中间永远有一个 $n=0$ 的直流谐波**——它对应介质的空间平均值，不能省；
- 卷积矩阵要求 $C_\varepsilon(\text{row,col})=\varepsilon_{n_r-n_c}$，完整索引集合才能保证 Toeplitz 结构自洽。

所以对这套代码，你不需要纠结奇偶：`Nspace=10` 和 `11` 都行，区别只是多一对谐波、多花一点时间。真正要检查的是收敛（5.2 节），而不是"用奇数还是偶数"。如果你自己写对称介质的半边截断，那是另一个独立实现，`stpwe_*` 不提供这个选项。

---

## 本章速查

| 你想干什么                            | 用什么                                                               |
| ------------------------------------- | -------------------------------------------------------------------- |
| 组装卷积矩阵系统                      | `stpwe_build_system(epsCoeff, [], Nspace, Mtime, g, Omega)`        |
| 固定$k$ 求复 $\omega$（动量带隙） | `stpwe_solve_omega(sys, k)`                                        |
| 固定$\omega$ 求复 $k$（频率带隙） | `stpwe_solve_k(sys, omega)`                                        |
| 在目标点挑一个干净模式                | `stpwe_select_mode(sys, kT, omegaT, omegaWindow)`                  |
| 把模式重构成电场图                    | `stpwe_reconstruct_field(mode, x, t, includeGrowth)`               |
| 静态参照能带                          | `stpwe_static_bands(kValues, epsCoeff0, Nspace, g, c0, nBands)`    |
| 任意材料的数值 Fourier 系数           | `stpwe_sample_fourier_coefficients` + `stpwe_lookup_coefficient` |

**一句话记住整章**：把 $\varepsilon(x,t)$ 和场展开成双 Fourier 级数，代入 1D Maxwell，乘积变卷积、求导变对角阵，于是 $k$ 当参数得到 $A(k)\Phi=\omega B\Phi$；$\omega$ 当参数得到 $A'(\omega)\Phi=k\Phi$。前者看动量带隙，后者看频率带隙。

## 下一步

- 想亲手跑一遍对比复谱的完整 demo → `demos/demo02_complex_frequency_and_momentum_gaps.m`（`startup_stm` 后直接运行）
- 想看懂论文图 2 的完整复现 → `demos/demo01_reproduce_fig2_stpwe.m`（本章例子的"官方完整版"）
- 想理解收敛审计怎么做 → `demos/demo12_ptc_convergence_audit.m`
- 需要每个函数的签名速查 → [docs/tool-reference.md](../docs/tool-reference.md)
- 想复习 MATLAB 结构体、函数句柄、`for` 循环 → [01-matlab-primer.md](01-matlab-primer.md)
