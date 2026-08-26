# 06 · FDTD 时域仿真精讲

> 目标读者：**物理基础 OK、已经会基本 MATLAB**（读过 [01-matlab-primer.md](01-matlab-primer.md)）的读者。
> 这一章讲工具箱的**第三套引擎**：`fdtd/fdtd1d_db.m`。你会看懂它怎么把 Maxwell 方程离散到时空网格上"蛙跳"，每个 `cfg` 字段是干什么的，以及怎么用它仿真一个**时间界面**。

读完这一章，你应该能做到：
1. 画出一维 Yee 网格的时空排布，说清楚为什么 E、H 要错开半个格点、半个时间步；
2. 说清为什么本工具箱用 **D/B**（而不是 E/H）做推进变量；
3. 看懂 `cfg` 全部字段的含义、默认值与使用建议；
4. 知道 `out` 里装了什么、什么时候该用 `storeFields=false` 的低内存模式；
5. 会用 CFL 条件判断"步长合不合法"，遇到 NaN/发散知道先查什么；
6. 独立搭一个"复高斯波包穿过单个时间界面"的仿真，并把数值透/反射系数和 Morgenthaler 解析值对比。

> 💡 本工具箱的物理约定（时间因子 $e^{ikx-i\omega t}$、归一化单位 $\varepsilon_0=\mu_0=c_0=1$、$\Lambda=1$）见 [04 章 §1.1](04-st-pwe.md)（归一化）与 [02 章](02-physics-background.md)（物理背景），读代码前建议先扫一遍。归一化波数 $\bar{k}=k\Lambda/(2\pi)$、归一化频率 $\bar{f}=\omega/(g c_0)$。

---

## 一、从"场方程"到"时空网格"：FDTD 的基本思想

你在麦克斯韦方程里学的 $\nabla\times$、$\partial_t$ 都是**连续**的。FDTD（Finite-Difference Time-Domain，时域有限差分）的核心想法只有一句话：

> **把连续的时空用一张网格盖住，把导数换成差分，然后一步一步往前推。**

别的引擎（ST-PWE、时间 TMM）都是先假设场长成某种"标准形状"（Fourier 级数、平面波），再求本征值。FDTD 不问形状，直接从初值出发、逐时间步**演化**。所以它最擅长：任意时空变化的材料、波包演化、"这道题我不知道解析解"的情况。

### 1.1 一维：波沿 $x$ 跑，场只和 $x,t$ 有关

一维下，取波沿 $x$ 传播、$E$ 沿 $y$ 偏振、$H$ 沿 $z$ 偏振。Maxwell 方程组退化成两个标量方程：

$$
\frac{\partial B}{\partial t} = -\frac{\partial E}{\partial x}, \qquad
\frac{\partial D}{\partial t} = -\frac{\partial H}{\partial x}.
$$

注意我们写的是 $D$ 和 $B$，而不是 $E$ 和 $H$——这是本工具箱最关键的选择，1.4 节专门讲。配套的本构关系是

$$
D = \varepsilon(x,t)\,E, \qquad B = \mu(x,t)\,H,
$$

其中 $\varepsilon(x,t)$ 可以随时间、空间任意变化（这就是"时空介质"进来的方式）。

### 1.2 一维 Yee 网格：E、H 在空间上错开半个格点

Yee 网格（Yee 1966）的核心技巧：**把 E 和 H 放在错开的格点上**，让中心差分两边对称、二阶精确。

```
E 网格（整数格点）:   x₀     x₁     x₂     x₃     ...    x_{Nx-1}
                      |      |      |      |      |
H 网格（半整数格点）:   xH₀    xH₁    xH₂    ...    xH_{Nx-2}
```

- E、D 定义在整数格点 $x_i$ 上，共 `Nx` 个，间距 $\Delta x$；
- H、B 定义在半整数格点 $x_{H,i}=x_i+\Delta x/2$ 上，正好落在相邻两个 E 格点的正中间。

为什么这样排？因为对 $E$ 求空间导数时，

$$
\left.\frac{\partial E}{\partial x}\right|_{x_{H,i}} \approx \frac{E(x_{i+1}) - E(x_i)}{\Delta x},
$$

右端两个点关于 $x_{H,i}$ 左右对称，所以这是**二阶精确**的中心差分。$H$ 的导数在 E 格点上同理。这是整个算法精度的地基。

> ⚠️ H 网格到底有几个点，取决于边界条件：
> - `boundary = 'periodic'`（周期边界）：最后一个 E 格点也参与环绕的 $\partial E/\partial x$，所以 **H 网格与 E 网格等长（`Nx` 个）**；
> - `boundary = 'sponge'`（吸收边界，默认）：两端不参与差分，所以 **H 网格少一个点（`Nx-1` 个）**。
> 这就是为什么后面 `cfg.Hhalf0` 的长度会和边界类型有关。

### 1.3 蛙跳（leapfrog）：时间上再错开半个步长

空间错开还不够，时间上也错开：**B 定义在半整数时刻，E 定义在整数时刻**，像两只脚交替踩地往前跳——所以叫"蛙跳"。

```
时间轴:  t=-dt/2    t=0      t=dt/2    t=dt     t=3dt/2   t=2dt
           |         |         |         |         |         |
已知量:    B⁰        E⁰        B¹        E¹        B²        E²
           |         |         |         |         |         |
本构关系:  H⁰=B⁰/μ   D⁰=ε·E⁰   H¹=B¹/μ   D¹=ε·E¹   H²=B²/μ
```

每次推进只需要当前时刻的场，不需要存历史（除了你想记录的），所以每一步的代价极小。**重要推论**：初始条件要分两个时刻给：

- `cfg.E0` —— 电场在 $t=0$ 时刻，定义在 E 网格上；
- `cfg.Hhalf0` —— 磁场在 $t=-\Delta t/2$ 时刻，定义在 H 网格上。

> 💡 如果初始场是行波（比如一个高斯波包），`Hhalf0` 不能随便给。行波满足 $H = nE$（阻抗关系），又要在时间上差半个步长，所以要写成 `Hhalf0 = n * profile(xH + v*dt/2)`，其中 $v=1/n$ 是波速、`+v*dt/2` 表示场在 $\Delta t/2$ 时间内前进的距离。第 8 节的完整例子就是这么写的。

### 1.4 为什么推进 D、B 而不是 E、H？

这是本工具箱区别于"教学版 FDTD"最核心的一点。回想 02 章的物理：在**时间界面**（介电常数突然变化的瞬间）上，电磁场要满足特定的连续性条件。对标准的"体介质快速开关"模型，连续的量是

$$
D(t^+)=D(t^-), \qquad B(t^+)=B(t^-),
$$

而 E、H 在界面处是**不连续的**——正因为不连续，才会产生时间折射、时间反射（Morgenthaler 1958）。

| 量 | 时间界面处 | 适不适合当"推进变量" |
|----|-----------|---------------------|
| $D$ | 连续（$D(t^+)=D(t^-)$） | 适合 |
| $B$ | 连续（$B(t^+)=B(t^-)$） | 适合 |
| $E=D/\varepsilon$ | 不连续（$\varepsilon$ 跳变） | 不适合 |
| $H=B/\mu$ | 不连续（$\mu$ 跳变） | 不适合 |

如果直接推进 E/H，每到界面时刻都要手工做"匹配"，又慢又容易错。反过来，推进 D/B 时它们**天然连续**，界面两侧平滑衔接，我们只需要在每个时间步用本构关系把 E、H 临时恢复出来：

$$
E = \frac{D}{\varepsilon(x,t)}, \qquad H = \frac{B}{\mu(x,t)}.
$$

材料怎么变都无所谓，只要本构关系是局域的、线性的。

> ⚠️ 精确地说，D/B 连续是**一种**跳变模型（"bulk"模型）。`tmm/temporal_interface_matrix_jump.m`（[源码](../tmm/temporal_interface_matrix_jump.m)）里还提供了 E/B、D/H、E/H 连续以及自定义 `jumpD`/`jumpB` 的模型（见 [05 章 §2.3](05-temporal-tmm.md#23-通用版本-temporal_interface_matrix_jump跳变律是第-5-个参数)）。而 `fdtd1d_db` 按构造**固定使用 D/B 连续模型**——它更新 D 和 B，本构关系永远是 $E=D/\varepsilon$。想比较不同跳变模型，只能走 TMM 那条路（第 05 章会讲）。

---

## 二、离散方程：三个公式撑起整个引擎

把 1.1 的两个偏微分方程直接换成差分，就是主循环的全部内容。记号：上标 $n$ 表示整数时刻 $t_n=n\Delta t$，$n\pm\tfrac12$ 表示半整数时刻。

**步骤 A —— 法拉第定律（由 $E^n$ 更新 $B^{n+1/2}$）：**

$$
B^{n+1/2}(i) = B^{n-1/2}(i) - \frac{\Delta t}{\Delta x}\big[E^{n}(i+1)-E^{n}(i)\big],
$$

然后恢复 $H$：

$$
H^{n+1/2}(i) = \frac{B^{n+1/2}(i)}{\mu\big(x_{H,i},\,t^{n+1/2}\big)}.
$$

**步骤 B —— 安培定律（由 $H^{n+1/2}$ 更新 $D^{n+1}$）：**

$$
D^{n+1}(i) = D^{n}(i) - \frac{\Delta t}{\Delta x}\big[H^{n+1/2}(i)-H^{n+1/2}(i-1)\big],
$$

然后恢复 $E$：

$$
E^{n+1}(i) = \frac{D^{n+1}(i)}{\varepsilon\big(x_i,\,t^{n+1}\big)}.
$$

注意两点：

- $\Delta t/\Delta x$ 是个**无量纲数**，它直接就是 Courant 数的一部分（第 5 节讲）；
- "curl"（旋度）在一维下就是 $\partial/\partial x$。周期边界下末尾用环绕（`circshift`），海绵边界下用 `diff`。

下面是主循环的**示意代码**（只保留周期边界路径，真实文件里还有时间界面处理和记录逻辑）。对照源码 [fdtd/fdtd1d_db.m](../fdtd/fdtd1d_db.m) 第 686–831 行：

```matlab
for step = 1:nSteps
    % ---- 步骤 A：法拉第定律，更新 B（n-1/2 -> n+1/2）----
    curlE = circshift(E,-1) - E;          % E(i+1)-E(i)，最后一点环绕回第一个点
    B = B - (dt/dx)*curlE;                % B^{n+1/2} = B^{n-1/2} - (dt/dx)*curlE
    B = B .* dampH;                       % 吸收层衰减（H 网格）
    tHalf = (step-0.5)*dt;
    H = B ./ muFun(xH, tHalf);            % 本构关系：H = B/mu

    % ---- 步骤 B：安培定律，更新 D（n -> n+1）----
    curlH = H - circshift(H,1);           % H(i)-H(i-1)，第一个点接最后一个点
    D = D - (dt/dx)*curlH;
    D = D + cfg.sourceD(x, step*dt, step); % 可选软源（见 3.2 节）
    D = D .* dampE;                       % 吸收层衰减（E 网格）
    tNow = step*dt;
    E = D ./ epsFun(x, tNow);             % 本构关系：E = D/eps
end
```

> 💡 `circshift(E,-1)` 是把 `E` **向左**循环移位一位：`[E(2) E(3) ... E(end) E(1)]`。所以 `circshift(E,-1)-E` 就是逐点的 $E(i+1)-E(i)$，末点自动环绕到首点。`circshift(H,1)` 向右移位一位，同理。
>
> 海绵边界下的差别（对照源码第 725–760 行）：`curlE` 用 `diff(E)`（长度正好 `Nx-1`，落在 H 网格上）；D 只更新内部点 `D(2:end-1) = D(2:end-1) - (dt/dx)*(H(2:end)-H(1:end-1))`，两个端点 `D(1)`、`D(end)` 不参与差分，靠 sponge 吸收层压下去。

---

## 三、配置结构体 cfg：一张表看懂全部字段

`fdtd1d_db` 只有一个输入，就是配置结构体 `cfg`；只有两个输出 `out`。用结构体打包几十个参数是 MATLAB 科学计算的常见风格（[01 章 §四](01-matlab-primer.md) 讲过）。

> 下面所有字段名、默认值、校验规则，都直接来自源码 [fdtd/fdtd1d_db.m](../fdtd/fdtd1d_db.m) 第 117–161 行的帮助注释和第 186–389 行的实际代码，保证和代码一一对应。完整字段表以 `help fdtd1d_db` 和本系列 [06 章 §三](06-fdtd.md#三配置结构体-cfg一张表看懂全部字段) 为权威参考，这里按**教学顺序**重排并给出使用建议。

### 3.1 四个必填字段

| 字段 | 类型 | 含义 | 使用建议 |
|------|------|------|---------|
| `cfg.x` | 数值行/列向量 | E/D 场的均匀空间网格。至少 5 个点，且必须是均匀的（代码用 `diff` 检查，不均匀直接报错） | 建议**每波长 ≥ 32 个格点**。先小网格跑通，再加密看收敛 |
| `cfg.dt` | 正实数标量 | 时间步长。必须满足 CFL 条件（第 5 节） | 均匀介质建议 `dt = 0.8*dx/c0`，更保守也可以 |
| `cfg.nSteps` | 非负整数标量 | 总时间步数。总仿真时间 = `nSteps*dt` | 想要仿真时长 `tEnd`，就取 `nSteps = ceil(tEnd/dt)` |
| `cfg.epsFun` | 函数句柄 `@(x,t)` | 相对介电常数函数。输入：空间坐标数组 `x` 和标量时间 `t`；输出：与 `x` 等长的 $\varepsilon_r$ 数组 | 材料"长什么样"完全由你定义。必须返回**和输入等长**的数组（见 3.3 的坑） |

### 3.2 可选字段全表

| 字段 | 默认值 | 类型 | 含义与使用建议 |
|------|--------|------|---------------|
| `cfg.muFun` | `@(xq,tq) ones(size(xq))` | 函数句柄 | 相对磁导率函数，接口与 `epsFun` 一样。大多数情况不用改，保持 $\mu=1$ 即可 |
| `cfg.E0` | 全零 | 数值向量 | 初始电场，$t=0$ 时刻、E 网格上。长度必须等于 `numel(cfg.x)`。不想从零场开始就给个波包（见第 8 节） |
| `cfg.Hhalf0` | 全零 | 数值向量 | 初始磁场，$t=-\Delta t/2$ 时刻、H 网格上。长度注意：**periodic 边界是 `Nx`，sponge 边界是 `Nx-1`**。行波初始条件记得做 $\Delta t/2$ 时移 |
| `cfg.boundary` | `'sponge'` | 字符串 | `'sponge'`（两端海绵吸收）或 `'periodic'`（空间环绕）。时空晶体研究常用 `'periodic'` 配合频谱滤波；波包演化用 `'sponge'` |
| `cfg.spongeCells` | `min(80, floor(Nx/8))` | 正整数 | 吸收层宽度（格点数）。太窄会反射回来污染结果；波包离边界远、跑得长就加大 |
| `cfg.spongeStrength` | `0.12` | 正实数 | 吸收强度。衰减因子是 `exp(-strength*s^3)`，`demo05` 用 `0.08` 偏轻，想更强可以到 `0.2` |
| `cfg.recordEvery` | `1` | 正整数 | 场记录间隔：每 N 步记录一次。总记录数 `nRecords = floor(nSteps/recordEvery)+1`（含 $t=0$ 的初始帧）。嫌 `out.E` 太大就调大它 |
| `cfg.storeFields` | `true` | 逻辑标量 | 是否存完整 E/H 历史矩阵。**`false` 进入低内存模式**：不分配全域历史，只记探针和最终快照（第 9 节） |
| `cfg.storeD` | `false` | 逻辑标量 | 是否独立存 E 网格上的 D 历史。与 `storeFields` 相互独立：即使 `storeFields=false`，也可以 `storeD=true` 只存 D |
| `cfg.probeIndices` | 空 | 整数向量 | 探针点索引（1..Nx 的整数）。**无论 `storeFields` 是什么，探针总被记录**，是低内存模式下看场的关键通道 |
| `cfg.precision` | `'double'` | 字符串 | `'double'` 或 `'single'`。`'single'` 约省 50% 内存、可能更快，精度降到 ~1e-3。大规模/长时仿真才需要 |
| `cfg.progressBar` | `false` | 逻辑标量 | 是否显示文本进度条（每约 2% 刷新一次，无图形界面也能用）。长仿真开着，心里有底 |
| `cfg.temporalInterfaces` | 空 | 严格递增实向量 | **介电常数突变时刻**（时间界面）。见第 6 节。注意：是普通数值**向量**，不是结构体数组 |
| `cfg.temporalInterfaceTolerance` | 自适应（见第 6 节） | 正实数 | 时间界面与整数时刻 $n\Delta t$ 的对齐容差。自己给时必须 $>0$ 且 $<\Delta t/4$ |
| `cfg.spectralFilterMask` | 空 | 数值/逻辑向量（长 `Nx`） | 周期网格 FFT 掩码，取值 $[0,1]$。**要求 `boundary='periodic'` 且至少给一个时间界面**。见 6.4 节 |
| `cfg.stabilityTimes` | `linspace(0, nSteps*dt, 9)` | 数值向量 | CFL 审计的采样时刻。第 5 节讲 |
| `cfg.maxWaveSpeed` | 无（自动采样） | 正实数 | 手动指定最大波速，**跳过采样审计**。色散/负折射率/有源介质必须用它，否则采样审计会跳过或误判 |
| `cfg.sourceD` | 无 | 函数句柄 `@(x,t,n)` | D 场**软源**：每个时间步更新完 D、加 sponge 衰减之前调用，`D = D + sourceD(x, tNow, step)`。软源不会阻挡反射波（硬源会） |

### 3.3 一个最常见的小坑：`epsFun` 必须返回与 `x` 等长的数组

新手写 `epsFun` 最常犯的错误：写成返回**标量**的函数。

```matlab
% 错误写法：tq 是标量，(tq<tSwitch) 是标量，整个表达式返回标量
cfg.epsFun = @(xq,tq) (tq < tSwitch)*nBefore^2 + (tq >= tSwitch)*nAfter^2;

% 正确写法：乘上 ones(size(xq)) 把它铺成和 xq 等长的数组
cfg.epsFun = @(xq,tq) ((tq < tSwitch)*nBefore^2 + (tq >= tSwitch)*nAfter^2) * ones(size(xq));
```

为什么？求解器会调用 `epsFun(x, 0)`，其中 `x` 是整条空间网格（上千个数）。`epsFun` 必须返回**同长度**的 $\varepsilon_r$ 数组，否则代码报 `epsFun or muFun returned an array with the wrong grid size.`。正确写法里的 `*ones(size(xq))` 就是把标量铺开（`tq` 永远是标量，所以 `tq < tSwitch` 是标量逻辑判断，整体是标量乘一维数组）。

> 凡是 `@(xq,tq)` 这种函数，写作习惯都是**末尾补 `*ones(size(xq))`**。你在 `tests/test_smoke.m` 和 `demos/demo06` 里都能看到这个模式。

---

## 四、输出结构体 out：仿真结束你能拿到什么

`fdtd1d_db` 返回一个结构体 `out`，字段如下（全部核对了源码第 843–868 行的组装代码）：

| 字段 | 含义 |
|------|------|
| `out.x` / `out.xH` | E 空间网格 / H 空间网格（相对 E 偏移 `dx/2`） |
| `out.t` | 记录时刻数组（`nRecords×1`，第 1 个是 $t=0$） |
| `out.E` | 电场历史矩阵（`nRecords × Nx`，每行一个时刻） |
| `out.H` | 磁场历史矩阵（`nRecords × Nx`，**已插值到 E 网格**，方便和 `E` 同格点画图） |
| `out.D` | 电位移历史（`nRecords × Nx`）；`storeD=false` 时是 `0×Nx` 空矩阵 |
| `out.probeIndices` / `out.probeX` | 探针索引 / 探针空间位置 |
| `out.probeE` / `out.probeH` | 探针点处的 E / H 时序（`nRecords × nProbe`） |
| `out.energy` | 瞬时电磁能量序列（`nRecords×1`）。公式见下 |
| `out.finalD` / `out.finalB` / `out.finalE` / `out.finalH` | **最终时刻**的全域场快照。注意 B 没有历史矩阵，只有这个最终快照 |
| `out.dx` / `out.dt` | 空间 / 时间步长 |
| `out.boundary` | 实际使用的边界类型 |
| `out.storeFields` / `out.storeD` | 存储模式回显 |
| `out.temporalInterfaces` / `out.temporalInterfaceTolerance` | 吸附到整数节点后的界面时刻 / 容差 |
| `out.spectralFilterMask` | 频谱滤波掩码（没给就是空） |
| `out.sampledMaxWaveSpeed` / `out.sampledCourant` | CFL 审计的最大波速 / Courant 数（非正折射率介质是 `NaN`） |
| `out.precision` | 数值精度字符串 |

能量序列的定义（源码第 654、818 行）是连续形式 $\int \tfrac12(E\cdot D + H\cdot B)\,dx$ 的离散近似：

$$
E_{\mathrm{em}}(t) = \frac{\Delta x}{2}\sum_i \mathrm{Re}\big[E^*(i)D(i) + H^*(i)B(i)\big],
$$

其中 H、B 用的是插值到 E 网格后的值，和 E 在同一些格点上求和。

> ⚠️ **时变介质里总能量一般不守恒**。时间调制会与场交换能量（这正是时间界面能"做功"、PTC 能放大的原因）。所以看到 `out.energy` 在界面处跳变、或者总体上升，先别急着怀疑数值有 bug——对静态介质能量才应当守恒（`test_smoke` 断言漂移 < 8%）。

---

## 五、CFL 稳定性：步长定生死

FDTD 有个著名的"生死线"——**Courant 条件**。一维 Yee 格式要求：

$$
S = \frac{v_{\max}\,\Delta t}{\Delta x} < 1, \qquad
v_{\max} = \max_{x,t}\left[\frac{1}{\sqrt{\varepsilon_r(x,t)\,\mu_r(x,t)}}\right].
$$

物理直觉：信息以波速 $v_{\max}$ 传播，一个时间步内波最多走 $v_{\max}\Delta t$，不能超过一个网格间距 $\Delta x$，否则信息会"追不上"。对真空（$\varepsilon_r=\mu_r=1$）就是 $\Delta t < \Delta x/c_0$。

**不满足会怎样？** 数值解指数增长，很快变成 `NaN` 或 `Inf`，或者场直接"爆表"。

### 5.1 代码的"采样 CFL 审计"

`fdtd1d_db` 在正式推进之前会先做一次**稳定性审计**（源码第 477–533 行）：

1. 如果给了 `cfg.maxWaveSpeed`，直接用这个值，跳过采样；
2. 否则在 `[0, nSteps*dt]` 里均匀取 9 个时刻（可用 `cfg.stabilityTimes` 自定义），每个时刻在整个空间网格上算 $1/\sqrt{\varepsilon_r\mu_r}$，取全局最大值作为 `sampledMaxWaveSpeed`；
3. 若 $\varepsilon_r,\mu_r$ 在某时刻**虚部不可忽略**（$>10^{-12}$）或**实部非正**（色散/损耗/增益/负折射率），审计返回 `NaN` 并**跳过检查**——标量 CFL 对这种介质本来就不适用；
4. 算 `sampledCourant = sampledMaxWaveSpeed*dt/dx`，若 $\ge 1$ 直接报错：

```
Error: Sampled one-dimensional CFL number is 1.023 >= 1.
       Reduce cfg.dt or provide a validated constitutive update.
```

审计结果会原样放进输出：`out.sampledCourant`、`out.sampledMaxWaveSpeed`。跑完看一眼它们，是诊断"为什么炸了"的第一手证据。

### 5.2 怎么调参

- **第一步永远是：把 `cfg.dt` 减半**（其他不动）。如果仿真变稳定了，说明是步长太大，按需逐步加回来；
- 或增大 `dx`（网格变粗），但那样空间精度下降，通常优先减 `dt`；
- 推荐**先小算例探路**：`Nx` 几百点、`nSteps` 几十步，跑通确认数量级对，再放大；
- 对色散、增益、负折射率介质，标量 CFL 不够，要自己用 `cfg.maxWaveSpeed` 指定一个（更严格）的上限。`fdtd1d_db` 本身不处理色散/有源介质——那些需要 ADE 或卷积状态。

> 💡 注意区分两个数：`dt/dx` 是"网格比"，而审计里真正判生死的是 `sampledCourant = v_max*dt/dx`。`demo05` 取 `dt = 0.80*dx`（真空下 `v_max=1` 时 `sampledCourant=0.8`），但介质是 $n=1.5$，`v_max = 1/n ≈ 0.67`，所以 `sampledCourant ≈ 0.53`，余量更大。`demo06` 取 `dt = 0.55*dx/c0`，更保守（那里面 $v_{\max}$ 还要乘上介质速度）。

---

## 六、时间界面：时变介质的关键时刻

FDTD 最亮眼的用法就是仿真**时间界面**——某个时刻整个空间的介电常数突然改变。这时 D/B 连续（1.4 节），E 在界面处不连续，产生时间折射 + 时间反射。

### 6.1 `cfg.temporalInterfaces` 是什么？

**一个严格递增的实数向量**（不是结构体数组，没有子字段）：

```matlab
cfg.temporalInterfaces = [tSwitch];      % 单个界面
cfg.temporalInterfaces = [3 7 11];       % 三个界面
```

每个时刻必须满足：

1. **与 E 的整数时间节点 $n\Delta t$ 对齐**，对齐误差 ≤ `cfg.temporalInterfaceTolerance`；
2. 落在 $[0, \texttt{nSteps}\cdot\texttt{dt})$ 内——**末时刻不允许**，因为界面后还要有半个 B 步来"接住"它；
3. 不同界面必须映射到不同节点（代码里 `round` 后不能重复）。

代码会把输入"吸附"到最近整数节点：`temporalInterfaceNodes = round(t/dt)`，实际用 `nodes*dt` 推进。如果某个时刻离整数节点太远，报错会提示最大对齐误差：

```
Error: Each cfg.temporalInterfaces value must align with an E-grid time
       n*cfg.dt within cfg.temporalInterfaceTolerance. ...
```

> ⚠️ 所以给 `tSwitch` 时最好保证它本来就是 `dt` 的整数倍，比如 `tSwitch = 12`、`dt = 0.025` → `12/0.025 = 480`，完美对齐。默认容差是 `min(max(128*eps(max(1,abs(nSteps*dt))), 1e-12*dt), dt/100)`（`max(1,abs(...))` 是防溢出保护，代码里先把 `timeScale = max(1,abs(nSteps*dt))` 再代入），非常小，就是为了防止"想给整数步却给了个浮点误差"。

### 6.2 事件感知：界面时刻怎么做特殊处理？

代码用一张查表数组 `temporalInterfaceIdByStep` 记住"第几步是界面"。关键细节（源码第 693–719 行）：**循环第 `step` 次更新 B 时，B 更新式的中心时刻是 `(step-1)*dt`**。所以节点 $n$ 的界面，是在第 `step = n+1` 次 B 更新时被"看到"的。

一旦命中界面，求解器就不再直接用当前的 `E` 去算旋度，而是用**两侧 ε 的平均恢复 E**：

```matlab
epsBefore = cfg.epsFun(x, tInterface - dt/4);
epsAfter  = cfg.epsFun(x, tInterface + dt/4);
EForBUpdate = 0.5*(D./epsBefore + D./epsAfter);
```

也就是头注释里写的

$$
E_{\text{interface}} = \tfrac12\left[\frac{D}{\varepsilon(t^-)} + \frac{D}{\varepsilon(t^+)}\right],
$$

其中两侧 $\varepsilon$ 在界面前后各 $\pm\Delta t/4$ 处采样（`temporalInterfaceSideOffset = dt/4`，避开突变点本身、又落在半个步长内）。两侧 $\varepsilon$ 必须有限且非零，否则报错。

### 6.3 为什么是"平均"？

如果不做任何特殊处理，界面那一步会只用界面一侧的 E 去算旋度，等于把一个连续过程当成"前一半用旧材料、后一半用新材料"，引入**一阶时间误差**。用两侧 E 的平均（本质是梯形平均），界面处理在时间上仍保持二阶精度。这是 `test_smoke.m` 里"事件感知 FDTD 必须和精确 D/B monodromy 对得上（误差 < 1e-3）"这条断言能通过的原因。

> 对应地，`epsFun` 在界面两侧必须**定义良好**：`tq < tSwitch` 给一个值、`tq >= tSwitch` 给另一个值（3.3 节那种写法）。如果在 $\pm\Delta t/4$ 处采样到 `NaN` 或 0，代码会报错。

### 6.4 频谱滤波掩码：压制高阶动量噪声

`cfg.spectralFilterMask` 是一个高级选项，但值得知道它存在（源码第 359–389、784–789 行）：

- **前置条件**：`boundary='periodic'`，且至少给一个 `temporalInterfaces`；
- **作用**：在每个时间界面的更新结束后，把 D 和 B **同时**投影回你指定的物理波矢支撑：

```matlab
D = ifft(fft(D).*spectralFilterMask);
B = ifft(fft(B).*spectralFilterMask);
E = D./epsNow;  H = B./muHalf;
```

- **物理动机**：有源时变介质会把舍入噪声里落在高阶动量带隙的非物理 $k$ 分量指数放大。对已知是窄带解析信号的周期问题，在每个真实时间界面处把场"拉回"物理支撑，能避免这种伪增长。D 和 B 必须**一起**投影，否则会破坏 Maxwell 状态的一致性。

新手阶段可以无视它；等你做"窄带脉冲 + 周期边界"的时空晶体时，它会派上用场（`test_smoke` 里有它的一条测试）。

---

## 七、吸收边界 sponge：把反射压到看不见

波包跑到计算域边缘会反射回来，污染结果。`'sponge'` 边界就是在两端各留 `spongeCells` 个格点当"海绵"，每步把边界附近的场乘一个小衰减因子。

衰减剖面的公式（源码第 542–583 行）：

$$
\text{damp}(s) = \exp\big(-\texttt{spongeStrength}\cdot s^3\big), \qquad
s = \frac{\texttt{spongeCells} - \text{dist}}{\texttt{spongeCells}},
$$

其中 `dist` 是该格点到**最近端边界**的距离（以 $\Delta x$ 为单位）。所以：

- $s=0$ 在吸收层**内边界**，`damp=1`，无衰减；
- $s=1$ 在计算域**最外边界**，`damp=\exp(-\texttt{spongeStrength})$，衰减最强。

三次方剖面保证衰减从零平滑过渡，让反射最小。每个时间步：

```matlab
B = B .* dampH;   % H 网格上的衰减
D = D .* dampE;   % E 网格上的衰减
```

> ⚠️ 这个 sponge **不是 CPML**（卷积完美匹配层），只是简化的乘法吸收。吸收强度有限，反射不可能完全消除。两条经验：
> - `spongeCells` 取 `min(80, floor(Nx/8))` 起步；波包离边界近、仿真时间长，就加大到 160–240；
> - `spongeStrength` 默认 `0.12`，`demo05`/`demo06` 用 `0.08`。太小吸不掉、太大本身可能引入反射，一般 0.08–0.3 之间试。
>
> 还有一个通用检查：把 sponge 区域以外的场画出来，看边界反射是否明显（时空图四角有没有"又冒出来"的波）。

---

## 八、完整可运行例子：单个时间界面的波分裂

下面这个脚本几乎就是 `demos/demo05_fdtd_temporal_interface.m`（[源码](../demos/demo05_fdtd_temporal_interface.m)）的核心，逐行加了注释。它让一个复高斯波包在 $n=1.5$ 的介质中前进，在 `tSwitch=12` 时刻介电常数跳到 $n=2.5$，然后对比 FDTD 数值的透射/反射系数与 Morgenthaler 解析值。

先在命令窗口跑 `startup_stm`（把 `fdtd/`、`tmm/` 等目录加进路径），然后运行：

```matlab
%% 单个时间界面的波分裂：FDTD vs Morgenthaler 解析解
% ---------- 1) 物理参数（归一化单位：c0=1, lambda0=1） ----------
lambda0 = 1;
k0 = 2*pi/lambda0;        % 载波波数
nBefore = 1.5;            % 界面之前的折射率
nAfter  = 2.5;            % 界面之后的折射率
tSwitch = 12;             % 时间界面发生的时刻
tEnd    = 30;             % 总仿真时长

% ---------- 2) 网格与时间步 ----------
dx = lambda0/32;          % 每波长 32 个格点
x  = 0:dx:60;             % 空间网格（E 网格）
dt = 0.80*dx;             % 网格比 0.8（真空 CFL=0.8；本介质采样 Courant≈0.53，见 5.2 节）
nSteps = ceil(tEnd/dt);   % 总步数（总时长 = nSteps*dt）

% ---------- 3) 初始行波：复高斯波包，沿 +x 传播 ----------
x0 = 17;  sigma = 3.5;
profile = @(xq) exp(-((xq-x0)/sigma).^2) .* exp(1i*k0*(xq-x0));
E0 = profile(x);                        % t=0 时刻的 E
xH = x(1:end-1) + dx/2;                 % H 网格（sponge 边界：Nx-1 个点）
vBefore = 1/nBefore;                    % 界面之前的波速 v = 1/n
Hhalf0 = nBefore*profile(xH + vBefore*dt/2);   % t=-dt/2 的 H，做 dt/2 时移

% ---------- 4) 组装 cfg ----------
cfg.x = x;
cfg.dt = dt;
cfg.nSteps = nSteps;
cfg.epsFun = @(xq,tq) ((tq < tSwitch)*nBefore^2 ...
                     + (tq >= tSwitch)*nAfter^2) * ones(size(xq));
cfg.muFun = @(xq,tq) ones(size(xq));
cfg.E0 = E0;
cfg.Hhalf0 = Hhalf0;
cfg.boundary = 'sponge';
cfg.spongeCells = 160;     % 吸收层宽度（给足，波包跑得远）
cfg.spongeStrength = 0.08;
cfg.recordEvery = 2;       % 每 2 步记一次，省内存

% ---------- 5) 运行 ----------
out = fdtd1d_db(cfg);
fprintf('sampledCourant = %.4f (应 < 1)\n', out.sampledCourant);

% ---------- 6) 界面后把场分解为前向/后向 ----------
[~, beforeId] = min(abs(out.t - (tSwitch - 2*out.dt)));  % 界面前 2 步的一帧
[~, afterId ] = min(abs(out.t - (tSwitch + 2*out.dt)));  % 界面后 2 步的一帧
Ebefore = out.E(beforeId, :);
Eafter  = out.E(afterId , :);
Hafter  = out.H(afterId , :);      % H 已插值到 E 网格，可直接同格点用
Eplus  = 0.5*(Eafter + Hafter/nAfter);    % 前向（时间折射）波
Eminus = 0.5*(Eafter - Hafter/nAfter);    % 后向（时间反射）波
tauNumeric = norm(Eplus )/norm(Ebefore);
rhoNumeric = norm(Eminus)/norm(Ebefore);

% ---------- 7) 解析 Morgenthaler 系数（D/B 连续模型） ----------
[~, tauExact, rhoExact] = temporal_interface_matrix(nBefore^2, 1, nAfter^2, 1);
fprintf('|tau| exact %.5f, FDTD %.5f;  |rho| exact %.5f, FDTD %.5f\n', ...
    abs(tauExact), tauNumeric, abs(rhoExact), rhoNumeric);

% ---------- 8) 画 |E(x,t)| 时空图 ----------
figure('Color','w');
imagesc(out.x, out.t, abs(out.E));
set(gca,'YDir','normal');
xlabel('x');  ylabel('Time');  title('|E(x,t)|: temporal refraction and reflection');
colorbar;
```

跑完你应该看到类似这样的输出：

```
sampledCourant = 0.5333 (应 < 1)
|tau| exact 0.48000, FDTD 0.48xxx;  |rho| exact 0.12000, FDTD 0.11xxx
```

`sampledCourant ≈ 0.533` 因为介质里 $v_{\max}=1/n_{before}=1/1.5$，乘以 `dt/dx=0.8` 得 $0.53$（不是 0.8，0.8 只是真空下的上限，见 5.2 节的说明）。

数值和解析应吻合到 ~1e-2 量级。图中 $t=12$ 处波包**一分为二**：前向波继续前进（$n$ 变大，波速变慢、波长变短、幅度改变），后向波折返——这就是时间反射。

> 💡 想看懂这几行：`Eplus/Eminus` 是把界面后的场按 $\pm x$ 行波方向分解（在 $n=n_{after}$ 的介质里，$H = n_{after}E_{+} - n_{after}E_{-}$，反解出 $E_{\pm}=\tfrac12(E\pm H/n_{after})$）。解析值 `temporal_interface_matrix`（[源码](../tmm/temporal_interface_matrix.m)）返回的是 **E 场系数**：D/B 连续模型下（本例两侧 $\mu=1$），`tau = 0.5*(epsBefore/epsAfter + sqrt(epsBefore/epsAfter))`、`rho = 0.5*(epsBefore/epsAfter - sqrt(epsBefore/epsAfter))`。代入 $1.5^2\to2.5^2$ 得 $|\tau|=0.48$、$|\rho|=0.12$。

**想升级到时空晶体？** `demos/demo06_fdtd_spacetime_wavepacket.m`（[源码](../demos/demo06_fdtd_spacetime_wavepacket.m)）把 `cfg.epsFun` 换成 `stm_fig2_epsilon(x,t,p)`（`p = stm_fig2_parameters()`，见 [core/stm_fig2_parameters.m](../core/stm_fig2_parameters.m)），再放一个探针 `cfg.probeIndices`，用 `stm_stft`（[源码](../core/stm_stft.m)）做短时傅里叶变换，就能看到 Floquet 边带（频率转换）。逻辑和本例子完全一样，只是介质从"突跳一次"变成"一直振荡"。

---

## 九、调试速查：遇到问题怎么办

FDTD 的报错和怪现象，绝大多数能归到下面几类：

| 症状 | 最可能的原因 | 排查/修复 |
|------|-------------|----------|
| 直接报 `Sampled one-dimensional CFL number is ... >= 1` | `dt` 太大，或介质里有高波速区 | `cfg.dt` 减半；或给 `cfg.maxWaveSpeed`；或查 `epsFun` 是否在某些时刻返回了极小的 $\varepsilon$（$v\propto 1/\sqrt{\varepsilon}$ 会变大） |
| 报 `epsFun or muFun returned an array with the wrong grid size` | `epsFun` 返回了标量/错误长度 | 按 3.3 节补 `*ones(size(xq))`；检查输出长度必须等于输入 `x` |
| 报 `Each cfg.temporalInterfaces value must align ...` | `tSwitch` 不是 `dt` 的整数倍 | 把 `tSwitch` 取成 `dt` 的整数倍，或加 `cfg.temporalInterfaceTolerance` |
| 场出现 `NaN`/`Inf` | 常见于 $E=D/\varepsilon$ 里 $\varepsilon=0$，或 CFL 超限 | 打印几个时刻的 `epsFun(x,t)`，确认没有 0、`NaN`；`dt` 减半；`spongeStrength` 别过大 |
| 静态介质里 `out.energy` 却不守恒 | 网格太粗、边界反射、或 `epsFun` 带了时间依赖 | 静态介质能量应几乎恒定（`test_smoke` 断言漂移 < 8%）。加密网格 + 加大 `spongeCells` 再试 |
| 时变介质里 `out.energy` 上升/跳变 | **可能是真实物理**：时间调制和场交换能量（PTC 参数放大） | 先和 TMM 交叉验证再下结论；若增长快得离谱，回头查 CFL 和 `epsFun` |
| `out.E`、`out.H` 是空矩阵（`0×Nx`） | `storeFields=false` 的低内存模式 | 想拿全域场就设 `storeFields=true`，或用 `probeIndices` + `out.finalE/finalD` 拿探针和最终快照 |
| 波包"原地几乎不动" | 初始条件没做好 | 检查 `Hhalf0` 的 $\Delta t/2$ 时移和阻抗匹配 `H = nE`；`test_fdtd_optimizations.m` 有这类诊断 |
| 时空图四角冒出反射波 | sponge 不够 | 加大 `cfg.spongeCells`；减小 `cfg.spongeStrength` 试试；确认波包仿真结束前没撞到吸收层 |

几条通用心法：

- **看审计值**：跑完先 `out.sampledCourant`、`out.sampledMaxWaveSpeed` 扫一眼，大部分稳定性问题一眼就暴露；
- **低内存跑法**：长时大网格用 `cfg.storeFields = false` + `cfg.probeIndices` + `cfg.recordEvery` 调大 +（可选）`cfg.precision = 'single'`。这时 `out.E/out.H` 是空的，但探针时序、能量、最终快照都在；
- **探针定位**：取某个空间位置的时序，用 `[~, id] = min(abs(out.x - xTarget)); probe = out.E(:, id);`（demo06 第 42–43 行就是这么干的）；
- **要动画**：`animate_fdtd1d(out)`（[源码](../fdtd/animate_fdtd1d.m)）能直接播放 `out` 的波包演化，还能存 `.gif`/`.mp4`；
- **先小后大**：永远先用小网格跑通、确认数量级，再放大参数。一次上大网格，发散了你都不知道是物理还是数值。

---

## 下一步

- 想把这套 D/B 连续模型的解析版搞透（时间界面矩阵、时间多层、能带）→ 读第 05 章时间 TMM（时间传输矩阵）
- 想搭自己的第一个时空晶体仿真 → 读 `demos/demo06` 和 [06 章 §八](06-fdtd.md#八完整可运行例子单个时间界面的波分裂)
- 想知道"三种引擎什么时候用哪个" → 读 [README.md](README.md) 里的学习路线图
- 函数太多记不住 → 随时查 [docs/tool-reference.md](../docs/tool-reference.md) 函数字典
- 动手练习：把第 8 节的例子改成**两个**时间界面（`cfg.temporalInterfaces = [12 20]`），观察能量曲线的两次跳变；再改成 `boundary='periodic'` 试试周期网格
