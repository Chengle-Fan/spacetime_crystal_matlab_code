# 01 · MATLAB 与科学计算入门

> 目标读者：**物理基础 OK、但没写过几行 MATLAB 代码**的读者。
> 这一章不追求教会你全部 MATLAB——只讲**读懂并运行本工具箱**真正需要的东西，并且把"科学计算"的思维方式铺垫好。

读完这一章，你应该能做到：
1. 熟悉 MATLAB 界面、把本项目加入路径、跑通自检测试；
2. 读懂本工具箱代码里最常见的语法：向量/矩阵、结构体、函数句柄、`for` 循环；
3. 能自己画图、改参数、看数值结果；
4. 理解"离散化 → 算 → 检查收敛"的科学计算流程。

---

## 一、MATLAB 界面与工作流

打开 MATLAB，你会看到几个关键窗口：

| 窗口 | 干什么 | 新手最常用操作 |
|------|--------|---------------|
| **命令窗口（Command Window）** | 直接敲代码、看输出 | 输入表达式按回车立即出结果 |
| **编辑器（Editor）** | 写脚本/函数 | 写代码、设断点、按 F5 运行 |
| **工作区（Workspace）** | 显示当前所有变量 | 双击变量查看数值 |
| **当前文件夹（Current Folder）** | 浏览文件 | 双击 `.m` 文件在编辑器打开 |

**一个最重要的概念：工作目录 + 路径（path）。**

- MATLAB 运行时，代码里写 `stpwe_build_system(...)` 这种名字，会先在**当前目录**找，再在**已加入路径的目录**里找。
- 本工具箱把函数分散在 `core/`、`tmm/`、`fdtd/`、`topology/` 等文件夹里。如果没把路径加好，运行任何 demo 都会报 **"Unrecognized function or variable"**。

**把本项目加入路径的方法（二选一，建议都做）：**

```matlab
% 在命令窗口，先 cd 到本项目根目录
cd 'd:\spacetime_crystal_code'

% 方法 A（v2 系列，推荐）：
startup_stm       % 加入 core/tmm/fdtd/topology/demos/tests

% 方法 B（v1 系列）：如果你想跑 examples/ 里的旧示例
stm_init          % 加入 core/tmm/fdtd/topology/examples/tests
```

> ⚠️ 记住这个规则：**跑 demos/ 用 `startup_stm`，跑 examples/ 用 `stm_init`**。
> 两个都跑一遍也无妨——它们加的是不同的子目录，不会冲突。详见 [03-methods-and-map.md](03-methods-and-map.md#五v1-vs-v2同一物理量两套名字)。

**验证环境是否就绪：**

```matlab
which stpwe_build_system    % 如果返回完整路径，说明 core/ 已在路径中
test_smoke                  % 跑自检测试
```

看到 `All smoke tests passed.` 就说明环境 OK 了。

---

## 二、向量与矩阵：MATLAB 的"母语"

MATLAB 里**一切数据都是数组**（array）。标量是 1×1 数组，向量是 1×n 或 n×1，矩阵是 m×n。这是它与 Python/C 最大的不同，也是它名字的来源（MATrix LABoratory）。

### 2.1 造数组

```matlab
% 行向量
v = [1 2 3 4];          % 逗号或空格分隔 => 1×4
w = [1, 2, 3, 4];       % 同上

% 列向量
c = [1; 2; 3; 4];       % 分号分隔 => 4×1

% 等差数列（本工具箱大量使用！）
k = linspace(-1, 1, 101);    % 从 -1 到 1 均匀取 101 个点
x = 0 : 0.1 : 1;             % 0, 0.1, 0.2, ..., 1.0
n = 0 : 2*Mtime;             % 0 到 2*Mtime 的整数

% 全零/全一/单位阵
Z = zeros(3, 4);             % 3×4 全零矩阵
O = ones(5, 1);              % 5×1 全一列向量
I = eye(4);                  % 4×4 单位阵
```

### 2.2 索引（下标从 1 开始！）

```matlab
v = [10 20 30 40];
v(2)          % => 20（注意：不是 0 开始，是 1 开始）
v(2:3)        % => 20 30（切片）
v(end)        % => 40（最后一个元素）
v(1:2:end)    % => 10 30（步长 2）

M = [1 2; 3 4];      % 2×2 矩阵：第一行 1 2，第二行 3 4
M(2, 1)       % => 3（第 2 行第 1 列）
M(:, 1)       % => [1; 3]（所有行，第 1 列）—— 冒号 = "全部"
M(1, :)       % => [1 2]
```

> 💡 **物理直觉提醒**：本工具箱里到处都是 `kBar`（归一化波矢）、`omegaBar`（归一化频率）这种向量，循环里常用 `for ib = 1:numel(kBar)` 逐个取标量。你看到 `xx(i)` 就是"取数组里第 i 个标量"。

### 2.3 算术：注意 `.` 的作用

```matlab
A = [1 2; 3 4];
A * A            % 矩阵乘法（线性代数意义的 A·A）
A .* A           % 逐元素乘法（1 4; 9 16）
A .^ 2           % 逐元素平方
exp(A)           % 逐元素 e^A（不是矩阵指数 expm！）
expm(A)          % 矩阵指数（线性代数意义的 e^A）—— tmm 里会见到
```

> **这是新手最常踩的坑**：`*` 和 `.*`、`^` 和 `.^` 完全不同。如果报错 "Matrix dimensions must agree" 或结果怪怪的，先检查是不是少了点号。

### 2.4 常用数学函数

```matlab
abs(x)     % 模/绝对值
angle(x)   % 复数的辐角（rad）
conj(x)    % 共轭
real(x), imag(x)   % 实部、虚部
exp, log, sqrt, sin, cos, tan, atan2
norm(v)    % 向量范数
sum(v), prod(v), max(v), min(v)
```

---

## 三、脚本 vs 函数：两种 `.m` 文件

项目里所有文件都以 `.m` 结尾。有两种：

### 3.1 脚本（Script）

- 没有 `function` 关键字，就是一个代码清单，**顺序执行**。
- 使用**当前工作区**的变量，运行后变量留在工作区。你在命令窗口逐行敲的代码，效果就是一个"隐形脚本"。
- ⚠️ 本项目**几乎没有真正的脚本文件**：连 `demos/demo01_reproduce_fig2_stpwe.m`、`run_all_demos.m` 都是 `function` 开头、只是**设计成无参数的函数入口**，调用方式跟脚本一样（直接敲函数名即可）。

```matlab
% 在命令窗口直接敲（脚本式用法，变量留在工作区）
p        = stm_fig2_parameters();                  % 参数结构体，见第四节
Nspace   = 10;
Mtime    = 1;
epsCoeff = @(m,n) stm_fig2_eps_coeff(m, n, p);     % 函数句柄，见第五节
sys = stpwe_build_system(epsCoeff, [], Nspace, Mtime, p.g, p.Omega);
```

"脚本式用法"适合在命令窗口临时试东西；写正式代码请用函数，避免变量名冲突。

### 3.2 函数（Function）

- 以 `function [输出] = 函数名(输入)` 开头，必须有 `end` 结尾（或第一行就是函数定义）。
- **有独立的工作区**：函数内部的变量不污染外部，外部变量也进不来（除非显式传入/传出）。
- 例子：`core/stpwe_build_system.m`、`fdtd/fdtd1d_db.m`。

```matlab
function sys = stpwe_build_system(epsCoeff, muCoeff, Nspace, Mtime, g, Omega)
% 第一行是函数签名：输入 6 个东西，输出 1 个结构体
% （下面几百行是函数体，用输入算出 sys）
end
```

> 💡 **识别技巧**：看一个 `.m` 文件第一行有没有 `function` 就知道它是函数还是脚本。函数内部的第一段注释（被 `%` 包围的说明）就是它的"帮助文档"，命令行敲 `help stpwe_build_system` 会显示出来。

### 3.3 命名规则（必须遵守）

- 函数名 = 文件名。`stpwe_build_system` 这个名字必须存在于 `stpwe_build_system.m`。
- MATLAB 名字**区分大小写**。
- 本项目约定：函数以 `stm_`、`stpwe_`、`temporal_`、`fdtd`、`zak_`、`fhs_` 等前缀开头，见 [03-methods-and-map.md](03-methods-and-map.md#三工具箱地图每个目录是干什么的) 和 [docs/tool-reference.md](../docs/tool-reference.md)。

---

## 四、结构体（struct）：本工具箱的数据"集装箱"

科学计算代码里，经常要把一堆相关的量打包传给函数。MATLAB 用**结构体**干这个，字段用 `.` 访问。

```matlab
% 造一个结构体
sys.epsCoeff = [1 0.2];     % 字段 epsCoeff
sys.muCoeff  = [1 0];
sys.G        = diag([-1 0 1]);   % 对角矩阵
sys.size     = 9;

sys            % 命令行看结构体内容
sys.epsCoeff   % 访问单个字段
```

**本工具箱到处是结构体**：
- `stpwe_build_system` **返回**一个结构体 `sys`，装着所有 Fourier 系数矩阵、对角矩阵等；
- `fdtd1d_db` 的**第一个输入**就是一个配置结构体 `cfg`（有几十个字段），**返回**一个输出结构体 `out`；
- `stm_fig2_parameters()` **无输入**，返回一个参数结构体 `p`（字段如 `Lambda`、`g`、`Omega`、`eps1` 等）。demos 里再用 `@(m,n) stm_fig2_eps_coeff(m,n,p)` 把它变成"取 Fourier 系数"的函数句柄。

```matlab
% demo 里的真实用法（示意）：
cfg.x = linspace(0, Lambda, 1000);
cfg.dt = 0.01;
cfg.nSteps = 1000;
cfg.epsFun = @(x) 1 + 0.1*cos(2*pi*x/Lambda);   % 函数句柄（见下一节）
out = fdtd1d_db(cfg);          % 输入一个 cfg，输出一个 out
out.energy; out.probeE;        % 从输出结构体里取数据
```

> 💡 结构体字段名就是"小变量"，`cfg.nSteps`、`cfg.dt` 读作"cfg 的 nSteps 字段"。看到 `out.E` 这种就明白：结果都打包在 `out` 里了。

---

## 五、函数句柄（function handle）：把函数当"值"传递

物理学里，你会说"介质分布是 $\varepsilon(x)$ 这样一个函数"。MATLAB 里想把"函数"本身当参数传出去，就用**函数句柄**，写法是 `@`：

```matlab
% 匿名函数（最常见的写法）
f = @(x) x.^2 + 1;      % 把 "x->x^2+1" 这个函数存进变量 f
f(3)                    % => 10
f([1 2 3])              % => 2 5 10（向量进来向量出去）

% 本项目里真实的样子：
epsFun = @(x) 1 + 0.1*cos(2*pi*x/Lambda);
%        ↑ 输入 x，返回一个数/数组
muFun  = @(x) 1;

cfg.epsFun = epsFun;    % 把函数句柄装进结构体，传给 FDTD
```

> 💡 **为什么需要它？** FDTD 要计算"每一格点的介电常数"。格点是离散的、由网格决定的，FDTD 代码不知道你的材料长什么样。于是你把材料描述（一个函数）作为"参数"传给它——函数句柄就是"函数的指针/引用"。
> 看到 `@(x) ...` 就懂：这是一个一次性定义的小函数，名字不重要，重要的是它的输入输出规则。

---

## 六、控制流：`for` / `if` / `while`

本工具箱的循环很简单，看一遍就会：

```matlab
% for 循环：遍历 k 的每个值
Nk = 101;
for ik = 1:Nk
    k = kBar(ik);                 % 取第 ik 个波矢
    sol = stpwe_solve_omega(sys, k);
    omega(:, ik) = sol.omega;     % 把结果存进矩阵第 ik 列
end

% if 分支
if nargin < 5      % 如果调用时给的输入参数少于 5 个
    Omega = 1;     % 用默认值
end

% while 循环（少用，但 FDTD 里判断"是否稳定"会碰到）
while n < cfg.nSteps
    ...
    n = n + 1;
end
```

**科学计算写循环的习惯**：尽量把向量化运算放在循环外面。比如 `kBar` 整体参与运算用 `.^`、`.*`，只有"每个 k 独立解一次本征值问题"这种无法向量化的才用 `for`。

---

## 七、画图：看懂 demo 的图

Demos 用到的画图命令集中在几个，先记这几个：

```matlab
figure;                       % 新建画布
plot(x, y, 'b-', 'LineWidth', 1.5);   % 蓝色实线
hold on;                      % 在同一张图继续画
plot(x, y2, 'r--');           % 红色虚线
hold off;
xlabel('$k\Lambda/(2\pi)$', 'Interpreter', 'latex');   % 用 LaTeX 写轴标签
ylabel('$\omega/(gc_0)$', 'Interpreter', 'latex');
legend('Im $\omega>0$', 'Im $\omega<0$', 'Interpreter', 'latex');
grid on;
xlim([-0.5 0.5]);             % 设定 x 轴范围
saveas(gcf, 'output/fig.png'); % 存图（很多 demo 存到 output/）
```

其他高频画图命令：`pcolor`（色块，看场分布）、`imagesc`、`colorbar`、`subplot` / `tiledlayout`（多子图）、`semilogy`（对数纵轴，看指数增长）。

> 💡 **看图的物理直觉**：本工具箱最常出现的图是"能带图"——横轴波矢 $k\Lambda/(2\pi)$，纵轴频率 $\omega/(gc_0)$。带隙 = 纵轴方向没有解的空白区域。**复频率图**则看 $\mathrm{Im}\,\omega$ 的正负：正则代表时间上的指数增长（时间带隙/动量带隙的特征）。这些图的意义在 [02-physics-background.md](02-physics-background.md) 详细讲。

---

## 八、调试：遇到报错怎么办

新手遇到报错是必然的，关键是**别慌，按套路走**：

1. **读第一行错误信息**。MATLAB 会告诉你在哪个文件第几行出错：
   ```
   Error in stpwe_build_system (line 42)
   sys.Bomega = ...
   ```
   "line 42" 就是断点——双击打开该文件跳过去。

2. **常见的几类报错**：
   | 报错 | 意思 | 排查方向 |
   |------|------|---------|
   | `Unrecognized function or variable 'xxx'` | 函数不在路径上 / 拼错 | `which xxx`，没找到就 `startup_stm` |
   | `Matrix dimensions must agree` | 数组形状不匹配 | 检查是不是忘写 `.`（`.*` vs `*`）|
   | `Index exceeds the number of array elements` | 下标越界 | 检查数组长度，1 起始 |
   | `Error using ^` | 把 `.^` 写成了 `^` | 数组逐元素运算要加点 |
   | `Not enough input arguments` | 函数调用少传了参数 | 查 [docs/tool-reference.md](../docs/tool-reference.md) 的签名 |

3. **设断点（breakpoint）**：在编辑器点行号左侧的红点，运行到那行会暂停。然后可以用"继续/单步/查看变量"看程序走到哪了、变量是什么值。

4. **随时打印中间量**：在可疑处加 `disp(variableName)` 或 `size(variableName)`，看形状对不对、数值奇不奇怪。

5. **用 `dbstop if error`**：遇到错误自动停在出错现场，比翻命令行历史好用。

---

## 九、科学计算思维：三件最重要的事

学完语法，还有三件比语法更重要的"软件思维"，也是本教程的核心心法：

### 9.1 一切皆离散

连续物理 → 计算机必须离散化：

| 连续量 | 离散化 | 本项目里的体现 |
|--------|--------|---------------|
| 空间坐标 $x$ | 网格 $\Delta x$ | FDTD 的 `cfg.x`，几千个格点 |
| 时间 $t$ | 时间步 $\Delta t$ | FDTD 的 `cfg.dt`、`cfg.nSteps` |
| 周期函数 | Fourier 截断（$N_{space}, M_{time}$） | ST-PWE 的 `Nspace`、`Mtime` |
| 波矢 $k$ | 采样点 `linspace(...)` | 能带图横轴 101 个点 |

**关键心法：离散化引入误差，你必须知道误差有多大。** 这就是"收敛性检查"——把网格加密（`Nspace` 20→40，或 `Δt` 减半），如果结果基本不变，说明你算得够准。

### 9.2 先跑通，再优化

- 第一次运行**用最小参数**（比如 `Nspace=2`、`Mtime=1`、网格很粗），确认"不报错、数量级对"。
- 再逐步加大，直到结果收敛。一次就上大参数，慢、且出错难找。

### 9.3 结果要"对得上"

数值结果必须通过至少一种独立检验：
- **解析检验**：单次时间界面的反射/透射系数 vs 解析的 Morgenthaler 公式（本工具箱有测试）；
- **交叉验证**：FDTD 的结果 vs TMM 的结果（`test_smoke` 里就有一组这样的对比）；
- **守恒量检验**：能量是否守恒、单值矩阵行列式是否 =1（时间 TMM 的 `det(U)=1`）。

> 这些"对得上"的检查，本工具箱已经写进 `tests/test_smoke.m` 了。你以后自己设计新算例，也要养成这个习惯。

---

## 十、动手练习（建议 30–60 分钟）

在命令窗口逐步做，体会语法：

1. **造数组**：`x = linspace(0, 1, 201);` 数一数 `numel(x)` 是不是 201。
2. **逐元素运算**：`y = cos(2*pi*x);` 再 `z = x .* y;`（尝试去掉点号 `x*y` 看报什么错）。
3. **画图**：`plot(x, y, 'LineWidth', 2); grid on;` 加标题 `title('cos')`。
4. **结构体**：`s.eps = 1; s.mu = 1; s` —— 看看工作区里 `s` 变成什么。
5. **函数句柄**：`f = @(x) 1 + 0.2*cos(2*pi*x);` 求 `f(0.25)`。
6. **循环**：写一个 `for i=1:5, disp(i^2), end` 循环。
7. **跑真代码**：`startup_stm;` 然后 `demo01_reproduce_fig2_stpwe`（如果嫌慢，先把 demo 里 `Nspace` 改小再跑）。
8. **查函数**：`help stpwe_build_system` 和 `doc stpwe_build_system` 分别试试。

做完这些，就可以进入下一章，开始理解"你算的到底是什么物理"。

---

## 下一步

- 想理解"时间晶体/时空晶体"的物理图景 → [02-physics-background.md](02-physics-background.md)
- 想先看工具箱整体布局 → [03-methods-and-map.md](03-methods-and-map.md)
- 中途遇到不懂的 MATLAB 语法 → 随时回翻本章，或在 [10-faq-troubleshooting.md](10-faq-troubleshooting.md) 里查 MATLAB 报错的解答（Q1/Q5）。
