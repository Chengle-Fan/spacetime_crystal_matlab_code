# 10 · 常见问题与排错手册

> 目标读者：**正在跑本工具箱、遇到报错或"结果看着不对劲"**的你。
> 这一章是"按症状查"的排错手册：每个问题都按 **症状 → 原因 → 修复 → 怎么验证修复成功** 四步给出。遇到问题先别慌，按 [一、使用说明](#一使用说明) 的索引找到对应问题，再按 [六、一条通用排查流程](#六一条通用排查流程) 兜底。

读完这一章，你应该能做到：
1. 遇到 `Unrecognized function` 时不乱跑入口，先用 `which` 定位；
2. 分清 `startup_stm` 与 `stm_init` 两个入口，知道 v1/v2 两套函数名；
3. 独立判断"复数频率、NaN/Inf、能量漂移、不收敛"是物理还是数值问题；
4. 用 `help`/`doc`/`which` + 文档查清楚任何函数的用法。

---

## 一、使用说明

本手册按主题分四组（路径与启动、数值与数学、物理与解释、工具与流程）。**快速索引表**（按症状查问题编号）：

| 你看到的症状 | 跳到 |
|------|------|
| 报 `Unrecognized function or variable 'xxx'` | [Q1](#q1-报-unrecognized-function-or-variable-xxx)、[Q2](#q2-stm_init-与-startup_stm-有什么区别用混了怎么办)、[Q3](#q3-跑-examples-旧脚本报-stm_preset_modulated_slab-未定义) |
| 分不清 `startup_stm` 和 `stm_init` | [Q2](#q2-stm_init-与-startup_stm-有什么区别用混了怎么办) |
| 遇到不认识的名字 / 函数签名想确认 | [Q4](#q4-v1v2-函数名对照表)、[Q15](#q15-函数签名用法怎么查) |
| 能带图上有复数频率（虚部非零） | [Q5](#q5-能带图上有复数频率虚部非零) |
| 报 `Matrix dimensions must agree` | [Q6](#q6-matrix-dimensions-must-agree) |
| 结果不收敛、能带随 `Nspace` 变 | [Q7](#q7-结果不收敛能带随-nspace-变化) |
| FDTD 出现 NaN/Inf（爆炸） | [Q8](#q8-fdtd-出现-naninf爆炸) |
| FDTD 能量单调涨/降 | [Q9](#q9-fdtd-能量单调涨降) |
| 本征值求解很慢 | [Q10](#q10-本征值求解很慢) |
| 分不清动量带隙 vs 频率带隙 | [Q11](#q11-动量带隙-vs-频率带隙怎么区分) |
| `demo07` 的 Zak 相位不是 0 或 π | [Q12](#q12-demo07-的-zak-相位为什么不是-0-或-π) |
| `demo08` 的 $C=1$ 是时空介质的 Chern 数吗 | [Q13](#q13-demo08-的-c1-是时空介质的-chern-数吗) |
| 时间界面跳跃律 D/B、E/B、E/H 选哪个 | [Q14](#q14-时间界面跳跃律-dbebdheh-选哪个) |
| 图/数据存到哪了 | [Q16](#q16-输出文件在哪) |
| 没有 MATLAB，Octave 能跑吗 | [Q17](#q17-octave-兼容吗) |
| 写论文怎么引用这个项目 | [Q18](#q18-论文里怎么引用) |
| 想加自己的模型 | [Q19](#q19-想加自己的模型) |

> ⚠️ **先记住一件事**：MATLAB 只会在**当前目录 + 已加入搜索路径（path）的目录**里找函数。本项目一半以上的新手报错，根源都是"该在路径上的目录没加进来"。详见 [01-matlab-primer.md](01-matlab-primer.md) 第一节。

---

## 二、路径与启动类

### Q1. 报 `Unrecognized function or variable 'xxx'`

**症状**

```matlab
>> stpwe_build_system(epsCoeff, [], 10, 1, g, Omega)
Unrecognized function or variable 'stpwe_build_system'.
```

**原因**

三种可能，按概率排列：

1. **没运行入口脚本**——`core/` 等目录没进 MATLAB 路径；
2. **入口加错了目录**——比如跑 v2 的 `demos/` 前只跑了 `stm_init`（它加 `examples/` 不加 `demos/`）；
3. **拼写或大小写错了**——MATLAB 函数名**区分大小写**，`stpwe_build_system` 写成 `STPWE_build_system` 就找不到。

**修复**

先用 `which` 定位，别瞎猜：

```matlab
which stpwe_build_system
```

- 如果返回完整路径，比如 `D:\spacetime_crystal_code\core\stpwe_build_system.m`，说明函数在路径上，问题出在**拼写**（检查大小写）；
- 如果返回 `'stpwe_build_system' not found.`，说明目录没进路径，先初始化：

```matlab
cd 'd:\spacetime_crystal_code'   % 切到项目根目录（换成你自己的路径）
startup_stm                       % v2 系列入口（加 demos/）
```

**怎么验证修复成功**

```matlab
which stpwe_build_system   % 应返回完整路径
test_smoke                  % 全部通过会打印 "All smoke tests passed."
```

> 💡 `which xxx` 是排查"找不到函数"的第一工具：它直接回答"这个函数在不在路径上、在哪个文件里"。想跑 v1 的 `examples/` 系列时，入口换成 `stm_init`，见 [Q2](#q2-stm_init-与-startup_stm-有什么区别用混了怎么办)。

---

### Q2. `stm_init` 与 `startup_stm` 有什么区别？用混了怎么办

**症状**

- 跑 `demo01_reproduce_fig2_stpwe` 报 `Unrecognized function or variable 'demo01_reproduce_fig2_stpwe'.`
- 或者跑 `example_dual_sweep_bands` 报 `Unrecognized function or variable 'example_dual_sweep_bands'.`

**原因**

项目有两个版本、两个入口，它们加入的目录**不完全相同**（这是本项目最容易踩的坑）。我直接读了两个源文件核实，区别只有一处：

| 入口文件 | 加入路径的子目录 | 服务对象 |
|------|------|------|
| [startup_stm.m](../startup_stm.m)（v2） | 根目录、`core/`、`tmm/`、`fdtd/`、`topology/`、**`demos/`**、`tests/` | 12 个论文复现 demo（v2） |
| [stm_init.m](../stm_init.m)（v1） | 根目录、`core/`、`tmm/`、`fdtd/`、`topology/`、**`examples/`**、`tests/` | 10 个旧示例脚本（v1） |

两个入口都：加根目录和 `core/`、`tmm/`、`fdtd/`、`topology/`、`tests/`，并且**自动创建 `output/` 目录**；返回值都是项目根目录的绝对路径（`rootDir`）。

**修复**

按"跑哪个系列、用哪个入口"选：

```matlab
startup_stm   % 要跑 demos/ 系列（demo01~demo12）就用这个
stm_init      % 要跑 examples/ 系列（example_xxx）就用这个
```

两个都跑一遍也无妨——它们加的是**不同的子目录**，不会冲突。习惯上建议：v2 是当前主力，默认 `startup_stm`。

**怎么验证修复成功**

```matlab
which demo01_reproduce_fig2_stpwe   % 需要 startup_stm，返回 demos\...m
which example_dual_sweep_bands      % 需要 stm_init，返回 examples\...m
```

---

### Q3. 跑 `examples/` 旧脚本报 `stm_preset_modulated_slab` 未定义

**症状**

```matlab
>> stm_preset_modulated_slab()
Unrecognized function or variable 'stm_preset_modulated_slab'.
```

**原因**

`stm_preset_modulated_slab` 定义在 [core/stm_preset_modulated_slab.m](../core/stm_preset_modulated_slab.m)。报这个错说明 **`core/` 没进 MATLAB 路径**——最常见的场景是：你 `cd` 进了 `examples/` 目录、又只把脚本里的**某几行**复制到新文件里跑（没带上脚本自带的初始化那行），此时 `core/` 不在路径上。

> 💡 其实 `examples/` 里的脚本自己很自觉：比如 [examples/example_dual_sweep_bands.m](../examples/example_dual_sweep_bands.m) 的第 21 行就是 `rootDir = stm_init();`。**完整地跑整个脚本**通常不会出这个错；出错多半是"只跑片段"或"没跑入口"。

**修复**

```matlab
stm_init       % 把 core/ 和 examples/ 都放进路径（v1 名字对应 v1 入口）
```

**怎么验证修复成功**

```matlab
which stm_preset_modulated_slab   % 应返回 core\stm_preset_modulated_slab.m
```

---

### Q4. v1/v2 函数名对照表

**症状**

- 读到 `demos/` 的脚本发现它用 `stm_fig2_parameters`，而 `examples/` 的旧脚本用 `stm_preset_modulated_slab`，怀疑自己记错了名字。

**原因**

同一份物理（"部分调制平板"原胞）在仓库里有**两套函数名并存**，数值完全一样，只是命名不同：v1 名被 `examples/` 和 [tests/stm_selftest.m](../tests/stm_selftest.m) 使用，v2 名被 `demos/` 和 [tests/test_smoke.m](../tests/test_smoke.m) 使用。都放在 `core/`，所以无论哪个入口都能调到。

**修复——对照表**（核心三对，其余函数名本身就没有新旧之分）

| 功能 | v1 名 | v2 名 |
|------|------|------|
| 返回参数结构体 $p$（$\Lambda=1, c_0=1, \varepsilon_1=2, \varepsilon_c=6, \mathrm{modDepth}=0.6, \bar\Omega=0.20$） | `stm_preset_modulated_slab()` | `stm_fig2_parameters()` |
| 计算时空 Fourier 系数 $\varepsilon_{m,n}$ | `stm_fourier_modulated_slab(m,n,p)` | `stm_fig2_eps_coeff(m,n,p)` |
| 返回实空间介电常数 $\varepsilon(x,t)$ | `stm_permittivity_modulated_slab(x,t,p)` | `stm_fig2_epsilon(x,t,p)` |

**怎么验证修复成功**

不认识任何函数名时，统一动作是查 [docs/tool-reference.md](../docs/tool-reference.md)——那是全部 113 个 `.m` 文件的函数字典（签名、输入输出、依赖）。例如在文件里搜 `stm_fig2_parameters`，就能看到它的备注："与 `stm_preset_modulated_slab` 返回相同的数值参数，仅命名/用途不同"。

---

## 三、数值与数学类

### Q5. 能带图上有复数频率（虚部非零）

**症状**

跑 `demo01` 或 `demo02` 后，能带图上某些 $k$ 处的散点颜色很深（虚部大），或者直接看到"本征频率是复数"。

**原因**

**这不是 bug，这是动量带隙（momentum gap）**——光子时间晶体 / 时空晶体的核心特征。按全库约定 $e^{ikx-i\omega t}$：

- $\mathrm{Im}\,\omega>0$：场在时间上**指数增长**；
- $\mathrm{Im}\,\omega<0$：场在时间上**指数衰减**；
- 动量带隙内两者**成对出现**（增长/衰减互为镜像），这就是参数放大（parametric amplification）。

`stpwe_solve_omega` 解的是固定实 $k$ 的广义本征问题，在带隙区间没有纯实数 $\omega$ 解，于是给出复频率。

**怎么读**

```matlab
sol = stpwe_solve_omega(sys, p.g * 0.2);   % 固定 k，解出全部 ω
imag(sol.omega)                            % 看哪些本征频率带虚部
```

- **允许带**：$\mathrm{Im}\,\omega\approx 0$（数值上小于脚本设定的虚部容差，`demo01` 用的量级是 $3\times10^{-3}$），对应普通传播模；
- **带隙**：$|\mathrm{Im}\,\omega|$ 明显非零，且该值**随截断增大趋于稳定**（见 [Q7](#q7-结果不收敛能带随-nspace-变化)）。

**怎么验证修复成功**

把"虚部非零的区间"和 [demos/demo02_complex_frequency_and_momentum_gaps.m](../demos/demo02_complex_frequency_and_momentum_gaps.m) 画出的动量带隙对照；若你的介质是纯空间周期（时间不调制），则不该出现复频率，出现就说明参数或时间谐波截断有问题。

---

### Q6. `Matrix dimensions must agree`

**症状**

```matlab
>> v = [1 2 3]; w = [4 5 6];
>> v * w
Error using  *
Incorrect dimensions for matrix multiplication. ...
   Matrix dimensions must agree.
```

**原因**

把**逐元素运算**（点乘 `.*`、点除 `./`、点幂 `.^`）写成了**矩阵运算**（`*`、`/`、`^`）。这是 MATLAB 新手第一坑，详见 [01-matlab-primer.md](01-matlab-primer.md) 第 2.3 节。本项目里能带计算到处都是数组，比如 `fBar = sol.omega ./ (p.g * p.c0)` 要的是逐元素除，写漏点号就会报错。

**修复**

```matlab
v .* w    % 逐元素乘 => 4 10 18
v ./ w    % 逐元素除 => 0.2500 0.4000 0.5000
v .^ 2    % 逐元素平方 => 1 4 9
```

**怎么验证修复成功**

`size` 检查两边形状一致，再跑一遍不再报错、数值量级合理：

```matlab
size(v)        % 1×3
size(w)        % 1×3
v .* w         % 正常出结果
```

> 💡 报错信息会精确指出在哪一行。看到 `Error using  *` 而它旁边有个变量是数组，先检查是不是该用 `.*`。

---

### Q7. 结果不收敛、能带随 `Nspace` 变化

**症状**

把 `Nspace` 从 10 改成 20，能带形状明显变了；或虚部阈值附近的判读忽大忽小。

**原因**

ST-PWE 是把场和介质都展开成 Fourier 级数再截断：空间保留 $n=-N_{\rm space}:N_{\rm space}$，时间保留 $m=-M_{\rm time}:M_{\rm time}$。**截断不够 → 结果没收敛**。特别要注意时间方向：

- **正弦调制**：时间上只有 $m=0,\pm1$ 三个分量，`Mtime=1` 就够（`demo01` 就是这么取的）；
- **方波（二元）调制**：波形不连续，Fourier 系数按 $1/|m|$ 慢衰减，需要很大的 `Mtime`（`demo03` 用 `Mtime=19`，`demo12` 扫描 $M\in\{3,5,9,13,19\}$ 展示收敛）。

**修复**

按"先跑通、再加密、到收敛"的顺序：

```matlab
% 粗截断，先跑通
sys1 = stpwe_build_system(epsCoeff, [], 10, 1, p.g, p.Omega);
% 翻倍截断，做收敛对比
sys2 = stpwe_build_system(epsCoeff, [], 20, 1, p.g, p.Omega);
```

**判断收敛的实操方法**：对同一批 $k$，分别用粗/细截断求 $\omega$，看最大变化量。频率变化量小于你关心的精度（比如 $10^{-3}\bar f$），就算收敛；否则继续加 `Nspace`（空间）或 `Mtime`（时间）。

**怎么验证修复成功**

两条能带曲线（`Nspace=10` 和 `20`）画在一起基本重合；正式一点的审计流程直接看 [demos/demo12_ptc_convergence_audit.m](../demos/demo12_ptc_convergence_audit.m)——它用精确的 2×2 单值矩阵当参考，量化"误差随 `Mtime` 怎么降"。

---

### Q8. FDTD 出现 NaN/Inf（爆炸）

**症状**

跑 `fdtd1d_db` 后 `out.E` 里全是 NaN 或 ±Inf，能量曲线暴涨。

**原因**

**CFL 条件被违反**。D/B-Yee 显式蛙跳格式要求 Courant 数

$$
S=\frac{v_{\max}\Delta t}{\Delta x}<1,\qquad v_{\max}=\max\frac{1}{\sqrt{\varepsilon_r\mu_r}},
$$

`dt` 取太大就会数值发散。

**修复**

```matlab
cfg.dt = 0.8 * cfg.dx;    % 安全做法：dt 取 dx 的一个分数（demo05 用 0.8、demo06 用 0.55）
out = fdtd1d_db(cfg);
```

- 减小 `cfg.dt`（或等价地增大 `cfg.dx`）；
- 求解器会做**采样 CFL 审计**并返回 `out.sampledCourant`——**判据：它应 < 1**；
- 如果材料有色散/负折射率/有源，采样审计可能返回 `NaN`（跳过检查），此时必须手动指定 `cfg.maxWaveSpeed` 让求解器用你给的最大波速判断稳定性。

**怎么验证修复成功**

```matlab
out.sampledCourant          % 应 < 1
max(abs(out.E(:)))          % 应有限（不再 NaN/Inf）
```

> 💡 时间调制会与场交换能量，动量带隙里场本来就指数增长——**增长 ≠ 爆炸**。区分方法：看是否在极短时间内涨到 NaN/Inf，且违反 CFL 判据；`out.sampledCourant` 是最直接的证据。

---

### Q9. FDTD 能量单调涨/降

**症状**

`out.energy` 曲线一路往上爬（或往下掉），你觉得"不对劲"。

**原因**

先分清**物理**还是**数值**：

1. **物理（最常见）**：时变介质里总能量**本来就不守恒**——调制在向场注入或抽取能量（[02 章 §5.1](02-physics-background.md#51-能量不守恒但物理没被打破) 明确写了这一点）。落在动量带隙里的波包，能量指数上升是正常现象；
2. **数值/边界**：波包跑出计算域前撞上边界，sponge 吸收层不够强，产生反射把能量"困"在域里。可调：`cfg.spongeCells`（默认 `min(80, floor(Nx/8))`）、`cfg.spongeStrength`（默认 `0.12`），或调整总时长 `cfg.nSteps` 让波包不撞边界。

> 💡 注意 `cfg.recordEvery` 只是"每多少步记录一次"的密度设置，不会造成单调漂移；但它太大会让能量曲线看起来"跳"。

**修复**

```matlab
cfg.spongeCells = 160;          % 加宽吸收层（demo05 就用 160）
cfg.spongeStrength = 0.08;      % 或加大强度
cfg.recordEvery = 2;            % 记录密度，视需要调
```

**怎么验证修复成功**

做一个"对照实验"：把时间调制关掉（`epsFun` 不随时间变）。若此时能量不再单调漂移（数值上应近似守恒，`test_smoke` 里均匀周期 FDTD 的能量漂移要求 < 8%），说明你的"单调涨/降"是时变介质的物理效应，不是 bug。

---

### Q10. 本征值求解很慢

**症状**

`demo01_reproduce_fig2_stpwe` 一跑好几分钟，或者自己写的扫描等得心焦。

**原因**

ST-PWE 对每个 $k$ 都要解一个 $2S\times 2S$ 的广义本征问题，其中

$$
S=(2N_{\rm space}+1)(2M_{\rm time}+1).
$$

`Nspace=20, Mtime=1` 时 $S=123$，矩阵 $246\times 246$，每点一次 `eig`；扫 181 个 $k$ 就是 181 次。截断越大、$k$ 点越多，越慢。

**修复**

1. **先跑通再加密**：从 `Nspace=10, Mtime=1`（$S=63$，矩阵 $126\times126$）开始，确认数量级对，再慢慢加；
2. **只扫需要的 $k$**：能带图只在 $k$ 的某个小区间有你要的带，就别把 `kBarSweep` 拉满整个 Brillouin 区，减少 `Nk`；
3. **用现成的 quality 参数**：`demo01` 和 `examples/` 脚本都带 `'quick'` / `'paper'` 档位：

```matlab
demo01_reproduce_fig2_stpwe('quick')   % Nspace=10, Nk=101 —— 先快速跑通
% demo01_reproduce_fig2_stpwe('paper') % Nspace=20, Nk=181 —— 确认方向后再跑重的
```

**怎么验证修复成功**

`tic/toc` 计时对比；并且确认 `'quick'` 与 `'paper'` 在你关心的频率区间结果一致（这正是 [Q7](#q7-结果不收敛能带随-nspace-变化) 的收敛检查）。

---

## 四、物理与解释类

### Q11. 动量带隙 vs 频率带隙怎么区分

**症状**

`demo02` 画出两类带隙，你分不清哪张图对应哪个。

**原因**

两类带隙对应**两个不同的本征值问题**（[04 章 §二](04-st-pwe.md#二两个本征问题用途不同)、[demos/demo02_complex_frequency_and_momentum_gaps.m](../demos/demo02_complex_frequency_and_momentum_gaps.m)）：

| | 动量带隙（momentum gap） | 频率带隙（frequency bandgap） |
|------|------|------|
| 固定什么 | 固定实波矢 $k$ | 固定实频率 $\omega$ |
| 什么变复 | 频率 $\omega$ 取复值（$k$ 实、$\omega$ 复） | 波矢 $k$ 取复值（$k$ 复、$\omega$ 实） |
| 物理来源 | 时间周期性 | 空间周期性 |
| 场的表现 | 时间上指数增长/衰减 | 空间上指数衰减（消逝波） |
| 对应函数 | `stpwe_solve_omega(sys, k)` | `stpwe_solve_k(sys, omega)` |

**记忆口诀**：**谁变成复数，带隙就随谁命名**——$\omega$ 复数 → 动量/时间带隙；$k$ 复数 → 频率/空间带隙。

**怎么验证修复成功**

跑 [demos/demo02_complex_frequency_and_momentum_gaps.m](../demos/demo02_complex_frequency_and_momentum_gaps.m)，左上/右上两图是固定 $k$（动量带隙，$\pm\mathrm{Im}\,\omega$ 成对），左下/右下是固定 $\omega$（频率带隙，$\mathrm{Im}\,k$ 消逝解）；把它们和 [examples/example_dual_sweep_bands.m](../examples/example_dual_sweep_bands.m) 的双扫描叠加图对照，重合处是传播模、分歧处是两类带隙。

---

### Q12. `demo07` 的 Zak 相位为什么不是 0 或 π

**症状**

[Demo07](../demos/demo07_zak_phase_stpwe.m) 打印出的 `Zak/π` 是一个中间值（比如 `0.34`），而不是教科书里常见的 0 或 π。

**原因**

**Zak 相位依赖单胞原点（unit-cell origin）**。`demo07` 的注释原文是：

> "The Zak phase depends on the unit-cell origin. The Fig. 2 cell starts at a layer boundary rather than an inversion center, so it need not appear as exactly 0 or pi in this coordinate convention."

翻译：Fig.2 的原胞从**层边界**（调制区起点，$x=3\Lambda/4$）开始，而不是从**反演中心**开始，所以在这一坐标约定下 Zak 相位不必等于 0 或 π。

另外还有一个**收敛诊断**要一起读：脚本输出 `info.minimumLinkMagnitude`（最小 Wilson link 幅值），若 < 0.1 会警告——单带 Zak 相位这时只能当收敛诊断，不能当稳健拓扑不变量。

**怎么验证修复成功**

把原胞原点平移半个周期（等价地，把材料函数整体平移）重算，Zak 相位会改变——这正说明了"Zak 相位不是坐标无关的绝对量"。读 `demo07` 输出里两个数：`Zak phase`（数值）与 `minimum link`（可信度）。

---

### Q13. `demo08` 的 $C=1$ 是时空介质的 Chern 数吗

**症状**

`demo08` 打印 `Rice-Mele pump FHS Chern number = 1.000000000000`，你以为是"时空晶体的 Chern 数 = 1"。

**原因**

**不是。** `demo08` 是一个**算法基准测试**，`demo08` 注释原文：

> "A Rice-Mele cycle provides a compact, independently quantized validation of the topology routine. The modulation phase is a synthetic dimension. It is not a Chern-number reproduction for demos 01--07's continuum Maxwell space-time medium."

它在一个 Rice–Mele 两带紧束缚模型上验证 [topology/fhs_chern_number.m](../topology/fhs_chern_number.m) 这个通用算法；"调制相位"是**合成维度**，跟 `demo01`–`07` 的连续 Maxwell 时空介质**没有直接关系**。[07 章 §5.1](07-topology.md#51-demo08c1-是-rice-mele-泵浦基准不是时空介质的-chern-数) 也强调了这一点："除非把目标 Maxwell/电路模型的本征态真正送入同一算法，否则这个整数不代表目标时空介质。"

**怎么验证修复成功**

看 `demo08` 的代码：它构造的是 $2\times2$ 哈密顿量 $H=[\text{mass},\ \text{offDiag};\ \text{conj(offDiag)},\ -\text{mass}]$ 的本征态，而不是 `stpwe_solve_omega` 给出的时空介质本征态。因此正确的表述是："`fhs_chern_number` 在 Rice–Mele 泵上通过了测试，$C=1$"。

---

### Q14. 时间界面跳跃律 D/B、E/B、D/H、E/H 选哪个

**症状**

算时间界面的透/反射系数，发现不同文献用不同的"连续量"，不知道本项目默认哪个、怎么换。

**原因**

**取决于你的微观开关机制**。标准"空间均匀、无源"的宏观时间界面常用 $D$、$B$ 连续；但真实开关线路可能保持电压、断开电荷、或注入冲量——连续的量就变了（[05 章 §2.3](05-temporal-tmm.md#23-通用版本-temporal_interface_matrix_jump跳变律是第-5-个参数)）。

本项目把跃迁写成 $D^+=q_D D^-$、$B^+=q_B B^-$，由 [tmm/temporal_interface_matrix_jump.m](../tmm/temporal_interface_matrix_jump.m) 提供：

| `model` 取值 | 含义 |
|------|------|
| `'DB'`（默认） | $D$、$B$ 连续（标准体模型；`temporal_interface_matrix` 就是它的别名） |
| `'EB'` | $E$、$B$ 连续 |
| `'DH'` | $D$、$H$ 连续 |
| `'EH'` | $E$、$H$ 连续 |
| `struct('jumpD', qD, 'jumpB', qB)` | 自定义跳变因子 |

**为什么必须显式声明**：见 [demos/demo11_coherent_time_interface.m](../demos/demo11_coherent_time_interface.m)——同样的材料跳变（$\varepsilon: 1.5^2 \to 2.5^2$），`'DB'` 和 `'EB'` 给出的 $\tau$、$\rho$ 不同，相干相消点也不同。

**怎么验证修复成功**

```matlab
[MMdb, tauDB, rhoDB] = temporal_interface_matrix_jump(1.5^2, 1, 2.5^2, 1, 'DB');
[MMeb, tauEB, rhoEB] = temporal_interface_matrix_jump(1.5^2, 1, 2.5^2, 1, 'EB');
[tauDB rhoDB; tauEB rhoEB]   % 看两个模型给不同系数
```

跑 `demo11`，对比 `'DB'` / `'EB'` 两行系数与两张输出曲线，就明白"开关机制必须先说清楚"。

---

## 五、工具与流程类

### Q15. 函数签名/用法怎么查

**症状**

拿到一个陌生函数（比如 `stpwe_track_band`），不知道输入输出。

**修复**

三个命令 + 两份文档：

```matlab
help  stpwe_track_band   % 命令行直接显示函数头注释（最快）
doc   stpwe_track_band   % 打开文档浏览器，带格式
which stpwe_track_band   % 看它在哪个文件，方便去读源码
```

再深一层：

- 全部 113 个 `.m` 文件的**签名/输入输出/依赖**都在 [docs/tool-reference.md](../docs/tool-reference.md)；
- 每个方法的完整用法在本系列对应章节（04/05/06/07），函数签名速查在 [docs/tool-reference.md](../docs/tool-reference.md)。

**怎么验证修复成功**

`help` 打印出来的第一行就是函数签名，核对输入参数个数与顺序；拿最小例子试跑一次不报 `Not enough input arguments`。

---

### Q16. 输出文件在哪

**症状**

demo/example 跑完了，图"消失了"。

**原因**

所有脚本把图片（和部分数据）保存到项目根目录下的 **`output/`** 目录。这个目录由 `startup_stm` 或 `stm_init` **自动创建**（两个入口里都有 `mkdir(outputDir)`）。比如 `demo01` 存 `output/demo01_fig2_stpwe_quick.png`，`demo12` 额外存 `output/demo12_ptc_convergence_data.mat`。

**修复**

```matlab
rootDir = startup_stm();          % 返回值就是项目根目录
dir(fullfile(rootDir, 'output'))  % 查看输出目录内容
```

**怎么验证修复成功**

在 MATLAB 当前文件夹面板里点进 `output/`，或 `dir` 列出 `.png`/`.mat` 文件，能看到刚跑完的脚本对应文件名。

---

### Q17. Octave 兼容吗

**症状**

没有 MATLAB 许可证，想用免费的 Octave 跑。

**原因/结论**

**部分兼容，不保证全部功能。** 本项目是纯 MATLAB 代码，[learn/README.md](README.md) 里的原话是："需要 MATLAB 许可证（没有的话 Octave 可部分兼容，但**不保证**所有功能可用）。" 代码里一些较新的绘图/导出函数（如 `exportgraphics`、`tiledlayout`）在旧环境不可用时大多有 `try/catch` 回退到 `print`，但整体以 MATLAB 为准。

**怎么验证修复成功**

在 Octave 里跑 [tests/test_smoke.m](../tests/test_smoke.m)，看能否看到 `All smoke tests passed.`；若个别数值断言或绘图失败，说明该部分依赖了 Octave 不兼容的特性。

---

### Q18. 论文里怎么引用

**症状**

要写论文/报告，不知道怎么给这套工具署名。

**原因/结论**

项目本身没有提供专门的"请引用本工具"字符串，但 [../README.md](../README.md) 末尾的"核心参考论文 / References"列了四篇底层文献：

- J. Park and B. Min, *Opt. Lett.* **46**, 484–487 (2021) —— ST-PWE 方法（`demo01`/`demo02` 复现对象）；
- D. Ramaccia, A. Alù, A. Toscano, F. Bilotti, *Appl. Phys. Lett.* **118**, 101901 (2021) —— 时间多层 TMM（`demo04` 复现对象）；
- G. R. Morgenthaler, *IRE Trans. Microwave Theory Tech.* **6**, 167 (1958) —— 时间界面解析解；
- T. Fukui, Y. Hatsugai, H. Suzuki, *J. Phys. Soc. Jpn.* **74**, 1674 (2005) —— FHS Chern 数离散化。

[docs/validation.md](../docs/validation.md)（英文验证记录）说明本仓库数值结果的三层验证情况，可作"可信度背书"引用。

**怎么验证修复成功**

复现了哪篇论文就引哪篇，并在正文注明"数值结果由本工具箱（spacetime_crystal_code）复现"；引用格式照 [../README.md](../README.md) 的参考文献列表抄写。

---

### Q19. 想加自己的模型

**症状**

现有的 Fig.2 介质、二元 PTC 都不够，想定义自己的 $\varepsilon(x,t)$ 并算能带/拓扑。

**修复**

按 `demos/` 脚本的结构照葫芦画瓢（四步走）：

```matlab
function my_new_model()
% 我的新模型：仿照 demos/demoXX 的最小模板
rootDir = startup_stm();            % 1) 初始化路径（v2 系）

% 2) 参数集中一处（"单一事实来源"的习惯）
Lambda = 1;  g = 2*pi/Lambda;
Omega  = 0.4*pi;                    % 你的调制频率

% 3) 材料函数 -> 建系统 -> 算
epsCoeff = @(m,n) my_eps_coeff(m, n, Lambda, Omega);   % 你的 Fourier 系数
sys = stpwe_build_system(epsCoeff, [], 10, 1, g, Omega);
% ... 你的扫描/求解/画图 ...

% 4) 至少用一个独立检查自检，再保存
saveas(gcf, fullfile(rootDir, 'output', 'my_new_model.png'));
end
```

要点：

1. **先读 [03 章 §二](03-methods-and-map.md#二选型决策树)**，确认你的模型属于三类中的哪一类（ST-PWE / 时间 TMM / FDTD），别混用；还想确认哪些物理在代码范围内，可以读 [docs/paper-map.md](../docs/paper-map.md)；
2. 函数名 = 文件名，放新目录时记得 `addpath`；
3. **用 `test_smoke` 的思想自检**：你的结果至少要通过一个独立检查——解析解（如 Morgenthaler 系数）、第二种数值方法（如 PWE vs TMM）、或守恒量（如 `det(U)=1`）。`tests/test_smoke.m` 里每条断言就是这么干的。

**怎么验证修复成功**

新脚本能跑通；并且"结果对得上"——把 `Nspace` 翻倍结果不变（收敛），再和至少一种独立方法对比误差在可接受范围。

---

## 六、一条通用排查流程

记不住上面所有条目没关系，遇到任何问题都可以走下面这条流程（这也是 `docs/` 里所有研究级检查的浓缩版）：

```
遇到报错 / 结果不对劲
        │
        ▼
读错误信息第一行：哪个文件、第几行？
        │
        ▼
是 "Unrecognized function or variable 'xxx'" 吗？
        │
        ├─ 是 ──► which xxx
        │            ├─ not found ──► 跑入口：startup_stm（demos 系）或 stm_init（examples 系）
        │            └─ 找到了 ──► 检查拼写/大小写、参数个数与顺序
        │
        ▼
查 docs/tool-reference.md 确认函数签名与参数顺序
        │
        ▼
是数值/数学类问题吗？
        │
        ├─ 能带随 Nspace 变 ──► 增大 Nspace/Mtime，直到结果不变（Q7）
        ├─ FDTD 出 NaN/Inf ──► 减小 cfg.dt，看 out.sampledCourant < 1（Q8）
        ├─ 能量单调漂移 ──► 区分物理（时变介质不守恒）vs 边界（sponge）（Q9）
        └─ 复频率 / 复数带隙 ──► 先想这是不是动量带隙的物理（Q5、Q11）
        │
        ▼
"结果对得上"吗？至少做一个独立检查：
解析解 / 第二种数值方法 / 守恒量（det(U)=1、能量守恒）
        │
        ▼
过了 ──► 记录参数与误差，跑 test_smoke 确认环境，继续下一个问题
```

> 最后一句忠告：**数值结果必须通过至少一种独立检验才可信**——这是 [docs/validation.md](../docs/validation.md) 里所有验证记录、以及 `tests/test_smoke.m` 里每条断言的共同哲学。把"跑通"和"跑对"分开：跑通只是不报错，跑对才值得写进论文。

---

*教程维护提示：本系列是物理约定的权威来源（时间因子 $e^{ikx-i\omega t}$、归一化单位 $\varepsilon_0=\mu_0=c_0=1$、$\Lambda=1$、$\bar k=k\Lambda/(2\pi)$、$\bar f=\omega/(gc_0)$），编写时请保持各章一致。*
