# 08 · 12 个演示脚本导读

> 目标读者：已经学完 [04-st-pwe.md](04-st-pwe.md)、[05-temporal-tmm.md](05-temporal-tmm.md)、[06-fdtd.md](06-fdtd.md)、[07-topology.md](07-topology.md) 四章、掌握了三种引擎和拓扑量基本用法的读者。
> 这一章把 `demos/` 目录下的 **12 个论文复现 demo** 逐个带你跑一遍：每个 demo 做什么物理、怎么跑、关键参数在哪改、图怎么读、容易踩什么坑。它不重复讲方法本身——遇到不懂的函数，回对应章节（04/05/06/07），或查 [docs/tool-reference.md](../docs/tool-reference.md) 的函数签名。

读完这一章，你应该能做到：
1. 一条命令跑完所有 demo，也知道每张图存在哪里；
2. 拿到任何一张 `output/*.png`，能说清它画的是什么、验证了什么；
3. 敢动手改参数，把"复现"变成"自己的算例"。

> ⚠️ 一个前提提醒：运行本章任何命令之前，请先确认你已经在 MATLAB 中执行过 `startup_stm`（详见 [01-matlab-primer.md](01-matlab-primer.md) 第一节）。虽然每个 demo 自己也会调用它，但养成先设置路径的习惯，能让你少遇到 `Unrecognized function or variable`。

---

## 一、运行方式：从"跑起来"开始

### 1.1 三种跑法

**方法 A：单跑一个 demo。** 12 个 demo 都是**函数**（不是脚本），函数名就是文件名，在命令窗口直接敲函数名即可：

```matlab
startup_stm                        % 先设置路径（会自动创建 output/ 目录）

demo05_fdtd_temporal_interface     % 跑第 5 个 demo（无参数）
demo01_reproduce_fig2_stpwe('quick')   % demo01 有参数：'quick' 或 'paper'
```

大多数 demo 不带参数，直接函数名回车。带参数的是少数：`demo01_reproduce_fig2_stpwe` 接受 `'quick'`（默认）或 `'paper'`；`run_all_demos` 接受 `true` / `false`。

**方法 B：一条命令跑全部。** 用根目录的 [run_all_demos.m](../run_all_demos.m)：

```matlab
run_all_demos(false)   % 默认（缺省就是 false）：跑 8 个"快"的
run_all_demos(true)    % 再补跑 4 个"重"的：demo01('paper')、demo06、demo07、demo12
```

注意 `run_all_demos` 的参数是一个**逻辑值** `runHeavy`，不是字符串——`true` 才会包含需要几分钟的论文级计算。从 [run_all_demos.m](../run_all_demos.m) 的源码看：

```matlab
if runHeavy
    demo01_reproduce_fig2_stpwe('paper');
    demo06_fdtd_spacetime_wavepacket();
    demo07_zak_phase_stpwe();
    demo12_ptc_convergence_audit();
else
    fprintf('Heavy demonstrations skipped. ...\n');
end
```

也就是说 `run_all_demos(true)` 里，`demo01` 是按 `'paper'` 档跑的（`Nspace=20`、181 个 k 点），比单独跑 `demo01_reproduce_fig2_stpwe`（默认 quick 档）慢得多。

**方法 C：想先看全貌，再逐个深挖。** 先 `run_all_demos(false)` 把快速图全部生成，再按第三节的分组一个个精读。

### 1.2 图存到哪里

每个 demo 内部都调用了 `startup_stm()` 拿到项目根目录 `rootDir`，然后通过 `fullfile(rootDir, 'output', ...)` 把 PNG 存到 **`output/`** 目录（该目录由 `startup_stm` 自动创建）。保存方式有两种，效果一样：

- `demo01` 用一个本地辅助函数 `save_demo_figure(fig, outputFile)`；
- 其余 demo 直接在代码里写 `try exportgraphics(...); catch print(...)`。

两种都是"先试较新的 `exportgraphics`（220 DPI），失败就退回 `print`"。所以你的 MATLAB 版本旧一点也能出图，只是清晰度可能略差。

> 💡 图表名就是 demo 名的小写加下划线，例如 `demo05_fdtd_temporal_interface.png`。个别 demo 还额外存 `.mat` 数据文件（见 `demo12`）。

---

## 二、总览表：12 个 demo 一张表

| demo 号 + 文件名 | 复现哪篇论文 / 哪个物理 | 引擎 | 主要输入参数 | 输出图文件名（`output/` 下） | 运行时间（约） |
|---|---|---|---|---|---|
| 01 `demo01_reproduce_fig2_stpwe` | Park & Min *Opt. Lett.* 46, 484 (2021) Fig.2：时空晶体 Floquet 能带 + 模式场型 | ST-PWE | `quality='quick'/'paper'`（`Nspace`=10/20，`Nk`=101/181），`Mtime=1` | `demo01_fig2_stpwe_quick.png` / `_paper.png` | 几十秒（quick）/ 几分钟（paper） |
| 02 `demo02_complex_frequency_and_momentum_gaps` | 同论文 Eq.(6) vs Eq.(7)：固定 k 复频率（动量带隙）、固定 ω 复波数（频率带隙） | ST-PWE | `kBar∈[0.12,0.23]`×111 点，`fBar∈[0.385,0.425]`×121 点，`Nspace=10,Mtime=1` | `demo02_complex_omega_and_k_gaps.png` | 约 1 分钟 |
| 03 `demo03_ptc_bands_pwe_vs_tmm` | 二元 PTC（ε:1↔4）能带：时间 Fourier PWE 与精确 2×2 TMM 交叉验证 | PWE + 时间 TMM | `epsA=1, epsB=4, dutyA=0.5, T=1, Mtime=19`，`k∈[0.02,1.45]`（181 点） | `demo03_ptc_pwe_vs_tmm.png` | 约 20–40 秒 |
| 04 `demo04_temporal_multilayer_tmm` | Ramaccia *APL* 118, 101901 (2021) 四类时间多层设计（任意堆叠/透明/周期透明/放大） | 时间 TMM | `fNorm∈[0.5,1.5]`（401 点）；4 组 `n`、`dt` | `demo04_temporal_multilayer_tmm.png` | 约几秒 |
| 05 `demo05_fdtd_temporal_interface` | Morgenthaler (1958) 单时间界面的时间折射/反射，数值 vs 解析 | FDTD | `nBefore=1.5, nAfter=2.5, λ0=1, dx=λ0/32`，sponge 吸收 | `demo05_fdtd_temporal_interface.png` | 约 5–15 秒 |
| 06 `demo06_fdtd_spacetime_wavepacket` | Fig.2 时空介质中波包传播，探针谱可见 Floquet 边带（频率转换） | FDTD | 36 元胞、每元胞 40 点，`tEnd=3T`，探针在 `x=14Λ` | `demo06_fdtd_spacetime_wavepacket.png` | 约 0.5–2 分钟 |
| 07 `demo07_zak_phase_stpwe` | 空间 Brillouin 区的双正交 Zak 相位：静态参照 vs 驱动（Floquet）带 | ST-PWE + 拓扑 | `Nk=161`；静态 `Nspace=12,Mtime=0`；驱动 `Nspace=10,Mtime=1` | `demo07_zak_phase_stpwe.png` | 约 1–3 分钟 |
| 08 `demo08_fhs_chern_thouless_pump` | Rice-Mele 泵浦上的 FHS Chern 数（=1），合成维度基准 | 拓扑（紧束缚） | `Nk=Np=61`，`t0=1, δ0=0.6, m0=1` | `demo08_fhs_chern_thouless_pump.png` | 约几秒 |
| 09 `demo09_finite_ptc_order_and_phase` | 有限时长 PTC：AB/BA 元胞无限能带相同、有限输出相位不同 | 时间 TMM | `epsA=3, epsB=1, T=1, nPeriods=8`，`k∈[0.35,0.82]`（401 点） | `demo09_finite_ptc_order_and_phase.png` | 约 10–30 秒 |
| 10 `demo10_temporal_domain_wall_mode` | 时间畴壁局域态：AB\|BA 动量带隙内增长/衰减本征态匹配 | 时间 TMM | `epsA=3, epsB=1, nLeft=nRight=8`，`k∈[0.45,0.72]`（1001 点） | `demo10_temporal_domain_wall_mode.png` | 约 10–30 秒 |
| 11 `demo11_coherent_time_interface` | 双输入相干时间界面：相对相位控制输出，D/B vs E/B 跳变律对比 | 时间 TMM | `ε_before=1.5², ε_after=2.5²`，`phase∈[-π,π]`（721 点） | `demo11_coherent_time_interface.png` | 约几秒 |
| 12 `demo12_ptc_convergence_audit` | 二元 PTC 的 PWE 收敛审计：误差随截断 M 下降 + 运行时成本（存 MAT） | PWE + 时间 TMM | `epsA=1, epsB=4, dutyA=0.5`，`M∈{3,5,9,13,19}`，`k∈[0.25,1.1]`（61 点） | `demo12_ptc_convergence_audit.png` + `demo12_ptc_convergence_data.mat` | 约 1–3 分钟 |

> ⏱️ 时间都是"量级估计"，取决于你的机器。判断一个 demo 是否属于"重"计算：看它是否在 `run_all_demos(true)` 的列表里（01/06/07/12）。

---

## 三、按主题精读

下面把 12 个 demo 按物理主题分成五组。每个 demo 的小节都按同一套结构讲：**做什么物理 → 怎么跑 → 关键参数在哪改 → 看图要点 → 常见坑**。

### 3.1 能带与带隙组（demo01 / demo02 / demo03 / demo12）

这组是"频域"主线：用 ST-PWE 或时间 TMM 把时空调制介质变成本征值问题，算出能带、带隙、增长/衰减率。方法细节见 [04-st-pwe.md](04-st-pwe.md) 和 [05-temporal-tmm.md](05-temporal-tmm.md)。

#### demo01 — 复现 Park-Min Fig.2 的 Floquet 能带与场型

**做什么物理。** 这是全库的"母体" demo。它用时空平面波展开（ST-PWE）求解 Park & Min *Opt. Lett.* 2021 论文 Fig.2 的时空周期介质——介电常数既随空间又随时间周期变化，频率不再是好量子数，能带变成 Floquet 能带，还可能打开**动量带隙**（动量带隙内 $\omega$ 取复值，对应场的时间指数放大/衰减）。它用到的介质是"Fig.2 介质"（[stm_fig2_parameters.m](../core/stm_fig2_parameters.m)）：

```
一个元胞 Λ=1 内：
┌──────────────────────┬──────────────┐
│  静态区 ε1=2         │ 调制区 εc=6  │
│  0 ≤ x < 3Λ/4        │ 3Λ/4≤x<Λ     │
│                      │ ε=6(1+0.6·sin Ωt) │
└──────────────────────┴──────────────┘
```

**怎么跑。** `demo01_reproduce_fig2_stpwe` 缺省就是 `'quick'` 档；想复现论文精度用 `'paper'` 档：

```matlab
demo01_reproduce_fig2_stpwe('quick')   % Nspace=10, Nk=101（约几十秒）
demo01_reproduce_fig2_stpwe('paper')   % Nspace=20, Nk=181（几分钟）
```

**关键参数在哪改。** 参数在脚本开头附近（[demo01_reproduce_fig2_stpwe.m](../demos/demo01_reproduce_fig2_stpwe.m)）：

```matlab
switch lower(quality)
    case 'paper'
        Nspace = 20;  Nk = 181;
    case 'quick'
        Nspace = 10;  Nk = 101;
end
Mtime = 1;               % 时间谐波截断：介质只有 m=0,±1，取 1 就够
fMax = 0.82;             % 只显示 fBar ≤ 0.82 的能带
weightThreshold = 0.035; % m=0 参与度低于此值的模被滤掉
imagTolerance = 3e-3;    % |Im(fBar)| 超过此值的模被滤掉
```

`Nspace`、`Mtime` 是 Fourier 截断，是"收敛"旋钮——调大更准更慢。`fMax`、`weightThreshold`、`imagTolerance` 是**显示筛选**旋钮，只影响画出来多少点，不影响数值本身。

**看图要点。** 输出是 3×4 的 `tiledlayout` 大图：

- 子图 (a)：时空介电常数曲面 $\varepsilon(x,t)$——直观看到"一个元胞后四分之一被时间调制"；
- 子图 (b)：**Floquet 能带结构**。灰色实线是调制深度置 0 时的静态带，灰色虚线是它的 Umklapp 复制（上下平移 $\pm\Omega$、$\pm 2\Omega$）；彩色散点是 ST-PWE 解，**颜色 = m=0 Floquet 扇区参与度**（[stpwe_solve_omega.m](../core/stpwe_solve_omega.m) 返回的 `sol.m0Weight`）。注意这个 demo 的色标是"参与度越高越蓝、越低越粉紫"（`cmapBand` 从 $[1,0,1]$ 渐变到 $[0,0,1]$）。**颜色偏蓝的带才是"物理上干净"的主带**；偏粉的多是多个时间扇区混合的模式；
- 子图 (c)–(h)：六个选定模式（`stpwe_select_mode` 在目标 $(k,\omega)$ 附近挑出来的）的电场时空图，红蓝 = 正负场（[stm_redblue.m](../core/stm_redblue.m)），竖黑线标元胞边界。

**常见坑。** `quality` 只能传 `'quick'` 或 `'paper'`，传别的会直接 `error`。另外 `stpwe_static_bands` 在 $k=0$ 处有一个特殊处理（删除广义本征问题的重复零特征值），这是已知行为，不是 bug——你改静态带数量时看到 $k=0$ 处一条带"跳行"不必惊讶。

#### demo02 — 复频率带隙 vs 复波数带隙（Eq.6 vs Eq.7）

**做什么物理。** 时间周期介质有两种互补的本征值问题，两种带隙画在复平面上长得很不一样：

| 本征值问题 | 固定 | 求 | 带隙里发生什么 |
|---|---|---|---|
| Park-Min Eq.(6) | 实 $k$ | 复 $\omega$ | **动量带隙**：$\mathrm{Im}\,\omega\neq0$，场时间指数增长/衰减（PTC 的参数放大带） |
| Park-Min Eq.(7) | 实 $\omega$ | 复 $k$ | **频率带隙**：$\mathrm{Im}\,k\neq0$，只能以渐逝 Bloch 解存在（对应普通光子晶体禁带） |

**怎么跑。** `demo02_complex_frequency_and_momentum_gaps` 无参数，跑完出一张 2×2 图。

**关键参数在哪改。** 两段扫描区间是核心旋钮（[demo02_complex_frequency_and_momentum_gaps.m](../demos/demo02_complex_frequency_and_momentum_gaps.m)）：

```matlab
kBarSweep = linspace(0.12, 0.23, 111);   % 固定 k 扫描范围（动量带隙附近）
fBarSweep = linspace(0.385, 0.425, 121); % 固定 ω 扫描范围（频率带隙附近）
```

筛选条件也值得看：固定 k 扫描用 `real(f)∈(0.075,0.125)`、`|Im(f)|<0.03`、`m0Weight>0.01`；固定 ω 扫描用 `Re(kBar)∈(-0.5,0.5)`、`|Im(kBar)|<0.12`。

**看图要点。** 四个子图：

- 左上 `Fixed k: momentum gap`：横轴 $k\Lambda/(2\pi)$、纵轴 $\mathrm{Re}\,\omega$，颜色 $=|\mathrm{Im}\,\omega|$。某段 k 上散点出现明显的彩色"山脊"（虚部非零）就是动量带隙；
- 右上 `Parametric growth/decay pair`：$\mathrm{Im}\,\omega$ vs $k$，注意 $\pm\mathrm{Im}\,\omega$ **成对出现在 0 两侧**——增长与衰减互为共轭，是参数谐振的特征；
- 左下 `Fixed ω: frequency gap`：横轴 $\mathrm{Re}\,k$、纵轴 $\omega$，颜色 $=|\mathrm{Im}\,k|$。某段 ω 上散点"断开"、颜色变深就是频率带隙；
- 右下 `Evanescent Bloch solutions`：$\mathrm{Im}\,k$ vs $\omega$，`xline(0.4)` 的竖虚线标出带隙中心 $\bar{f}=0.4$。

**常见坑。** 这张图全是"被筛选过的散点"，所以**空白的区域不一定是没解，可能是被阈值滤掉了**。看到奇怪的空洞，先回去看 `imagTolerance` 一类的阈值。

#### demo03 — 二元 PTC：PWE 与 TMM 交叉验证

**做什么物理。** 光子时间晶体（PTC）是只在时间上周期的介质，最简单的模型是二元方波调制：$\varepsilon(t)$ 在 $1.0$ 与 $4.0$ 之间按占空比 0.5 跳变。同一个准频率能带，用两条独立路线各算一遍并叠在一起对比：

- **TMM（精确参照）**：[temporal_crystal_bands.m](../tmm/temporal_crystal_bands.m) 用 2×2 单周期矩阵 `U`（[temporal_crystal_monodromy.m](../tmm/temporal_crystal_monodromy.m)）的本征值 $\lambda=e^{-i\omega_F T}$ 得到准频率；
- **PWE**：把方波展开成时间 Fourier 级数截断到 `Mtime`，用 `stpwe_build_system` + `stpwe_solve_omega`，再用 [stpwe_fold_frequency.m](../core/stpwe_fold_frequency.m) 把 $\omega$ 折回第一时域 Brillouin 区。

方波不连续、Fourier 系数按 $1/|m|$ 衰减，收敛远慢于正弦调制，所以要取很大的 `Mtime=19`。

**怎么跑。** `demo03_ptc_bands_pwe_vs_tmm` 无参数。

**关键参数在哪改。** 介质参数和截断在开头（[demo03_ptc_bands_pwe_vs_tmm.m](../demos/demo03_ptc_bands_pwe_vs_tmm.m)）：

```matlab
epsA = 1.0;  epsB = 4.0;  dutyA = 0.5;  T = 1;
durations = [dutyA, 1-dutyA]*T;
Omega = 2*pi/T;
Mtime = 19;   % 方波要取很大才收敛；改小体会误差变大
```

PWE 这边因为介质空间均匀，`Nspace=0`；Fourier 系数用 `temporal_binary_eps_coeff(m, epsA, epsB, dutyA)`，再乘以 `double(n==0)`（只有 $n=0$ 的空间谐波）。

**看图要点。** 1×2 布局：

- 左图 `Principal Floquet bands`：黑线 = TMM（精确），彩色散点 = PWE（颜色 = m=0 参与度）。**散点应几乎完全落在黑线上**；若在能带拐弯处偏离，说明 `Mtime` 不够；
- 右图 `Momentum-gap growth/decay`：$\mathrm{Im}\,\omega_F/\Omega$ vs $k$。**$\mathrm{Im}\,\omega_F\neq0$ 的区间就是动量带隙**，PTC 的场在这里指数增长；红圈（PWE）与黑线（TMM）的重合程度就是收敛性证据。

**常见坑。** 脚本后半段（"Figure 2: how the binary PTC is modulated"）会**弹出第二个图形窗口**画 $\varepsilon(t)$ 方波波形（这段代码是激活的），但**保存语句**（`exportgraphics` 到 `demo03_modulation_waveform.png`）被注释掉了——所以 `output/` 里**只有一张 PNG**（`demo03_ptc_pwe_vs_tmm.png`）。别在 `output/` 里找不存在的第二张图。另外，demo03 与 demo12 用同一介质（εA=1, εB=4, duty=0.5），可以连着看：demo03 固定 `Mtime=19` 做交叉验证，demo12 扫描 `Mtime` 做收敛审计。

#### demo12 — PTC 收敛审计（还存档 .mat）

**做什么物理。** 这是"数值方法可信度审计"：以精确 2×2 monodromy 为参照，扫描时间 Fourier 截断 $M \in \{3,5,9,13,19\}$，量化"方波调制下 PWE 的准频率误差随 M 怎么降"，同时记录每个 M 的运行时成本。这是做论文"误差—成本"表格的标准范式，也回答了上一节 demo03 里"为什么 `Mtime` 要取 19"。

**怎么跑。** `demo12_ptc_convergence_audit` 无参数。跑完除了 PNG 还存一个 **`.mat` 数据文件**，把参数、误差、运行时打包在一起：

```matlab
parameters = struct('epsA',epsA,'epsB',epsB,'dutyA',dutyA,'T',T,'kNorm',kNorm);
save(dataFile,'parameters','Mvalues','errors','medianError', ...
    'maximumError','runtime');
```

**关键参数在哪改。** 扫描集合和 k 网格：

```matlab
Mvalues = [3 5 9 13 19];
kNorm = linspace(0.25,1.1,61);
reference = temporal_crystal_bands(kValues,[epsA epsB],[1 1],[dutyA 1-dutyA]*T);
```

`errors(branch,k,M)` 记录两条分支、每个 k、每个 M 的 $|\Delta\omega|/\Omega$（PWE 折频后按最近分支匹配到精确参照）。

**看图要点。** 2×2：

- 左上：`semilogy` 画 median / max 误差 vs M——**误差应随 M 单调下降**。方波调制的下降近似 $1/M$ 量级、不光滑，与正弦调制的指数收敛形成鲜明对比（这正是"方波要取大 Mtime"的量化理由）；
- 右上：每 k 的最坏分支误差热图（纵轴 M、横轴 k）——误差最大的 k 通常落在带隙边界/能带交叉处，可据此判断哪些物理区域需要更高截断；
- 左下：精确参照的 $\mathrm{Im}\,\omega_F/\Omega$（增长/衰减参考）；
- 右下：运行时 vs M——"要花多少时间换精度"。

**常见坑。** 别忘了同时存在 `.mat` 文件——用 `load('output/demo12_ptc_convergence_data.mat')` 可以把 `parameters`、`Mvalues`、`errors` 等变量直接读回工作区。另外误差匹配用 `omitnan` 跳过无解点，所以 `errors` 里会有 `NaN`，自己分析时不要直接 `mean`。

---

### 3.2 时间多层与界面组（demo04 / demo05 / demo11）

这组围绕**时间界面**和**时间多层**：折射率在某个时刻（或一串时刻）突变，波发生"时间折射 + 时间反射"。这是所有 PTC 物理的最小单元。方法细节见 [05-temporal-tmm.md](05-temporal-tmm.md) 和 [06-fdtd.md](06-fdtd.md)。

#### demo04 — 复现 Ramaccia APL 四类时间多层设计

**做什么物理。** "时间多层"（temporal multilayer）= 折射率随时间台阶式跳变的一串空间均匀介质层。每个时间界面用 Morgenthaler 匹配矩阵连接（本工具箱取 D、B 连续，见 [temporal_interface_matrix.m](../tmm/temporal_interface_matrix.m)），层内经历相位积累。通过设计各层折射率 $n_i$ 与持续时间 $dt_i$，可以实现四种设计（对应论文 Table I 与三个派生设计）：

| case | 折射率序列 $n$ | 时长 $dt$ | 设计意图 |
|---|---|---|---|
| 1 任意堆叠 | [2.5 4.4 2.0 1.6] | [0.5 1.0 3.5 0.25] | 复现论文 Table I |
| 2 f₀ 处透明 | [2 7 1.5 5] | `n/(2·na)` ⇒ $\delta=\pi$ | 对 $f_0$ 透明 |
| 3 周期透明 | [6 1.5 6 1.5 6] | `n/(2·na)` ⇒ $\delta=\pi$ | 周期堆叠透明 |
| 4 时间放大 | [6 1.5 6 1.5 6] | `n/(4·na)` ⇒ $\delta=\pi/2$ | 能量注入、放大 |

`dt = n/(2·na)` 使每层内相位 $\delta = \omega_{slab} dt = \pi$（"等距传播时间"），对应透明；`dt = n/(4·na)` 使 $\delta=\pi/2$，对应放大。

**怎么跑。** `demo04_temporal_multilayer_tmm` 无参数，秒级出图。

**关键参数在哪改。** 四个 case 都定义在脚本开头的 `cases` 结构体数组里（[demo04_temporal_multilayer_tmm.m](../demos/demo04_temporal_multilayer_tmm.m)），`cases(q).n` 和 `cases(q).dt` 是核心旋钮。求解循环：

```matlab
spectrum = temporal_tmm_spectrum(omega0, epsInitial, 1, ...
    epsSlabs, ones(size(epsSlabs)), cases(q).dt, epsFinal, 1);
```

其中 `epsInitial=na^2`、`epsSlabs=n.^2`、`epsFinal=nb^2`，第二个和最后一个 `1` 是磁导率（μ≡1）。

**看图要点。** 2×2，每张子图横轴 $f/f_0$：蓝实线 $= |E_b^+/E_a^+|$（前向透射），红虚线 $= |E_b^-/E_a^+|$（后向/时间反射），$f/f_0=1$ 处竖虚线标 $f_0$。判读标准：

- **透明**：$f=f_0$ 处后向系数降到 ~0、前向保持 ~1；
- **放大**：case 4 的前向系数整体 $>1$，说明调制在向场注入能量。

**常见坑。** 注意相位约定的提醒：**Ramaccia 论文沿用自己的 $e^{+i\omega t}$ 相位约定**，与本库其他部分（$e^{ikx-i\omega t}$）不同。比较**幅度**不受影响，但如果你把它的复相位拿去做比较，必须先统一相位约定。本 demo 只画幅度，所以没问题。

#### demo05 — FDTD 模拟单个时间界面（波分裂）

**做什么物理。** 全空间折射率在同一时刻（$t_{switch}=12$）从 $n=1.5$ 突跳到 $2.5$。入射的复高斯波包在时间界面处一分为二：前向波继续前进（速度、波长改变），后向波折返。数值得到的透射/反射幅度要与解析的 Morgenthaler 系数对比。这是把 [06-fdtd.md](06-fdtd.md) 的引擎 `fdtd1d_db` 用起来的最直观算例。

**怎么跑。** `demo05_fdtd_temporal_interface` 无参数。命令行会打印类似：

```
Temporal interface: |tau| exact 0.48000, FDTD 0.48xxx; |rho| exact 0.12000, FDTD 0.12xxx
```

**关键参数在哪改。** `cfg` 结构体是核心（[demo05_fdtd_temporal_interface.m](../demos/demo05_fdtd_temporal_interface.m)）：

```matlab
lambda0 = 1;  k0 = 2*pi/lambda0;
nBefore = 1.5;  nAfter = 2.5;  tSwitch = 12;  tEnd = 30;
dx = lambda0/32;                 % 每波长 32 点
dt = 0.80*dx;                    % CFL 安全（< 1/√2 附近）
x0 = 17;  sigma = 3.5;           % 波包中心与宽度

cfg.x = x;  cfg.dt = dt;  cfg.nSteps = ceil(tEnd/dt);
cfg.epsFun = @(xq,tq) temporal_eps(xq,tq,tSwitch,nBefore^2,nAfter^2);
cfg.boundary = 'sponge';  cfg.spongeCells = 160;  cfg.spongeStrength = 0.08;
cfg.recordEvery = 2;
out = fdtd1d_db(cfg);
```

**看图和数值的对应关系。** 解析系数来自 `[~,tauExact,rhoExact] = temporal_interface_matrix(nBefore^2,1,nAfter^2,1)`（D/B 连续的 Morgenthaler 系数）；数值分解用 `Eplus=0.5*(Eafter + Hafter/nAfter)`、`Eminus=0.5*(Eafter - Hafter/nAfter)` 把界面后的场拆成前向/后向，再算模比 `tauNumeric=norm(Eplus)/norm(Ebefore)`。**这两个数应对到 $10^{-2}$ 量级**——对不上就说明网格不够细或边界吸收没起效。

**看图要点。** 2×2 布局里三块内容（第一块横跨两列）：

- 顶部（横跨两列）：$|E(x,t)|$ 时空图。`yline(tSwitch,'w--')` 白虚线是界面时刻，清晰看到波包一分为二；
- 左下：界面前后两帧的 `real(E)` 叠加，黑/蓝/红分别画 "Before / Forward / Backward"；
- 右下：`out.energy/out.energy(1)` vs 时间——**在 $t_{switch}$ 处能量跳变**，这就是"时间调制对波做功"的直接证据，也是"时间界面产生放大"的种子。

**常见坑。** `cfg.E0` 和 `cfg.Hhalf0` 必须分别在 $t=0$ 与 $t=-dt/2$ 两个交错时刻给出（Yee 蛙跳格式要求），所以 `Hhalf0` 用了 `nBefore*profile(xH + vBefore*dt/2)` 这种带半个时间步时移的构造。你改波包形状时，别只改 `E0` 忘了同步改 `Hhalf0`。

#### demo11 — 双输入相干时间界面（相对相位控制）

**做什么物理。** 让**两束对向传播的相干波**在同一个全局时间开关处碰撞。时间界面把输入态写成二分量 $[E^+;\, E^-]$，界面矩阵 $M=\begin{bmatrix}\tau & \rho \\ \rho & \tau\end{bmatrix}$ 是 2×2。调节两束输入光的相对相位 $\varphi$，可以让某支出射波相干相消（输出幅度极小）或相干增强。更重要的是：**不同微观开关机制给出的 $\tau,\rho$ 不同，相消点也不同**，所以脚本对比了两种理想跳变律——`'DB'`（D、B 连续）与 `'EB'`（E、B 连续）。

**怎么跑。** `demo11_coherent_time_interface` 无参数，秒级出图。命令行打印 DB 与 EB 各自的相消最小幅度和相位位置（例如 `DB minimum ~1e-4 at 0.xxx*pi`）。

**关键参数在哪改。** 界面矩阵来自 [temporal_interface_matrix_jump.m](../tmm/temporal_interface_matrix_jump.m)：

```matlab
[Mdb,tauDB,rhoDB] = temporal_interface_matrix_jump( ...
    epsBefore,muBefore,epsAfter,muAfter,'DB');
[Meb,tauEB,rhoEB] = temporal_interface_matrix_jump( ...
    epsBefore,muBefore,epsAfter,muAfter,'EB');
```

```matlab
phase = linspace(-pi,pi,721);
ratioDB = abs(tauDB/rhoDB);
inputDB = [ones(size(phase)); ratioDB*exp(1i*phase)];
outputDB = Mdb*inputDB;   % 第二路幅度取 |τ/ρ| 以保证能相消
```

**看图要点。** 2×2：

- 左上：`bar` 画 `real(tau)`、`real(rho)` 的对比——**同一材料跳变，不同跳变律给出不同系数**；
- 右上 / 左下：DB / EB 两种律下 $|E^+_{out}|$（蓝实线）、$|E^-_{out}|$（红虚线）vs 相对相位 $/\pi$，黑点标相消点。**对比两张图：谷底相位不同**，说明相消条件依赖开关模型；
- 右下：$\sum|output|^2/\sum|input|^2$ vs 相位——比值部分 $>1$（放大）、部分 $<1$（减幅），说明时间调制净地交换能量。

**常见坑。** "时间界面的物理"第一步就是**显式声明跳变律**。`'EB'` 律对应 `jumpD=epsAfter/epsBefore`（[temporal_interface_matrix_jump.m](../tmm/temporal_interface_matrix_jump.m) 源码里 `'EB'` case），不是"E 连续就是万能"——实际开关电路可能保持电压、注入电荷或空间局域开关，那需要更大的状态空间。这个 demo 之外，别把 `'DB'` 当成普适结论。

---

### 3.3 时域演化组（demo06）

#### demo06 — Fig.2 介质中的时空波包（探针谱可见频移）

**做什么物理。** 这是 demo01 的"时域实况转播"：让一个初始波包在 **Fig.2 时空晶体介质**中真实演化（调制全程存在，`cfg.epsFun = @(xq,tq) stm_fig2_epsilon(xq,tq,p)`），并在 $x=14\Lambda$ 放一个探针。时空晶体不是简单散射波包，而是把它分解到多个 Floquet-Bloch 模式上，产生**频率转换**——探针的短时傅里叶变换（STFT）谱上能清楚看到 $\pm\Omega$ 的 Floquet 边带。

**怎么跑。** `demo06_fdtd_spacetime_wavepacket` 无参数。这是 `run_all_demos(true)` 才跑的"重" demo（约 0.5–2 分钟），因为 FDTD 步数多 + 谱图计算。

**关键参数在哪改。** 网格与波包（[demo06_fdtd_spacetime_wavepacket.m](../demos/demo06_fdtd_spacetime_wavepacket.m)）：

```matlab
cells = 36;  pointsPerCell = 40;
dx = p.Lambda/pointsPerCell;              % Λ/40 = 0.025
dt = 0.55*dx/p.c0;                        % CFL 安全
tEnd = 3*p.T;                             % 调制周期 T=5，共 3 个周期

k0 = 0.175*p.g;                           % 落在 Fig.2 能带里的波矢
x0 = 8*p.Lambda;  sigma = 2.5*p.Lambda;
nEffective = sqrt(0.75*p.eps1 + 0.25*p.epsc);   % 有效折射率
```

探针位置 `[~,probeId] = min(abs(out.x-14*p.Lambda))`，谱图用 [stm_stft.m](../core/stm_stft.m)：`stm_stft(probe, dtRecord, 128, 12, 512)`（窗长 128、跳步 12、FFT 长度 512）。

**看图要点。** 2×2 布局里三块内容（第一块横跨两列）：

- 顶部（横跨两列）：$|E(x,t)|$ 时空演化图，白点线标探针位置——看波包传播、畸变、产生次级波；
- 左下：**探针谱图**。纵轴 $f\Lambda/c\in[0,1]$，颜色 $=20\log_{10}$ 的功率（dB，色标 $-60\sim0$）。主谱线之外出现**上下间隔 $\Delta f=\Omega/2\pi$ 的平行边带脊**，就是 Floquet 边带/频率转换的直接证据；
- 右下：瞬时能量/初始能量——调制与波包的能量交换（振荡或净增）。

**常见坑。** 谱图里 `20*log10(abs(S)/max(abs(S))+1e-8)` 用了 $+10^{-8}$ 防除零；`caxis([-60 0])` 只显示主峰以下 60 dB。如果某条边带太弱看不到，先确认 `k0` 是否落在能带里（用 demo01 的能带图对照），或把 `tEnd` 加长让频谱更清晰。

---

### 3.4 拓扑组（demo07 / demo08）

这组算拓扑不变量。先跑 demo08（干净的紧束缚基准），再回到 demo07（连续时空介质的 Zak 相位），顺序更顺。方法细节见 [07-topology.md](07-topology.md)。

#### demo08 — Rice-Mele 泵浦上的 FHS Chern 数

**做什么物理。** 这是**算法基准**，不是对连续时空介质的 Chern 数复现。它在经典 Thouless 泵浦——一个随时间周期绝热调制的 Rice-Mele 两带模型——上用 FHS（Fukui-Hatsugai-Suzuki）离散 Brillouin 区公式算 Chern 数。调制相位被当作**合成维度**，与动量 $k$ 一起张成 $(k,\varphi)$ 二维环面，Chern 数就是环面上 Berry 曲率的积分。**预期结果 $C=1.000000000000$（12 位精度）**——这是验证拓扑程序取相正确的"金标准"。

**怎么跑。** `demo08_fhs_chern_thouless_pump` 无参数，秒级完成。命令行打印 `Rice-Mele pump FHS Chern number = 1.000000000000.`。

**关键参数在哪改。** 模型在脚本开头（[demo08_fhs_chern_thouless_pump.m](../demos/demo08_fhs_chern_thouless_pump.m)）：

```matlab
Nk = 61;  Np = 61;
kGrid = -pi + (0:Nk-1)*2*pi/Nk;      % 注意：不含重复右端点
phaseGrid = (0:Np-1)*2*pi/Np;
t0 = 1;  delta0 = 0.6;  mass0 = 1.0;
% 每个 (k,φ)：t1=t0+δ0·cosφ, t2=t0−δ0·cosφ, mass=m0·sinφ
offDiagonal = t1 + t2*exp(-1i*k);
H = [mass, offDiagonal; conj(offDiagonal), -mass];
[V,D] = eig(H,'vector');
states(:,1,ik,ip) = V(:,order(1));   % 取基态
```

然后 `[chern,curvature] = fhs_chern_number(states)`。

**看图要点。** 1×3：

- 左：$(k,\varphi)$ 平面的 FHS Berry 曲率热图，标题直接打 `C=...`（应为 1.0000）；
- 中：**展开后的 Zak 相位 $/(2\pi)$ vs 调制相位**——应是一条从 0 到 1 的单调直线（每周期泵浦一个量子）；
- 右：能隙热图（标题打印最小能隙）——相位扫过处能隙保持打开，这是拓扑相变的必要条件。

**常见坑。** [fhs_chern_number.m](../topology/fhs_chern_number.m) 要求 `states` 形状是 `[Nbasis, Noccupied, Nk, Np]`，并且**两个网格方向都闭环、不包含重复端点**。`kGrid = -pi + (0:Nk-1)*2*pi/Nk` 正好满足（最后一格是 $2\pi-\Delta$，不是 $2\pi$）。如果你把 k 网格改成 `-pi:...:pi`，函数会报 link 奇异错误或 Chern 数不再是整数。若 Chern 数出错，先怀疑网格点数和采样一致性。

#### demo07 — 空间 BZ 的双正交 Zak 相位

**做什么物理。** 对 Fig.2 时空晶体的 ST-PWE 能带，沿**空间 Bloch 动量 $k$** 做双正交 Wilson loop，得到 Zak 相位，并对比**静态参照**（`modDepth=0`）与**驱动（Floquet）**两种情形。注意脚本开头的大段注释反复强调：**这是空间方向的 Zak 相位，不是空间均匀 PTC 的"时间 Zak 相位"**（时间 Zak 相位的 Bloch 坐标是准频率，本库 v2 没有把这一量冒充为已完成，见 [07 章 §5.2](07-topology.md#52-demo07算的是空间-bloch-动量的-zak不是时间-zak-相位)）。

**怎么跑。** `demo07_zak_phase_stpwe` 无参数，约 1–3 分钟（`run_all_demos(true)` 才跑）。

**关键参数在哪改。** 两条路径的参数（[demo07_zak_phase_stpwe.m](../demos/demo07_zak_phase_stpwe.m)）：

```matlab
kGrid = -p.g/2 + (0:Nk-1)*p.g/Nk;          % 不含重复端点；Nk=161
omegaSeed = 0.344*p.g*p.c0;

% 静态：modDepth=0, Nspace=12, Mtime=0
optionsStatic.omegaWindow = [0.25 0.70]*p.g*p.c0;
optionsStatic.maxImag = 1e-7*p.g*p.c0;

% 驱动：Nspace=10, Mtime=1
optionsDriven.omegaWindow = [0.15 0.70]*p.g*p.c0;
optionsDriven.maxImag = 0.03*p.g*p.c0;
optionsDriven.minM0Weight = 0.005;
```

```matlab
bandStatic = stpwe_track_band(sysStatic,kGrid,omegaSeed,optionsStatic);
[zakStatic,infoStatic] = zak_phase_biorthogonal( ...
    bandStatic.R,bandStatic.L,sysStatic.Bomega, ...
    stpwe_bz_sewing_matrix(sysStatic));
```

关键诊断是 `info.minimumLinkMagnitude`：若 `infoDriven.minimumLinkMagnitude < 0.1`，脚本会警告"单带 Zak 相位只能当收敛性诊断，不是稳健不变量"。

**看图要点。** 2×2：

- 左上：追踪到的能带 $\mathrm{Re}\,\omega$ vs $k$（黑=静态，蓝=驱动）；
- 右上：驱动带的 $\mathrm{Im}\,\omega$——不为 0 处即动量带隙；
- 左下：m=0 参与度沿 $k$ 的分布；
- 右下：$|\text{Wilson link}|$ vs $k$（黑=静态，蓝=驱动），**标题直接打印 $\mathrm{Zak}/\pi$ 数值**。

判读：若 $\mathrm{Zak}/\pi\approx0$ 或 $\approx1$，是平凡/非平凡二值；若落在中间值，通常意味着**元胞原点不对称**——Fig.2 元胞从一个层边界开始，不处于反演中心，所以不必是 0 或 π（脚本注释已说明）。

**常见坑。** 别把这张图的"Zak 相位"拿去和 PTC 文献里的"时间 Zak 相位"对号入座——坐标系完全不同。另外 Wilson link 的模应处处接近 1；一旦某处大幅塌陷，说明该带与其他带接近/闭合，单带 Zak 相位失去拓扑含义，此时应改用复合子空间 Wilson loop 并报告最小能隙（见 [07 章 §5.4](07-topology.md#54-什么时候单带不变量会失效)）。

---

### 3.5 有限时长与畴壁组（demo09 / demo10）

这组用时间 TMM 研究"有限个数时间周期"和"两种 PTC 交界"的物理。建议按 demo09 → demo10 顺序看，因为它们共用同一套 AB/BA 介质（`epsA=3, epsB=1`）：demo09 是"体响应"，demo10 是"界面模式"。

#### demo09 — 有限时长 PTC：带相同、相位不同

**做什么物理。** 二元 PTC 的 AB 元胞（$\varepsilon:3\to1$）与 BA 元胞（$1\to3$）在时间上互为循环置换，单周期矩阵相似，所以**无限晶体的能带完全相同**（乘法器、准频率、带隙全一样）。但真实实验只有有限个周期（`nPeriods=8`），时间折射/反射输出的**复数相位**取决于元胞起点（时间原点），因此 AB/BA 可以被实验区分。这是一个"无限能带简并 + 有限观测去简并"的教学案例，也提醒你：**别把这个相位差误读成拓扑不变量**（它不是 Zak 不变量）。

**怎么跑。** `demo09_finite_ptc_order_and_phase` 无参数。命令行打印 `AB/BA maximum Floquet-multiplier difference: ~1e-16`——证明能带真的相同。

**关键参数在哪改。** 有限晶体响应用 [temporal_finite_crystal_response.m](../tmm/temporal_finite_crystal_response.m) 计算（[demo09_finite_ptc_order_and_phase.m](../demos/demo09_finite_ptc_order_and_phase.m)）：

```matlab
epsA = 3;  epsB = 1;  T = 1;  durations = [0.5 0.5]*T;  nPeriods = 8;
epsBackground = 2;  muBackground = 1;
kNorm = linspace(0.35,0.82,401);  kValues = kNorm*2*pi/T;

responseAB = temporal_finite_crystal_response(kValues, ...
    [epsA epsB],[muA muB],durations,nPeriods, ...
    epsBackground,muBackground,[1;0]);    % [1;0] = 纯前向输入
```

**看图要点。** 2×2：

- 左上：`Tr(U)/2` vs $k$（能带判别式），$\pm1$ 红虚线标带隙边界——$|\mathrm{Tr}(U)/2|>1$ 处为动量带隙（放大区）；
- 右上：`gain = sqrt(|E⁺|²+|E⁻|²)`（semilogy），AB 蓝实线、BA 红虚线——**两条线几乎重合**（放大相同）；
- 左下：`arg(E⁻/E⁺)/π` vs $k$（时间反射/折射相对相位）——**两条线分叉**；
- 右下：`phaseDifference/π` vs $k$（AB 与 BA 的有序敏感可观测量）。

重点对比右上（重合）与左下/右下（分叉）：**增益不敏感于元胞顺序，相位敏感**——这不是矛盾，正是"无限能带相同、有限观测不同"的体现。

**常见坑。** 相位差用了 `angle(exp(1i*(...)))` 折回 $[-\pi,\pi]$，所以看这张图要理解"包裹"效果。另外 `temporal_finite_crystal_response` 要求 `inputDirectional` 是 2×1 列向量（[1;0] 表示纯前向），漏了方向约定结果会错。

#### demo10 — 时间畴壁局域态

**做什么物理。** 把两个**顺序相反的二元 PTC**（左侧 AB、右侧 BA）拼起来，中间就是"时间畴壁"。两侧体态谱相同（AB、BA 同带），但界面两侧的 Floquet 本征空间"错位"。在共同的动量带隙内，把 **AB 的增长本征态**（$|\lambda|>1$）与 **BA 的衰减本征态**（$|\lambda|<1$）匹配起来，就能得到一个**指数局域在时间界面处的界面态**——包络在界面最大、向两侧指数衰减。

**怎么跑。** `demo10_temporal_domain_wall_mode` 无参数。命令行打印最佳界面态波矢 `k/(2π/T)` 与本征空间失配度。

**关键参数在哪改。** 核心求解（[demo10_temporal_domain_wall_mode.m](../demos/demo10_temporal_domain_wall_mode.m)）：

```matlab
mode = temporal_domain_wall_mode(kValues, ...
    [epsA epsB],[1 1],durations, ...
    [epsB epsA],[1 1],durations,nLeft,nRight);
% nLeft=nRight=8；kNorm=linspace(0.45,0.72,1001)
```

[temporal_domain_wall_mode.m](../tmm/temporal_domain_wall_mode.m) 的流程：每个 $k$ 用 `temporal_crystal_monodromy` 算左右单周期矩阵 → 取 AB 最大 $|\lambda|$（增长）与 BA 最小 $|\lambda|$（衰减）的本征向量 → 失配度 `mismatch=|det([v_grow,v_decay])|` → 用 `fminbnd` 精修 $k$ → 再在 $\pm n$ 个元胞上正向/反向演化出包络 `stateNorm`。

**看图要点。** 2×2：

- 左上：共同动量带隙内 $|\lambda_{grow}|$（蓝）与 $|\lambda_{decay}|$（红虚线）跨过 $|\lambda|=1$（黑虚线），品红虚线标最优 $k$；
- 右上：`semilogy` 画失配度 vs $k$，红圈标最小值（界面态位置）；
- 左下：**界面包络 $\|[D,B]^\top\|$ vs 时间元胞索引**——在界面（索引 0，黑虚线）处峰值、两侧指数衰减，这就是局域界面态的直接证据；
- 右下：ε 层的台阶图——左 8 个元胞 AB、右 8 个 BA，界面在 0 处。

**常见坑。** `temporal_domain_wall_mode` 要求两侧周期相等（内部会检查 `abs(sum(durationsLeft)-sum(durationsRight))`），且要求扫描范围内**确实存在共同动量带隙**——如果 `bestMismatch` 不有限，函数会直接报错 `No common momentum gap was found`。看到这个错，把 `kNorm` 范围扩大或检查介质参数。

---

## 四、建议的初学顺序

从"时间界面的直觉"出发，逐步加难。下面的顺序结合了每章的配套章节：

| 步 | demo | 一句话理由 | 配套章节 |
|---|---|---|---|
| 1 | **demo01** | 全库母体：ST-PWE 能带 + 场型，建立时空晶体图景 | [04-st-pwe.md](04-st-pwe.md) |
| 2 | **demo05** | 最直观：FDTD 看波分裂，数值对解析 | [06-fdtd.md](06-fdtd.md) |
| 3 | **demo03** | 双方法交叉验证 PTC 能带与动量带隙 | [05-temporal-tmm.md](05-temporal-tmm.md) + [04-st-pwe.md](04-st-pwe.md) |
| 4 | **demo06** | FDTD 实况验证 demo01 的频移预言（需要 demo05 的 FDTD 基础） | [06-fdtd.md](06-fdtd.md) |
| 5 | **demo02** | 复频谱补全：动量带隙 vs 频率带隙 | [04-st-pwe.md](04-st-pwe.md) |
| 6 | **demo12** | 收敛审计：为什么方波要取大 Mtime | [04-st-pwe.md](04-st-pwe.md) + [05-temporal-tmm.md](05-temporal-tmm.md) |
| 7 | **demo04** | 时间多层传递函数设计，巩固 TMM | [05-temporal-tmm.md](05-temporal-tmm.md) |
| 8 | **demo11** | 双输入相干控制 + 跳变律辨析 | [05-temporal-tmm.md](05-temporal-tmm.md) |
| 9 | **demo09** | 无限带相同、有限相位不同 | [05-temporal-tmm.md](05-temporal-tmm.md) |
| 10 | **demo10** | 时间畴壁界面态（承接 demo09 的 AB/BA 介质） | [05-temporal-tmm.md](05-temporal-tmm.md) |
| 11 | **demo07** | 空间 Zak 相位（需要双正交概念，建议先有 demo08 铺垫） | [07-topology.md](07-topology.md) |
| 12 | **demo08** | FHS Chern 数"金标准"（概念独立、跑得最快，可提前到第 7 步左右） | [07-topology.md](07-topology.md) |

> 💡 demo08 其实全库最"干净"、跑得最快，把它排在最后只是为了让"拓扑"概念收尾时一鼓作气；想提前了解 Chern 数也可以先跑它再回头做 demo07。

---

## 五、一句话汇总：每个 demo 验证了什么

| demo | 一句话：它验证了什么 |
|---|---|
| 01 | 时空平面波展开能复现 Park-Min Fig.2 的 Floquet 能带与模式场型（颜色 = m=0 扇区参与度） |
| 02 | 固定 k 的复频率给出**动量带隙**，固定 ω 的复波数给出**频率带隙**，两者是互补的本征值问题 |
| 03 | 二元 PTC 的 PWE 散点精确落在精确 TMM 黑线上 → 两套方法互相验证 |
| 04 | 通过设计时间多层的 `n` 与 `dt` 可实现任意/透明/周期透明/放大四种传递函数 |
| 05 | FDTD 的单时间界面波分裂幅度与 Morgenthaler 解析系数吻合，且界面处能量跳变 |
| 06 | 时空介质中波包的探针谱出现间隔 $\Omega$ 的 Floquet 边带 → 频域预言被时域演化证实 |
| 07 | 空间 BZ 的双正交 Wilson loop 给出静态/驱动两种情形下的 Zak 相位，并报告 link 可靠性 |
| 08 | FHS 公式在 Rice-Mele 泵浦上给出精确量子化的 Chern 数 1 → 拓扑程序取相正确 |
| 09 | AB/BA 元胞无限能带相同（λ 差 ~1e-16），但有限时长输出相位不同 → 相位是有序敏感可观测量 |
| 10 | AB 增长态与 BA 衰减态在共同动量带隙内匹配 → 时间界面处存在指数局域的界面态 |
| 11 | 双输入相对相位可控制输出（相干相消/增强），且相消点依赖 D/B 还是 E/B 跳变律 |
| 12 | 方波调制下 PWE 误差随截断 M 下降、同时记录运行时 → 数值方法可信度审计范式 |

---

## 下一步

- 想自己动手改参数、设计新晶体 → [09-exercises.md](09-exercises.md)
- 遇到报错 → [10-faq-troubleshooting.md](10-faq-troubleshooting.md)
- 想看每个函数的输入/输出字典 → [docs/tool-reference.md](../docs/tool-reference.md)
- 想知道"这套代码能复现哪些论文、哪些只是定性覆盖" → [docs/paper-map.md](../docs/paper-map.md)
