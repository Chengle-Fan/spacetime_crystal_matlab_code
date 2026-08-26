# 03 · 三种数值方法与工具箱地图

> 目标读者：已经会一点 MATLAB（至少读过 [01-matlab-primer.md](01-matlab-primer.md)），想弄明白"**这个工具箱到底怎么组织的**"。
> 读完这一章，你应该能做到：
> 1. 说出工具箱的三套数值引擎各自是什么、什么时候用哪个；
> 2. 看到任意一个 `.m` 文件名，能根据前缀猜到它是干什么的；
> 3. 分清 v1 / v2 两套命名，不再被 `stm_preset_modulated_slab` 和 `stm_fig2_parameters` 搞晕；
> 4. 知道该用 `startup_stm` 还是 `stm_init`，以及一键运行脚本怎么用。

这一章基本不写代码，只画地图。具体每个引擎怎么实现，留给 [04-st-pwe.md](04-st-pwe.md)、[05-temporal-tmm.md](05-temporal-tmm.md)、[06-fdtd.md](06-fdtd.md) 三章分别讲。

---

## 一、三种引擎总览

这个工具箱研究的是"介电常数既随空间、又随时间变化"的材料。把 Maxwell 方程喂给计算机，有不止一种吃法。工具箱提供了三套独立的数值引擎，名字和目录如下：

| 引擎 | 目录 | 数学出发点 | 离散化的取舍 | 适用场景 | 主要局限 | 代表性 demo |
|------|------|-----------|-------------|---------|---------|------------|
| **ST-PWE**（时空平面波展开） | `core/` | 把场和介质都展开成时空 Fourier 级数，Maxwell 方程变成**本征值问题**：固定 $k$ 求 $\omega$（广义本征值 $A(k)\Phi=\omega B\Phi$），或固定 $\omega$ 求 $k$（标准本征值 $A'(\omega)\Phi=k\Phi$） | 截断：只保留 $n\in[-N_{\text{space}},N_{\text{space}}]$、$m\in[-M_{\text{time}},M_{\text{time}}]$ 的谐波 | 时空**双周期**介质：能带、复频率/复波数、场型重构 | 需要介质的 Fourier 系数（解析或数值采样）；截断必须收敛；方波这类间断调制收敛慢，需要大 `Mtime` | `demo01`、`demo02`、`demo07` |
| **时间 TMM**（时间传输矩阵） | `tmm/` | 把时间演化切成一段段"常数台阶"，每段内解析求解，段与段之间用 2×2 **传输/单值矩阵**精确拼接 | 时间方向分段常数（每段内是精确的指数，矩阵乘积没有离散误差） | 时间界面、时间多层膜、**二元时间晶体**能带、畴壁态、有限周期响应 | 要求介质的时间演化是分段常数（或可近似为分段常数）；不擅长任意 $x,t$ 连续调制的介质 | `demo03`、`demo04`、`demo09`、`demo10` |
| **FDTD**（时域有限差分） | `fdtd/` | 在空间网格 $\Delta x$、时间步 $\Delta t$ 上直接"蛙跳"推进 Maxwell 方程（D/B-Yee 格式），是**初值问题** | 空间网格 + 时间步，受 CFL 稳定性条件约束 | 波包演化、时域波形、任意时空变化介质、时域交叉验证 | 数值色散与稳定性要盯着；时间调制会与场交换能量，总能量一般**不守恒**；负介电常数/色散材料需要额外机制 | `demo05`、`demo06`、`demo11`、`demo12` |

> 一句话记忆：
> - 想知道"**能带长什么样**" → ST-PWE；
> - 介质是"**时间上的台阶**" → TMM，而且精确；
> - 想看"**波包随时间怎么跑**" → FDTD。

三个引擎都遵守同一套物理约定：时间因子 $e^{ikx-i\omega t}$、$\varepsilon_0=\mu_0=c_0=1$、$\Lambda=1$，归一化波矢 $\bar k = k\Lambda/(2\pi)$、归一化频率 $\bar f = \omega/(g c_0)$。这些约定的完整清单见 [02 章](02-physics-background.md)，做研究前必读。

---

## 二、选型决策树

新手最常问的一句话是："我这个介质，该用哪个引擎？" 下面这张图帮你做决定：

```
          你的介质在空间/时间上怎么变化？
                        │
   ┌────────────────────┼───────────────────────────┐
   ▼                    ▼                           ▼
 时空都周期            只有时间周期、空间均匀        （其它情况）
（时空晶体）            （光子时间晶体 PTC）          要做时域演化 / 任意介质
   │                    │                           │
   ▼                    ▼                           ▼
  ST-PWE              时间演化是分段常数吗？         FDTD
 (core/)              ├── 是 ──→ 时间 TMM（tmm/）   (fdtd/)
   │                  │           精确、最快         │
   │                  └── 否 ──→ ST-PWE 或 FDTD      │
   │                    （连续调制如正弦）           │
   │                                                │
   └──────────────────────┬─────────────────────────┘
                          ▼
              交叉验证：换一个独立方法再算一遍
              （PWE ↔ TMM，或 FDTD ↔ TMM）
```

读法：**先看介质是不是时空双周期**。是，就用 ST-PWE 算能带；然后问"时间演化是不是分段常数台阶"。是，时间 TMM 给出**精确**结果，而且最快。最后，如果你关心的是"波包在时间上怎么演化"，或者介质是任意时空变化的，那就上 FDTD。

**最重要的一条理念：交叉验证。** 任何数值结果都不该只由一种方法"自证清白"。工具箱里到处是这种互验：

- `demo03`（[demo03_ptc_bands_pwe_vs_tmm.m](../demos/demo03_ptc_bands_pwe_vs_tmm.m)）用**时间 Fourier PWE 和精确 2×2 D/B 单值矩阵（TMM）两种独立方法**，算同一个二元时间晶体的准频率能带，然后互相匹配、互相印证；
- `test_smoke` 里也有一组 FDTD vs TMM 的对比（见 [docs/validation.md](../docs/validation.md) 的误差表）。

这也是[本文章节一](#一三种引擎总览)反复强调的：**三类模型不能混用**。空间均匀+时间周期（用 `tmm/temporal_crystal_*`）、时空均周期（用 `core/stpwe_*`）、有限时空演化（用 `fdtd/fdtd1d_db`）是三个不同的数学问题，别把 A 引擎的算例硬塞给 B 引擎。

---

## 三、工具箱地图：每个目录是干什么的

把整个仓库摊开，目录结构长这样（每个目录一行说明）：

```
spacetime_crystal_code/
├── README.md                # 项目主页：学习路径 + 文档导航 + 版本说明
├── startup_stm.m            # v2 入口：把 demos/ 等加入路径（见第六节）
├── stm_init.m               # v1 入口：把 examples/ 等加入路径（见第六节）
├── run_all_demos.m          # 一键运行 12 个 v2 演示脚本
├── stm_run_examples.m       # 一键运行 v1 示例脚本
│
├── core/                    # 引擎 1：ST-PWE（时空平面波展开），18 个函数
├── tmm/                     # 引擎 2：时间传输矩阵，12 个函数
├── fdtd/                    # 引擎 3：一维 D/B-Yee FDTD（fdtd1d_db）
├── topology/                # 拓扑不变量：Zak 相位、Chern 数、能带跟踪
│
├── demos/                   # 12 个论文复现演示脚本（v2，重点学习对象）
├── examples/                # 10 个示例脚本（v1，旧版，demo 的"前身"）
├── tests/                   # 自检测试：test_smoke（v2）和 stm_selftest（v1）
│
├── docs/                    # 参考文档（manual / tool-reference / conventions 等）
├── reproduction/            # 两个论文复现工程（进阶，见下）
│   ├── Observation of temporal reflection/   # Nature Physics 2023：时间反射实验复现
│   └── Topology of photonic crystals/        # Optica 2018：光子时间晶体拓扑复现
│
├── output/                  # 所有脚本保存图片/数据的目录（运行入口时自动创建）
├── Simple_FDTD_NathanZechar/  # 第三方教学 FDTD（1D/2D/3D），可选学习，非本库原创
└── learn/                   # 本教程
```

几个要点：

- **三个引擎目录 `core/`、`tmm/`、`fdtd/` 是"发动机"，`demos/` 和 `examples/` 是"试车场"**。两类目录里**都是函数文件**（以 `function` 开头）；demos/examples 特意设计成无参数入口，直接敲函数名即可一键运行（详见 [01-matlab-primer.md](01-matlab-primer.md) 第三节）。
- **`docs/` 是"参考书"**：陌生函数查 [docs/tool-reference.md](../docs/tool-reference.md)（全部 113 个 `.m` 文件的函数字典），物理约定与模型边界见本系列 [02](02-physics-background.md)、[03](03-methods-and-map.md)、[07](07-topology.md) 章。
- **`reproduction/` 两个子工程是"最高价值区"**：不是演示，而是把两篇顶刊论文的图逐幅复现并自动验收。各自有独立的 `README.md`（[时间反射](<../reproduction/Observation of temporal reflection/README.md>) 和 [PTC 拓扑](<../reproduction/Topology of photonic crystals/README.md>)）和各自的 `run_all_reproductions.m` 调度脚本。它们偏进阶，第一遍学可以跳过。
- **`output/` 是"出图的地方"**：脚本算完的 PNG、`.mat` 都存这里。你每次运行入口脚本，它会自动创建这个目录。

---

## 四、函数命名约定：看前缀猜用途

工具箱几万个代码字符里藏着一条规律：**函数名前缀 = 它属于哪套引擎/哪类工作**。看到名字就大概知道它在干什么，这是快速读懂代码的最重要捷径。

| 前缀 | 主要目录 | 含义 | 例子 |
|------|---------|------|------|
| `stm_` | `core/` | 介质模型与参数（参数结构体、$\varepsilon(x,t)$、时空 Fourier 系数） | `stm_fig2_parameters`、`stm_preset_modulated_slab`、`stm_fig2_eps_coeff`、`stm_interval_fourier` |
| `stpwe_` | `core/` | 时空平面波展开（建系统矩阵、求解、选模、重构场） | `stpwe_build_system`、`stpwe_solve_omega`、`stpwe_solve_k`、`stpwe_reconstruct_field` |
| `temporal_` | `tmm/` | 时间传输矩阵（界面、晶体能带、多层、畴壁） | `temporal_interface_matrix`、`temporal_crystal_bands`、`temporal_multilayer_tmm`、`temporal_domain_wall_mode` |
| `fdtd1d_` | `fdtd/` | 一维时域差分 | `fdtd1d_db` |
| `zak_`、`fhs_`、`stpwe_bz_` | `topology/` | 拓扑不变量与能带工具 | `zak_phase_biorthogonal`、`fhs_chern_number`、`stpwe_bz_sewing_matrix` |

两个补充说明：

1. **`stm_` 里还有一小撮"工具函数"**，比如 `stm_redblue`（配色）、`stm_stft`（短时傅里叶变换）。它们不建模也不算物理，只是画图和分析的帮手，看到别慌。
2. **`stm_fig2_*` 是"论文 Figure 2 模型的专用函数"**：`stm_fig2_parameters`（参数）、`stm_fig2_eps_coeff`（$\varepsilon$ 的时空 Fourier 系数）、`stm_fig2_epsilon`（实空间 $\varepsilon(x,t)$）。它们精确对应 Park & Min 论文的 Figure 2 算例。下一节你会看到，它们其实是 v1 老函数换了个 v2 名字。

完整的函数清单（每个函数的输入、输出、依赖）见 [docs/tool-reference.md](../docs/tool-reference.md)。

---

## 五、v1 vs v2：同一物理量，两套名字

这一节是本章最重要的内容，因为它是新手最容易踩的坑。

### 5.1 为什么会存在两套？

简单说：**代码是分两个版本写的，后来合并进了同一个仓库**。仓库的 [README.md](../README.md)「版本说明」写得很清楚：

- **v1（原始版）**：根目录 + `examples/`，参数封装用 `stm_preset_modulated_slab` 等；
- **v2（增强版）**：新增 `demos/` 及 12 个论文复现演示脚本，扩展了 `core/`、`tmm/`，参数封装统一改成 `stm_fig2_parameters`，还加了 v1 没有的新功能（有限 PTC、时间畴壁、相干时间界面、收敛审计）。

所以 `examples/` 里的脚本在 v2 里**几乎都有一一对应的 `demos/demoXX_*`**（比如 `example_ptc_pwe_vs_tmm` ↔ `demo03_ptc_bands_pwe_vs_tmm`），只是换了一套函数名。两套并存，物理上是同一个东西。

### 5.2 对照表

下面这张表请收藏。遇到不认识的名字，先对号入座：

| 物理量 | v1 命名 | v2 命名 | 说明 |
|--------|---------|---------|------|
| 参数结构体 $p$ | `stm_preset_modulated_slab()` | `stm_fig2_parameters()` | **返回完全相同的数值参数，仅命名/用途不同**（这是两个源文件自己的文档备注，见 [core/stm_fig2_parameters.m](../core/stm_fig2_parameters.m) 与 [core/stm_preset_modulated_slab.m](../core/stm_preset_modulated_slab.m)） |
| $\varepsilon$ 的时空 Fourier 系数 $\varepsilon_{m,n}$ | `stm_fourier_modulated_slab(m,n,p)` | `stm_fig2_eps_coeff(m,n,p)` | 同一套解析公式（`m=0` 直流、`m=±1` 基频），只改了名字 |
| 实空间介电常数 $\varepsilon(x,t)$ | `stm_permittivity_modulated_slab(x,t,p)` | `stm_fig2_epsilon(x,t,p)` | 同上 |
| 示例/演示脚本 | `examples/`（10 个） | `demos/`（12 个） | v1 旧脚本 vs v2 新脚本 |
| 入口 | `stm_init` | `startup_stm` | 见第六节 |
| 自检 | `stm_selftest` | `test_smoke` | v2 的自检覆盖更广 |

这两组参数函数现在都在 `core/` 目录里，**返回的是同一组数值**。你可以亲自验证——先随便跑 `startup_stm` 或 `stm_init`（两个入口都会把 `core/` 加进路径），然后在命令窗口：

```matlab
% v2 命名
p2 = stm_fig2_parameters();
% v1 命名
p1 = stm_preset_modulated_slab();

% 逐个字段比较，完全一样
p2.Lambda, p1.Lambda
p2.epsc,  p1.epsc
p2.OmegaBar, p1.OmegaBar

% Fourier 系数也一样
c2 = stm_fig2_eps_coeff(1, 2, p2);      % v2
c1 = stm_fourier_modulated_slab(1, 2, p1); % v1
c2 - c1     % 结果是 0（或 1e-17 量级的浮点噪声）
```

### 5.3 我该用哪套命名？

**建议统一用 v2 命名**（`startup_stm` + `stm_fig2_*`），因为：
- v2 的演示脚本（`demos/`）是重点学习对象，也是 `test_smoke` 用的命名；
- v2 功能更全、文档更完整。

遇到 v1 名字（`stm_preset_modulated_slab`、`stm_fourier_modulated_slab`），知道它是 v2 的同物别名就行。另外提醒一句：v1 时代的示例都还用着 v1 名字，所以你在参考书里看到 `stm_preset_modulated_slab()` 时，别以为它和 v2 是两回事。

> 铁律：**遇到不认识或记不清的函数名，先查 [docs/tool-reference.md](../docs/tool-reference.md)**。它按目录列了全部 113 个 `.m` 文件的签名和说明，是工具箱的"字典"。

---

## 六、两个入口：`stm_init` vs `startup_stm`

工具箱有**两个"开机脚本"**，它们干的事几乎一样，但**加进路径的目录不同**。这是新手 100% 会遇到的坑，我逐行读过两个文件，把它们的差别讲清楚。

### 6.1 它们各自加了哪些目录

**`startup_stm.m`（v2 入口）** 把以下目录加进 MATLAB 搜索路径：

```matlab
addpath(rootDir);
addpath(fullfile(rootDir, 'core'));
addpath(fullfile(rootDir, 'tmm'));
addpath(fullfile(rootDir, 'fdtd'));
addpath(fullfile(rootDir, 'topology'));
addpath(fullfile(rootDir, 'demos'));   % ← v2 特有
addpath(fullfile(rootDir, 'tests'));
```

**`stm_init.m`（v1 入口）** 加的是：

```matlab
addpath(rootDir);
addpath(fullfile(rootDir, 'core'));
addpath(fullfile(rootDir, 'tmm'));
addpath(fullfile(rootDir, 'fdtd'));
addpath(fullfile(rootDir, 'topology'));
addpath(fullfile(rootDir, 'examples'));   % ← v1 特有
addpath(fullfile(rootDir, 'tests'));
```

差别就一行：**`startup_stm` 加 `demos/`，`stm_init` 加 `examples/`**。其余 `core/`、`tmm/`、`fdtd/`、`topology/`、`tests/` 两个入口都加。它们还会自动创建 `output/` 目录（见 [startup_stm.m](../startup_stm.m) 和 [stm_init.m](../stm_init.m)）。

### 6.2 用哪个？一句话规则

> **跑 `demos/`（v2 演示）→ 用 `startup_stm`；跑 `examples/`（v1 示例）→ 用 `stm_init`。**

为什么必须分清？因为 MATLAB 只在**当前目录 + 已加入路径的目录**里找函数。如果你只跑了 `startup_stm`，然后去敲 `example_complex_gaps`，MATLAB 会报 "Unrecognized function or variable"，因为 `examples/` 不在路径上；反过来只跑 `stm_init` 就找不到 `demo02`。

```matlab
% 推荐流程：在 MATLAB 里 cd 到本仓库根目录后
startup_stm      % v2 系列，跑 demos/ 用这个
test_smoke       % 自检，看到 "All smoke tests passed." 就说明环境 OK
```

两个都跑一遍也没冲突——它们加的是不同子目录。想跑哪个系列，确保对应的入口跑过即可。更详细的路径机制（`which` 验证、报错排查）见 [01-matlab-primer.md](01-matlab-primer.md) 第一节。

---

## 七、一键运行：`run_all_demos` 与 `stm_run_examples`

不想一个一个敲 demo 名字？工具箱准备了两个"总开关"：

### 7.1 v2：`run_all_demos`

```matlab
run_all_demos(false)   % 快速：8 个 demo（demo02,03,04,05,08,09,10,11）
run_all_demos(true)    % 完整：额外追加 demo01('paper')、demo06、demo07、demo12
```

它会先自动调用 `startup_stm`，然后按合理顺序跑 8 个快速 demo；`run_all_demos(true)` 才追加几个耗时的"论文级精度"算例（重算论文图 2、长 FDTD 波包、Zak 相位、收敛审计）。具体清单见 [run_all_demos.m](../run_all_demos.m)。

### 7.2 v1：`stm_run_examples`

```matlab
stm_run_examples(false)   % 快速：5 个示例（complex_gaps, ptc_pwe_vs_tmm, ...）
stm_run_examples(true)    % 完整：追加 paper 级 fig2、fdtd_wavepacket、zak_phase
```

它内部会调用 `stm_init`，所以**用的是 v1 命名**。清单见 [stm_run_examples.m](../stm_run_examples.m)。

### 7.3 图去哪了？

所有脚本算完都把图存进 **`output/`** 目录（入口脚本会自动创建它）。比如 `demo01`（快速模式）存出的 PNG 就叫 `output/demo01_fig2_stpwe_quick.png`（文件名里带 `quick`/`paper` 说明精度档位）。改完参数重跑，同一文件名会被覆盖。

### 7.4 自检

- `test_smoke`（v2）：十多项数值断言，全过打印 `All smoke tests passed.`。它调用 `startup_stm`，用 v2 命名（[tests/test_smoke.m](../tests/test_smoke.m)）。
- `stm_selftest`（v1）：`test_smoke` 的旧版子集，调用 `stm_init`，用 v1 命名（[tests/stm_selftest.m](../tests/stm_selftest.m)）。

环境装好、怀疑"我是不是改坏东西了"的时候，跑一下 `test_smoke` 是最快的安心丸。

---

## 八、建议的代码阅读顺序

工具箱 22 个脚本（12 demo + 10 example）不可能一晚上全读完。按下面的顺序，把"引擎"和"示例"配对着看，效率最高：

```
第一步（建立主线）   demo01 论文图2 能带   ←→  04 章 ST-PWE 精讲
                    （顺便跑 demo02 看两类复谱）

第二步（交叉验证）   demo03 PWE vs TMM    ←→  05 章 时间 TMM 精讲
                    demo05 FDTD 时间界面  ←→  06 章 FDTD 精讲

第三步（拓扑进阶）   demo07 Zak 相位
                    demo08 Chern 数      ←→  07 章 拓扑精讲
                    demo10 时间畴壁

第四步（综合走读）   08 章 演示脚本导读：把 12 个 demo 逐个过一遍
```

对应关系记好：**04 章 ↔ `core/`（ST-PWE）、05 章 ↔ `tmm/`（时间 TMM）、06 章 ↔ `fdtd/`（FDTD）、07 章 ↔ `topology/`（拓扑量）**。每学完一个引擎，就回去看一眼它的 demo，代码和章节互相印证。

v1 的 `examples/` 建议等 v2 熟悉之后再对照着看——它们是 v2 demo 的"前身"，价值主要在理解版本演进（见第五节）。

---

## 九、章节 ↔ 目录 ↔ 引擎对照表

最后把整张地图收进一张表，方便你随时回来定位：

| 学习章节 | 对应目录 | 引擎/主题 | 重点文件 |
|---------|---------|-----------|---------|
| [04-st-pwe.md](04-st-pwe.md) | `core/` | ST-PWE | `stpwe_build_system`、`stpwe_solve_omega`、`stpwe_solve_k` |
| [05-temporal-tmm.md](05-temporal-tmm.md) | `tmm/` | 时间传输矩阵 | `temporal_interface_matrix`、`temporal_crystal_bands` |
| [06-fdtd.md](06-fdtd.md) | `fdtd/` | 时域有限差分 | `fdtd1d_db` |
| [07-topology.md](07-topology.md) | `topology/` | 拓扑不变量 | `zak_phase_biorthogonal`、`fhs_chern_number` |
| [08-demos-walkthrough.md](08-demos-walkthrough.md) | `demos/` | 12 个论文复现 demo | `demo01`–`demo12` |
| （本文） | `examples/`、`docs/` | v1 旧示例 + 参考书 | `stm_init`、[04 章](04-st-pwe.md) [05 章](05-temporal-tmm.md) [06 章](06-fdtd.md)、[docs/tool-reference.md](../docs/tool-reference.md) |

想复习物理背景，先读 [02-physics-background.md](02-physics-background.md)；想自己动手写代码，跳到 [09-exercises.md](09-exercises.md)；跑挂了，看 [10-faq-troubleshooting.md](10-faq-troubleshooting.md)。

---

## 下一步

- 想知道"ST-PWE 到底怎么把 Maxwell 方程变成矩阵" → [04-st-pwe.md](04-st-pwe.md)
- 中途忘了某个函数 → 查 [docs/tool-reference.md](../docs/tool-reference.md)
- 不确定某类介质该不该用某个引擎 → 回看本文章节二。
