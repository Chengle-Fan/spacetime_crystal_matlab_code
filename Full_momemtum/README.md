# 论文图 2、图 3 的传输线模版复现

本目录使用 `../Transmission line` 中的 `tl_build_model`、`tl_pwe_*`、
`tl_tmm_bands` 和 `tl_fdtd1d`，尝试复现 arXiv:2604.17408v1 中可由一维
传输线模型计算的面板：图 2(d–j) 和图 3(d–i)。图 2/3(a–c) 是样品照片与
PCB 三维示意图，不是数值求解器输出，因此不在复现范围内。

## 运行

在 MATLAB 中进入本目录并执行：

```matlab
run_reproduction
```

默认执行完整的 PWE/TMM/FDTD 计算，并把 PNG 和 MAT 结果写入 `outputs/`。
第一次检查环境时可以使用快速网格：

```matlab
results = reproduce_figures_2_3(struct( ...
    'quick',true,'runFields',true,'runPhaseMap',true, ...
    'saveOutputs',true));
```

也可以只验证能带求解器：

```matlab
results = reproduce_figures_2_3(struct( ...
    'runFields',false,'runPhaseMap',false,'saveOutputs',false));
```

## 从论文录入的参数

| 图 | 拓扑 | 参数 |
|---|---|---|
| Fig. 2 | SSPP | `w=4 mm`, `C0=19.8 pF`, `DeltaC=2.38 pF`, `f_r0=371 MHz`, `f_m=575/675/725 MHz`, `Vdc=14 V`, `Cvar=12 pF` |
| Fig. 3 | CROW | `w=12 mm`, `L0=9.2 nH`, `C0=27.8 pF`, `DeltaC=5 pF`, `f_m=700 MHz`, `Cblock=200 pF`, `L_s=165 nH`, `f_col=315 MHz`, `Vdc=10.8 V`, `Cvar=17 pF` |

代码沿用模版的 SI 单位和 `exp(i*k*x-i*omega*t)` 约定，因此
`Im(omega)>0` 表示增长。论文把条带间距记作 `w`，但没有另行给出晶格常数；
复现按模版的显式假设取 `a=w`。

## 输出与验证内容

- `figure2_bands.png`：SSPP 静态色散，以及三个调制频率下的复 Floquet 带；
- `figure2_fields.png`：三个调制频率下的有限链电流强度（`|I|^2` 是
  `|Hz|^2` 的电路代理）；
- `figure3_bands.png`：CROW 静态色散、`Q-kappa` 相图和论文所称的
  700 MHz 全动量带隙检验；
- `figure3_fields.png`：三个激励位置的 CROW 谐振支路电流强度；
- `reproduction_results.mat`：参数、带隙范围、PWE/TMM 误差、场增长率及绘图数据。

每个动态能带都用 PWE 和 TMM 作独立交叉核对；SSPP 以密集 PWE 为主，
含有限 `Cblock` 的四状态 CROW 则以不会误选低频寄生分支的 TMM 为主。有限链部分
直接调用 `tl_fdtd1d`，在约 32 ns 开启时间调制并计算到 70 ns。图中的电流
只是局域磁场的电路代理，不能解释为 PCB 三维全波磁场。

## 论文未报告的信息与限制

论文没有给出完整互感、串/并联损耗、源波形/脉宽、端口阻抗、有限样品精确
单元数及泵浦逐单元误差。本复现因此采用以下可追踪假设：`S=0`、`Rs=Gp=R0=0`、
50 ohm 匹配端口、短时高斯电流源、Fig. 2 使用 97 个节点、Fig. 3 使用 33 个
节点。它们只用于验证模版能否定性产生“带隙随泵频扩展”和“CROW 宽动量
不稳定带及强局域化”，不能当作实验 BOM 或逐像素复刻。

论文对 Fig. 3 还存在一个需要保留的建模差异：正文用
`f_col=1/(2*pi*sqrt(L0*C0))`（相当于 `Cblock=Inf`）解释 315 MHz 起点，
但实际又给出与 `L0` 串联的 `Cblock=200 pF`。模版完整有限电容模型预测静态
起点约 336 MHz，而不是 315 MHz；有限 `Cblock` 在给定
`DeltaC=5 pF, f_m=700 MHz` 下把不稳定区扩展到约 `0...0.83 pi/a`，但仍未
覆盖区边界。独立灵敏度试算表明，把幅度提高到约 `6 pF` 才能在同一电路
模型中得到全区不稳定。代码在静态图中同时画出理想无限电容和 200 pF 修正，
并把这一区别写入结果结构，没有通过暗改参数来隐藏它。因此当前验证结论是：
模版稳定、不同求解器互相吻合，并能复现 Fig. 2 的定性趋势与 Fig. 3 的强局域
增长；但仅凭论文已报告参数，不能定量复现其 Fig. 3(f) 的全动量带隙声明。
