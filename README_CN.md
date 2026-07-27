# STM：时空介质科研 MATLAB 工具包

**STM (Space-Time Media)** 是一套用于一维时空周期性介质电磁仿真的 MATLAB 工具包，将以下数值方法统一在一致的符号与归一化体系下：

1. Park 与 Min 的一维时空平面波展开法（ST-PWE），包括固定波数求复准频率、固定频率求复 Bloch 波数、模式重构和论文图 2 复现。
2. Ramaccia 等人的时间多层传输矩阵法（temporal TMM），包括任意时间层、等传播距离设计及周期性光子时间晶体的单周期矩阵。
3. 直接推进 $D/B$ 的一维 Yee-FDTD，用于突变时间界面和连续时空调制中的电磁波演化。
4. 双正交 Wilson 环 Zak 相位，以及把调制相位当作合成维度的 FHS Chern 数。

全部程序均为 MATLAB 源码，不依赖额外工具箱。建议 MATLAB R2020b 或更新版本。

## 1. 快速开始

把 MATLAB 当前目录切换到本文件所在文件夹，然后运行：

```matlab
stm_init
stm_selftest
```

依次运行较快的示例：

```matlab
stm_run_examples(false)
```

运行包括论文质量 PWE 和较长 FDTD 在内的全部示例：

```matlab
stm_run_examples(true)
```

所有图片默认保存在 `output` 文件夹。

## 2. 代码结构

| 文件 | 用途 |
| --- | --- |
| `examples/example_floquet_bands_and_fields.m` | 重构 Park-Min 论文图 2 的介电常数、Floquet 能带和六个模式场 |
| `examples/example_complex_gaps.m` | 对照论文 Eq. (6) 与 Eq. (7)，显示复频率动量带隙和复波数频率带隙 |
| `examples/example_ptc_pwe_vs_tmm.m` | 对同一个二元光子时间晶体分别做 PWE 与精确单周期 TMM，交叉验证 |
| `examples/example_tmm_multilayer.m` | Ramaccia 论文的任意、透明、周期透明和放大时间多层结构 |
| `examples/example_fdtd_interface.m` | 单个突变时间界面的前向/后向波包分裂，并与 Morgenthaler 系数比较 |
| `examples/example_fdtd_wavepacket.m` | 在 Park-Min 图 2 介质中传播波包，绘制 $E(x,t)$、边带谱和能量交换 |
| `examples/example_zak_phase.m` | 连续介质广义本征问题的双正交 Zak 相位与可靠性诊断 |
| `examples/example_chern_pump.m` | Rice-Mele 调制周期的 FHS Chern 数，作为量子化拓扑算法基准 |
| `core/stpwe_build_system.m` | 建立任意 $\epsilon(x,t)$、$\mu(x,t)$ 的时空卷积矩阵 |
| `core/stpwe_solve_omega.m` | 固定 $k$ 求复 $\omega$ |
| `core/stpwe_solve_k.m` | 固定 $\omega$ 求复 $k$ |
| `core/stpwe_sample_fourier_coefficients.m` | 从任意材料函数数值积分得到二维时空 Fourier 系数 |
| `core/stm_preset_modulated_slab.m` | Park-Min 模型的归一化材料参数预设 |
| `core/stm_permittivity_modulated_slab.m` | Park-Min 模型的实空间介电常数分布 |
| `core/stm_fourier_modulated_slab.m` | Park-Min 模型的解析时空 Fourier 系数 |
| `core/stm_interval_fourier.m` | 空间子区间 Fourier 积分工具函数 |
| `tmm/temporal_multilayer_tmm.m` | 任意时间多层结构的总传输矩阵 |
| `tmm/temporal_crystal_bands.m` | 从 $D/B$ 单周期演化矩阵得到时间晶体准频率 |
| `fdtd/fdtd1d_db.m` | 通用一维时变介质 D/B-Yee 求解器 |
| `topology/zak_phase_biorthogonal.m` | 非 Hermitian/广义本征问题的单带 Wilson 环 |
| `topology/fhs_chern_number.m` | 单带或复合占据子空间的规范不变 Chern 数 |

## 3. 统一的场与 Fourier 约定

ST-PWE 使用

$$
\begin{pmatrix}E\\H\end{pmatrix}
=e^{i(kx-\omega t)}
\sum_{m,n}\widetilde{\Psi}_{m,n}
e^{i(ng x-m\Omega t)} .
$$

材料参数写为

$$
p(x,t)=\sum_{m,n}p_{m,n}e^{i(ng x-m\Omega t)} .
$$

因此 Fourier 系数是

$$
p_{m,n}=\frac{1}{\Lambda T}
\int_0^\Lambda\int_0^T
p(x,t)e^{-ingx+im\Omega t}\,dt\,dx .
$$

在此约定下：

- $\operatorname{Im}\omega>0$ 表示 $e^{-i\omega t}$ 随时间增长；
- $\operatorname{Im}\omega<0$ 表示衰减；
- $\operatorname{Im}k>0$ 或 $<0$ 的物理解读取决于空间传播方向和所选边界条件。

Park-Min 图 2 的默认归一化为

$$
\Lambda=c=1,\quad
\bar{k}=\frac{k\Lambda}{2\pi},\quad
\bar{f}=\frac{\omega\Lambda}{2\pi c},\quad
\bar{\Omega}=0.2 .
$$

## 4. ST-PWE 的两种本征问题

`stpwe_solve_omega` 实现固定 $k$ 的广义本征问题

$$
A(k)\Phi=\omega B\Phi ,
$$

适合寻找动量带隙中的复频率和参量增长/衰减模。

`stpwe_solve_k` 实现固定 $\omega$ 的普通本征问题

$$
\mathcal{K}(\omega)\Phi=k\Phi ,
$$

适合寻找频率带隙中的复 Bloch 波数。

运行：

```matlab
example_complex_gaps
```

可同时看到两类复谱。二者不能简单互换：一个描述时间不稳定性，另一个描述空间倏逝。

### 4.1 复现论文图 2

快速截断：

```matlab
example_floquet_bands_and_fields('quick')
```

论文质量截断：

```matlab
example_floquet_bands_and_fields('paper')
```

`paper` 使用 20 阶正负空间谐波和一阶正负时间谐波，即 $41\times3=123$ 个时空 Fourier 状态、246 阶 Maxwell 广义本征矩阵。

原始 `time.m` 中的着色量实际是 $m=0$ Fourier 扇区占比，并不严格等于论文所述的"投影到所有静态本征模后的能量权重"。本代码把它明确标记为 `m=0 sector participation`，避免把近似量误称为精确静态带权重。能带本征值、复带隙和场重构不依赖这一着色近似。

### 4.2 改成任意时空材料

若可以解析写出系数，只需提供函数：

```matlab
epsCoeff = @(m,n) my_epsilon_coefficient(m,n);
muCoeff  = @(m,n) double(m==0 && n==0);
sys = stpwe_build_system(epsCoeff,muCoeff,12,3,g,Omega);
solution = stpwe_solve_omega(sys,k);
```

若材料分布较复杂，可直接数值取系数：

```matlab
Lambda = 1;
T = 5;
epsFun = @(x,t) 3 + 0.4*cos(2*pi*x/Lambda-2*pi*t/T);

table = stpwe_sample_fourier_coefficients( ...
    epsFun,-6:6,-12:12,Lambda,T,512,512);
epsCoeff = @(m,n) stpwe_lookup_coefficient(table,m,n);

sys = stpwe_build_system(epsCoeff,[],12,6, ...
    2*pi/Lambda,2*pi/T);
```

采样区间不包含重复终点，避免 FFT/积分时重复一个周期边界。

### 4.3 使用内置材料预设

工具包内置了 Park-Min 模型作为材料预设，可直接使用：

```matlab
p = stm_preset_modulated_slab();                          % 获取预设参数
epsCoeff = @(m,n) stm_fourier_modulated_slab(m,n,p);     % 获取 Fourier 系数函数句柄
sys = stpwe_build_system(epsCoeff, [], 12, 1, p.g, p.Omega);
```

你也可以参考 `stm_fourier_modulated_slab.m` 的实现，编写自己的材料 Fourier 系数函数。

## 5. 时间传输矩阵与时间晶体能带

Ramaccia 模型假定每个时间层在空间上均匀，并在整个空间同时切换。时间界面保持 $D$ 与 $B$ 连续。代码采用列向量，因此按实际时间顺序写成

$$
\mathrm{TM}=
\mathrm{MM}_{M+1}\mathrm{DM}_{M}\cdots
\mathrm{MM}_{2}\mathrm{DM}_{1}\mathrm{MM}_{1}.
$$

`example_tmm_multilayer` 给出四类频率响应。TMM 部分沿用该论文的 $e^{+i\omega t}$ 符号；ST-PWE/FDTD 使用 $e^{-i\omega t}$。这只会改变相位符号，不改变绘制的幅度、带隙位置或增长率绝对值。

对周期时间晶体，更稳健的方法是直接推进连续的状态 $[D,B]^T$：

$$
\frac{d}{dt}\begin{pmatrix}D\\B\end{pmatrix}
=-ik
\begin{pmatrix}
0&1/\mu(t)\\
1/\epsilon(t)&0
\end{pmatrix}
\begin{pmatrix}D\\B\end{pmatrix}.
$$

若单周期矩阵的本征值为 $\lambda$，则

$$
\lambda=e^{-i\omega_FT},\qquad
\omega_F=\frac{i}{T}\log\lambda .
$$

`example_ptc_pwe_vs_tmm` 会把这一 $2\times2$ 精确单周期结果与时间 Fourier PWE 结果画在同一张图上，是修改代码后最重要的交叉验证之一。

## 6. 时变介质 FDTD

`fdtd1d_db` 更新的是

$$
B^{n+1/2}=B^{n-1/2}-\Delta t\,\partial_xE^n,
$$

$$
D^{n+1}=D^n-\Delta t\,\partial_xH^{n+1/2},
$$

然后在对应时间用

$$
E=D/\epsilon(x,t),\qquad H=B/\mu(x,t)
$$

恢复场。这样在 $\epsilon$ 或 $\mu$ 突变时，$D$ 与 $B$ 不会被错误地重置。

最小配置示例：

```matlab
cfg.x = 0:0.01:20;
cfg.dt = 0.005;
cfg.nSteps = 2000;
cfg.epsFun = @(x,t) (2+0.5*cos(0.4*pi*t))*ones(size(x));
cfg.muFun = @(x,t) ones(size(x));
cfg.E0 = exp(-((cfg.x-5)/1.5).^2).*exp(1i*2*pi*cfg.x);
cfg.boundary = 'sponge';
cfg.recordEvery = 4;

result = fdtd1d_db(cfg);
imagesc(result.x,result.t,abs(result.E));
set(gca,'YDir','normal');
```

`sponge` 是便于研究原型的渐消边界，不是严格 CPML。做高精度反射系数时，应增加空间长度、检查渐消层收敛，或替换为经过验证的 CPML。

### FDTD 稳定性

时间步必须由整个仿真期间的最大相速度决定。归一化一维非色散介质中，保守条件是

$$
\Delta t < \frac{\Delta x}
{\max_{x,t}\left[1/\sqrt{\epsilon_r(x,t)\mu_r(x,t)}\right]} .
$$

若允许 $\epsilon_r<1$、接近零、为负值、具有材料色散或有源增益，当前简单更新式不再充分，需要 ADE、卷积或专门的色散 FDTD。

## 7. 拓扑不变量

### 7.1 ST-PWE 的双正交 Zak 相位

固定 $k$ 的 ST-PWE 是广义且一般非 Hermitian 的本征问题。右、左本征矢满足

$$
A R_j=\omega_jBR_j,\qquad
L_j^\dagger A=\omega_jL_j^\dagger B .
$$

离散 Wilson 链接使用

$$
U_j=
\frac{L_j^\dagger B R_{j+1}}
{|L_j^\dagger B R_{j+1}|}.
$$

Zak 相位为

$$
\gamma=-\arg\prod_jU_j .
$$

`stpwe_bz_sewing_matrix` 负责在 $k=-g/2$ 与 $k=+g/2$ 之间移动空间 Fourier 指标。若忽略这一 sewing，计算结果会依赖端点基底，而不是物理 Bloch 态。

运行：

```matlab
example_zak_phase
```

请同时检查：

- 目标带是否在整个 Brillouin 区与其他带分离；
- 相邻 Wilson 链接是否远离零；
- 增加空间/时间谐波后 Zak 相位是否收敛；
- 改变 $k$ 网格后结果是否稳定；
- 单胞原点是否与要采用的反演中心一致。

Zak 相位本身随单胞原点改变。只有在明确的空间对称性和原点约定下，才可把 0 或 $\pi$ 的量子化直接解释为拓扑分类。

### 7.2 FHS Chern 数

`fhs_chern_number` 的右本征态接收尺寸为

```matlab
[Nbasis, Noccupied, Nk, Nparameter]
```

的数组。Hermitian 问题只需传右本征态；广义或非 Hermitian 问题还可传入同尺寸的左本征态和重叠度量，使链接变为
$\det[L_i^\dagger M R_j]$。它用链接矩阵行列式处理单带或复合占据子空间，并在两个方向施加周期边界。网格中不要重复 $2\pi$ 端点。

`example_chern_pump` 把调制相位作为第二个周期坐标，对 Rice-Mele 泵得到整数 Chern 数 $C=1$。这是将调制相位、延迟相位或其他周期控制参数当作合成维度时可直接复用的模板。

## 8. 科研计算的收敛顺序

建议每次修改模型后依次检查：

1. `stm_selftest` 是否通过。
2. 材料 Fourier 系数是否满足实函数的共轭关系
   $p_{-m,-n}=p_{m,n}^*$。
3. 增大空间截断 `Nspace`，确认主要带和复带隙收敛。
4. 增大时间截断 `Mtime`。时间方波比正弦调制需要更多谐波。
5. 用 `example_ptc_pwe_vs_tmm` 一类的独立 TMM 结果校验时间 PWE。
6. FDTD 依次减半 $\Delta x$ 与 $\Delta t$，检查带边、增长率和散射幅度。
7. 拓扑量同时检查能隙、Wilson 链接、网格和 Fourier 截断。

只看一张"像论文"的彩图不足以证明数值正确；复谱、独立方法交叉验证和收敛测试应一起保留。

## 9. 适用范围与尚未包含的物理

当前工具包针对一维、各向同性、局域、非色散的时空介质。以下问题需要扩展：

- 二维/三维矢量 Maxwell ST-PWE；
- 各向异性、双各向异性或磁电耦合；
- Lorentz/Drude 色散和因果时变色散；
- 严格 CPML 与开放系统准正常模；
- Floquet 异常拓扑相的完整 Rudner winding；
- 简并复合带的非 Abelian 双正交 Wilson 环；
- 非线性调制与泵浦耗尽。

尤其要注意：时间 TMM 不能直接替代空间也周期变化的 ST-PWE；反之，有限 Fourier 截断的 PWE 也不应被当作突变时间层的"精确解"。

## 10. 对应论文

1. J. Park and B. Min, "Spatiotemporal plane wave expansion method for arbitrary space-time periodic photonic media," *Optics Letters* **46**, 484–487 (2021). DOI: 10.1364/OL.411622.
2. D. Ramaccia, A. Alù, A. Toscano, and F. Bilotti, "Temporal multilayer structures for designing higher-order transfer functions using time-varying metamaterials," *Applied Physics Letters* **118**, 101901 (2021). DOI: 10.1063/5.0042567.
