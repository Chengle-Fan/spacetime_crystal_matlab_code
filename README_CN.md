# 时空介质科研 MATLAB 代码包 v2

v2 根据 `Photonic space-time crystal` 资料库的 20 份 PDF（17 项独立工作）重新审查了物理边界、代码覆盖和科研路线。建议先阅读：

- `START_HERE_CN.md`：正确的入门顺序；
- `PAPER_CODE_MAP_CN.md`：20 份文件逐项覆盖矩阵；
- `CONVENTIONS_AND_LIMITS_CN.md`：边界条件、拓扑定义与适用域；
- `RESEARCH_ROADMAP_CN.md`：可以继续形成课题的开发路线。

核心数值骨架把以下工具放在同一套符号和归一化下：

1. Park 与 Min 的一维时空平面波展开法（ST-PWE），包括固定波数求复准频率、固定频率求复 Bloch 波数、模式重构和论文图 2 复现。
2. Ramaccia 等人的时间多层传输矩阵法（temporal TMM），包括任意时间层、等传播距离设计及周期性光子时间晶体的单周期矩阵。
3. 直接推进 $D/B$ 的一维 Yee-FDTD，用于突变时间界面和连续时空调制中的电磁波演化。
4. 双正交 Wilson 环 Zak 相位，以及把调制相位当作合成维度的 FHS Chern 数。

v2 另外增加了有限时间晶体、AB/BA 单胞相位响应、时间畴壁模、双输入相干时间界面、可配置跃迁律、FDTD CFL/探针模式和 PWE 收敛审计。空间 Zak 和 Rice–Mele Chern 示例被明确标为相应模型的算法基准，不再代替光子时间晶体的严格时间 Zak 相。

全部程序均为 MATLAB 源码，不依赖额外工具箱。建议 MATLAB R2020b 或更新版本。

## 1. 快速开始

把 MATLAB 当前目录切换到本文件所在文件夹，然后运行：

```matlab
startup_stm
test_smoke
```

依次运行较快的示例：

```matlab
run_all_demos(false)
```

运行包括论文质量 PWE 和较长 FDTD 在内的全部示例：

```matlab
run_all_demos(true)
```

所有图片默认保存在 `output` 文件夹。

## 2. 代码结构

| 文件 | 用途 |
| --- | --- |
| `demos/demo01_reproduce_fig2_stpwe.m` | 重构 Park-Min 论文图 2 的介电常数、Floquet 能带和六个模式场 |
| `demos/demo02_complex_frequency_and_momentum_gaps.m` | 对照论文 Eq. (6) 与 Eq. (7)，显示复频率动量带隙和复波数频率带隙 |
| `demos/demo03_ptc_bands_pwe_vs_tmm.m` | 对同一个二元光子时间晶体分别做 PWE 与精确单周期 TMM，交叉验证 |
| `demos/demo04_temporal_multilayer_tmm.m` | Ramaccia 论文的任意、透明、周期透明和放大时间多层结构 |
| `demos/demo05_fdtd_temporal_interface.m` | 单个突变时间界面的前向/后向波包分裂，并与 Morgenthaler 系数比较 |
| `demos/demo06_fdtd_spacetime_wavepacket.m` | 在 Park-Min 图 2 介质中传播波包，绘制 $E(x,t)$、边带谱和能量交换 |
| `demos/demo07_zak_phase_stpwe.m` | 连续介质广义本征问题的双正交 Zak 相位与可靠性诊断 |
| `demos/demo08_fhs_chern_thouless_pump.m` | Rice-Mele 调制周期的 FHS Chern 数，作为量子化拓扑算法基准 |
| `demos/demo09_finite_ptc_order_and_phase.m` | 同谱 AB/BA 时间单胞的有限周期输出与复相位 |
| `demos/demo10_temporal_domain_wall_mode.m` | 共同动量带隙中增长/衰减 Floquet 态的时间畴壁匹配 |
| `demos/demo11_coherent_time_interface.m` | 双输入相干碰撞以及 D/B 与 E/B 跃迁律对照 |
| `demos/demo12_ptc_convergence_audit.m` | 二元时间晶体 PWE 对精确单周期矩阵的收敛审计 |
| `core/stpwe_build_system.m` | 建立任意 $\epsilon(x,t)$、$\mu(x,t)$ 的时空卷积矩阵 |
| `core/stpwe_solve_omega.m` | 固定 $k$ 求复 $\omega$ |
| `core/stpwe_solve_k.m` | 固定 $\omega$ 求复 $k$ |
| `core/stpwe_sample_fourier_coefficients.m` | 从任意材料函数数值积分得到二维时空 Fourier 系数 |
| `tmm/temporal_multilayer_tmm.m` | 任意时间多层结构的总传输矩阵 |
| `tmm/temporal_crystal_bands.m` | 从 $D/B$ 单周期演化矩阵得到时间晶体准频率 |
| `tmm/temporal_interface_matrix_jump.m` | 显式 D/B、E/B、D/H、E/H 或自定义时间跃迁律 |
| `tmm/temporal_finite_crystal_response.m` | 有限周期 PTC 的复方向输出 |
| `tmm/temporal_domain_wall_mode.m` | 两个 PTC 的增长/衰减本征子空间匹配 |
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
demo02_complex_frequency_and_momentum_gaps
```

可同时看到两类复谱。二者不能简单互换：一个描述时间不稳定性，另一个描述空间倏逝。

### 4.1 复现论文图 2

快速截断：

```matlab
demo01_reproduce_fig2_stpwe('quick')
```

论文质量截断：

```matlab
demo01_reproduce_fig2_stpwe('paper')
```

`paper` 使用 20 阶正负空间谐波和一阶正负时间谐波，即 $41\times3=123$ 个时空 Fourier 状态、246 阶 Maxwell 广义本征矩阵。

原始 `time.m` 中的着色量实际是 $m=0$ Fourier 扇区占比，并不严格等于论文所述的“投影到所有静态本征模后的能量权重”。本代码把它明确标记为 `m=0 sector participation`，避免把近似量误称为精确静态带权重。能带本征值、复带隙和场重构不依赖这一着色近似。

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

## 5. 时间传输矩阵与时间晶体能带

Ramaccia 模型假定每个时间层在空间上均匀，并在整个空间同时切换。时间界面保持 $D$ 与 $B$ 连续。代码采用列向量，因此按实际时间顺序写成

$$
\mathrm{TM}=
\mathrm{MM}_{M+1}\mathrm{DM}_{M}\cdots
\mathrm{MM}_{2}\mathrm{DM}_{1}\mathrm{MM}_{1}.
$$

`demo04_temporal_multilayer_tmm` 给出四类频率响应。TMM 部分沿用该论文的 $e^{+i\omega t}$ 符号；ST-PWE/FDTD 使用 $e^{-i\omega t}$。这只会改变相位符号，不改变绘制的幅度、带隙位置或增长率绝对值。

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

`demo03_ptc_bands_pwe_vs_tmm` 会把这一 2×2 精确单周期结果与时间 Fourier PWE 结果画在同一张图上，是修改代码后最重要的交叉验证之一。

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

## 7. 拓扑不变量与名称边界

### 7.1 ST-PWE 的空间双正交 Zak 相位

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
demo07_zak_phase_stpwe
```

这里的闭合路径是空间 Brillouin 区的 $k$，因此它是空间 Zak 相。Lustig 等人 PTC 文献中的“时间 Zak 相”以准频率为时间倒空间坐标，需要固定准频率求 $k$ 并重构一周期时间 Bloch 模；当前例程不能互换使用。`demo10` 给出时间畴壁态的独立本征态匹配，但不把该匹配量称为 Zak 相。

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

`demo08_fhs_chern_thouless_pump` 把调制相位作为第二个周期坐标，对 Rice-Mele 泵得到整数 Chern 数 $C=1$。这是通用算法的独立基准；只有在把目标 Maxwell/电路模型的本征态实际送入算法并验证能隙后，才能把整数解释为目标介质的拓扑不变量。

## 8. 科研计算的收敛顺序

建议每次修改模型后依次检查：

1. `test_smoke` 是否通过。
2. 材料 Fourier 系数是否满足实函数的共轭关系
   $p_{-m,-n}=p_{m,n}^*$。
3. 增大空间截断 `Nspace`，确认主要带和复带隙收敛。
4. 增大时间截断 `Mtime`。时间方波比正弦调制需要更多谐波。
5. 用 `demo03` 一类的独立 TMM 结果校验时间 PWE。
6. FDTD 依次减半 $\Delta x$ 与 $\Delta t$，检查带边、增长率和散射幅度。
7. 拓扑量同时检查能隙、Wilson 链接、网格和 Fourier 截断。

只看一张“像论文”的彩图不足以证明数值正确；复谱、独立方法交叉验证和收敛测试应一起保留。

## 9. 适用范围与尚未包含的物理

当前代码针对一维、各向同性、局域、非色散的时空介质。以下问题需要扩展：

- 二维/三维矢量 Maxwell ST-PWE；
- 各向异性、双各向异性或磁电耦合；
- Lorentz/Drude 色散和因果时变色散；
- 严格 CPML 与开放系统准正常模；
- Floquet 异常拓扑相的完整 Rudner winding；
- 简并复合带的非 Abelian 双正交 Wilson 环；
- 非线性调制与泵浦耗尽。

库中的传输线、时变超表面、合成光子网格、共形移动边界和量子紧束缚时空晶体还需要独立状态方程。它们不应通过替换一个 `epsFun` 就被称为已复现。逐篇判断见 `PAPER_CODE_MAP_CN.md`。

尤其要注意：时间 TMM 不能直接替代空间也周期变化的 ST-PWE；反之，有限 Fourier 截断的 PWE 也不应被当作突变时间层的“精确解”。

## 10. 核心直接复现论文

1. J. Park and B. Min, “Spatiotemporal plane wave expansion method for arbitrary space-time periodic photonic media,” *Optics Letters* **46**, 484–487 (2021). DOI: 10.1364/OL.411622.
2. D. Ramaccia, A. Alù, A. Toscano, and F. Bilotti, “Temporal multilayer structures for designing higher-order transfer functions using time-varying metamaterials,” *Applied Physics Letters* **118**, 101901 (2021). DOI: 10.1063/5.0042567.

其余 18 份库文件的“覆盖、部分覆盖、独立分支或重复版本”关系见 `PAPER_CODE_MAP_CN.md`。
