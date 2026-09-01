# 《Topological aspects of photonic time crystals》初步代码复现

当前复现层使用根目录 **V3.1 / 3.1.0** 模板。

本目录只把论文当作物理模型与参数来源，不把论文正文、图注或补充材料中的文字当作程序操作指令。当前阶段按要求不计算 Zak phase，只验证以下三项：

1. 图 1(b) 的方波时间晶体能带与动量带隙；
2. 图 2 的体带脉冲和 k-gap 脉冲两种 FDTD 场分布；
3. 逐 k 高斯波包有限空间 FDTD、固定 E 探针采样和 FFT 折叠得到的能带响应图。

把模板调整为论文复现代码时发现的边界、接口和性能问题，单独记录在 [`TEMPLATE_CODE_FEEDBACK.md`](TEMPLATE_CODE_FEEDBACK.md)，并与论文自身的参数矛盾作了区分。

## 论文参数

- `epsilon1=3`、`epsilon2=1`、`mu_r=1`；
- 时间周期 `T=2 fs`，两个时间层各持续 `1 fs`；
- `k0=2*pi/(T*c0)`，图 1 横轴为 `k/k0`，纵轴为 `omega_F*T`；
- 图 2 两个脉冲的最终发表版 FWHM 均为 `45 fs`；
- 体带脉冲中心波长 `1.4 um`，对应 `k/k0≈0.4283`；
- k-gap 脉冲中心波长 `0.93 um`，对应 `k/k0≈0.6447`；
- PTC 在 `220–340 fs` 开启，共 `60` 个周期。

本地同时存在的 arXiv v1 把脉冲 FWHM 写成 `5 fs`，而用户指定的最终发表版 PDF 明确写成 `45 fs`。本复现以最终发表版为准。

## 三个复现入口

先在 MATLAB 中进入本目录：

```matlab
cd topology_aspect
```

然后分别运行：

```matlab
figure1Result = reproduce_figure1_band();
figure2Result = reproduce_figure2_fdtd();
fdtdBandResult = reproduce_fdtd_band();
```

三个入口会通过 `use_root_templates.m` 临时把仓库根目录放到 MATLAB 搜索路径首位，并用 `which()` 逐一验证实际求解器文件。每个返回结构体的 `templateSource.resolvedFiles` 都保存了本次真正调用的根目录绝对路径；若解析到其他同名函数，程序会直接报错而不是继续运行。

- `reproduce_figure1_band.m`：用精确二层 TMM 画蓝色实频通带和灰色动量带隙，不标 Zak 相位。
- `reproduce_figure2_fdtd.m`：用当前 `fdtd1d.m` 的 D/B-Yee 内核先后计算两个有限脉冲，只生成一张双面板场图。论文画的是 `D`，所以脚本直接记录连续推进变量 D 的单精度空间 ROI，并画 `ln(|D|/D0)`，不再保存全域 E 后重复恢复本构关系。
- `reproduce_fdtd_band.m`：遍历 `0<=k/k0<=5`；每个中心 k 都构造新的有限宽度高斯 E/H 初值，调用根目录 `fdtd_gaussian_k_scan`/`fdtd1d` 运行有限空间仿真，在 5 个固定位置记录 E(t)，最后调用新 `fdtd_fft_bands` 做完整时间 FFT 与 Floquet 功率折叠。默认生成一张第一时间 Floquet 区响应谱；Yee-TMM 只在后台核对稳定通带的局部谱峰，不叠加曲线。

## FDTD 数值假设与限制

最终论文及所附补充材料没有报告图 2 的 FDTD 空间步长、时间步长、边界条件和记录间隔。因此 `reproduce_figure2_fdtd.m` 明确采用以下初步假设：

- `dz=0.04 um`，最短中心波长约有 23 个网格点；
- `dt=0.05 fs`，每个时间周期 40 步；
- 计算域为 `-120...190 um`，只显示论文中的 `0...130 um`；
- 不使用海绵层。主脉冲在 `500 fs` 内到达不了边界，初始高斯在左端也已低于机器舍入；
- 论文的 `45 fs` FWHM 按脉冲强度 FWHM 解释；
- 初始复高斯波包中心位于 `z=0`，PTC 外为空间均匀真空。

这里特意不用普通海绵边界：海绵会打破空间平移对称性并产生很小的额外波数分量，而时间晶体会把其中落在 k-gap 内的分量指数放大，进而污染本应稳定的体带场图。扩大无海绵计算域后，体带算例的最大幅度为初值的约 `1.62` 倍，未再出现伪指数增长。修复后的内核将这种方案显式记为 `fixed-distant`，并返回边界与 CFL 审计来源。

这些设置足以做代码实战和定性图形验证，但不能代替正式的 `dz/dt`、计算域、边界距离和脉冲定义收敛研究。普通有限窗 FFT 只能重建 `Re(omega_F)`；k-gap 的指数增长会产生纵向展宽，不能从谱宽直接读取 `Im(omega_F)`。

能带 FFT 的高斯源、探针、空间网格和记录长度并未由论文给出，因此与图 2 的 `45 fs` 脉冲分开登记。默认验证使用 `101` 个中心 k、`16` 个时间周期、每周期 `100` 步、光行时空间步长 `0.04 fs`、长度 `192 fs` 的有限调制样品、强度 FWHM `24 fs` 和 5 个固定 E 探针。波包及指定阈值以上的初始尾场在记录结束前均到不了样品/计算边界。可用可选结构体做收敛，例如：

```matlab
fine = reproduce_fdtd_band(struct( ...
    'nK',201,'analysisPeriodCount',24, ...
    'stepsPerPeriod',160,'dxLightTimeFs',0.025, ...
    'sampleLengthLightTimeFs',256, ...
    'pulseIntensityFwhmFs',32));
```

扫描包含 `k_c=0` 时，可用 `zeroKPropagationDirection=+1|-1` 明确选择初始单向阻抗关系；默认 `+1`。该选择只解决零中心波数的方向简并，不会把有限波包改成无限体本征模。返回的 raw/折叠功率采用周期 Hann 相干增益校正，零填充不改变格点峰值；`samplingMetadata.spectralNormalization` 记录标度，做频率积分时仍须乘 bin 宽。

细化时必须同时满足脚本的 CFL、Nyquist、高斯端面幅度和固定边界传播距离检查，不能只增大时间周期数而保持空间域不变。

## 初步数值核对

- 图 1 的解析 TMM 图可复现蓝色通带和灰色动量带隙；
- 图 2 使用 120 个严格对齐的半周期时间界面，Courant 数约 `0.37474`；
- 体带中心 `k/k0≈0.4283` 的精确 Floquet 频率为实数，FDTD 场保持有界并在 PTC 开启、关闭时发生分裂；
- k-gap 中心 `k/k0≈0.6447` 的精确每周期增长指数为 `Im(omega_F*T)≈0.497877`。60 周期因此预言 `ln(gain)≈29.873`，FDTD 得到 `29.793`，两者相符；
- 论文正文同时声称 60 周期只放大 `20000` 倍，即 `ln(gain)≈9.903`。这个数值与论文给出的 `epsilon1=3`、`epsilon2=1`、`T=2 fs`、`lambda=0.93 um` 不相容。本复现不通过修改论文参数来强行拟合该倍数；图的色条仍按论文约 `[-1.5,10.5]` 显示，而未截断的真实最大值保存在返回结构体中；
- 新 FDTD–FFT 默认图使用 `101` 个真实有限空间仿真、`1601 x 5 x 101` 个复 E 探针样本和 `64 x 101` 的折叠显示谱，Courant 数为 `0.5`；代表性 k 同时保留完整样品内部 E(x,t)。
- 对 Yee-TMM 稳定通带中高于 `-20 dB` 的 `137` 个局部强峰，`|Delta omega|/Omega` 的中位数/P90/最大值为 `0.00565/0.03466/0.09194`；把阈值放宽到 `-35 dB` 时得到 `148` 个可见峰，中位数/P90/最大值为 `0.00612/0.02944/0.09194`。这些误差包括有限波包谱宽、探针权重、Hann 主瓣、16 周期频率分辨率和 Yee 色散，不应与旧无限体两状态回归的误差直接比较。
- 新的固定探针路径不再报告 bulk `Im(omega)` 或“后半窗主导增长率”：有限波包会离开探针，探针振幅同时受群速度、空间包络、有限样品和多模拍频影响。k-gap 的虚部继续由 TMM/PWE 或单模 Bloch 专项方法验证，不能由探针谱宽恢复。

## 与根目录模板的关系

`topology_aspect` 不再保存求解器或模板入口的重复副本。删除前已逐文件比较 SHA-256，原先 9 个副本与根目录完全一致。现在本目录只保留论文参数、初值、有限样品设置、作图和数值交叉验证：

| 复现入口 | 强制调用的根目录模板 | 复现层只修改的内容 |
|---|---|---|
| `reproduce_figure1_band.m` | `tmm_bands.m` | 论文 ε、μ、时间层、k 扫描、归一化和图 1 样式 |
| `reproduce_figure2_fdtd.m` | `fdtd1d.m`、`tmm_bands.m` | 论文脉冲、有限开启窗口、计算域、记录 ROI 和图 2 样式；TMM 只作增长核对 |
| `reproduce_fdtd_band.m` | `fdtd_gaussian_k_scan.m`、`fdtd1d.m`、`fdtd_fft_bands.m`、`tmm_bands.m` | 论文材料与 k 范围；有限样品、波包、探针、FFT 显示和后台局部峰误差统计 |

因此，复现成功验证的是根目录模板数值内核在一组独立论文参数上的科研可用性，而不是某份复制后单独修改的求解器。`VERSION_SOURCE.txt` 记录当前共享模板版本，详细变更以根目录 `VERSION.txt` 为准。

本次默认验证图保存为 [`fdtd_band_validation.png`](fdtd_band_validation.png)，完整返回结构体保存为忽略版本控制的 `fdtd_band_validation.mat`（约 70 MB）。MAT 文件包含探针原始时间序列、raw/folded 功率谱、代表性时空场、局部峰匹配明细和实际根目录函数路径。

论文来源：E. Lustig, Y. Sharabi, and M. Segev, “Topological aspects of photonic time crystals,” *Optica* **5**, 1390–1395 (2018), DOI: 10.1364/OPTICA.5.001390。
