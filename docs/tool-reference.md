# 函数工具书 (Function Reference)

> 本工具书是**函数级参考**，逐一记录本仓库全部 MATLAB 文件的函数签名、输入参数、输出参数与功能说明。每类方法的完整用法与物理约定见 [learn/ 教程](../learn/README.md) 的对应章节。
> 返回 [README.md](../README.md)

---

## 说明 / About

- **覆盖范围**：本仓库 113 个 `.m` 文件，共 **154 个文档条目**（含嵌套子函数 / 局部辅助函数）。文档基于逐文件阅读源代码生成，函数签名与参数以实际代码为准，不含臆测。
- **阅读方式**：每个函数是一个 `### 函数名` 三级标题，可直接用编辑器（VSCode / MATLAB）的「大纲 / Outline」面板按函数名跳转；每个目录为一章（`##` 二级标题）。
- **符号约定**：时间因子 $e^{i(kx-\omega t)}$、时间界面上 $D$ 与 $B$ 连续、$\operatorname{Im}\omega>0$ 表示时间增长等约定见 [learn/02 章](../learn/02-physics-background.md)（物理背景）与 [learn/05 章 §2](../learn/05-temporal-tmm.md)（界面连续量），各条目的「备注」中也会标注局部约定。
- **依赖关系**：跨目录的调用关系在每个条目的「备注」中说明；入口脚本的运行顺序见 `run_all_demos.m` 与各子项目的 `run_all_reproductions.m`。

## 目录 / Contents

| 章 | 目录 | 文件数 | 条目数 | 说明 |
| ---: | --- | ---: | ---: | --- |
| 1 | `core/` | 18 | 18 | ST-PWE 时空平面波展开核心（Fourier 域求解） |
| 2 | `tmm/` | 12 | 12 | 时间传输矩阵法（精确求解时间多层结构） |
| 3 | `fdtd/` 与 `tests/` | 6 | 9 | 一维 D/B-Yee 时域有限差分与自检 |
| 4 | `topology/` 与根目录 | 8 | 9 | 拓扑不变量（Zak 相位、Chern 数）与顶层脚本 |
| 5 | `demos/` | 12 | 14 | 12 个论文复现演示脚本（v2） |
| 6 | `examples/` | 10 | 14 | 10 个示例脚本（v1） |
| 7 | `reproduction/Observation of temporal reflection/`（`otr_*`） | 18 | 18 | 时间反射实验复现核心函数 |
| 8 | `reproduction/Observation of temporal reflection/`（`fig*`） | 13 | 44 | 时间反射实验复现脚本 |
| 9 | `reproduction/Topology of photonic crystals/` | 11 | 11 | 光子时间晶体拓扑复现脚本 |
| 10 | `Simple_FDTD_NathanZechar/` | 5 | 5 | 第三方 FDTD 参考实现（外部代码） |

> 注：`Simple_FDTD_NathanZechar/` 为第三方开源代码（见其 `license.txt`），此处仅作参考收录，非本工具箱原创。

---
## core/ — ST-PWE 核心（时空平面波展开）

### `stm_fig2_eps_coeff`

**文件：** `core/stm_fig2_eps_coeff.m`

**类型：** 函数

**函数签名：**

```matlab
function coeff = stm_fig2_eps_coeff(m, n, p)
```

**简介：** 解析计算 Park & Min 图 2 时空晶体中相对介电常数 $\varepsilon(x,t)$ 的时空傅里叶系数 $\varepsilon_{m,n}$。

**功能：**

按照约定 $\varepsilon(x,t)=\sum_{m,n}\varepsilon_{m,n}\,e^{i n g x - i m \Omega t}$，返回指定谐波阶 $(m,n)$ 的系数。积分在空间上被拆成静态区 $[0,x_{\rm ModStart})$ 与调制区 $[x_{\rm ModStart},x_{\rm ModEnd})$ 两段，分别用 `stm_interval_fourier` 求解析积分。由于调制是纯正弦的 $\varepsilon_{\rm c}(1+\mathrm{modDepth}\sin\Omega t)$，时间上只有 $m=0,\pm1$ 三个非零项：$m=0$ 给出直流分量 $\varepsilon_1\,\mathrm{staticPart}+\varepsilon_c\,\mathrm{modPart}$，$m=\pm1$ 给出基频分量 $\mp i\,(\varepsilon_c\,\mathrm{modDepth}/2)\,\mathrm{modPart}$，其余阶返回 0。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `m` | double | 时间谐波阶数（整数） |
| `n` | double | 空间谐波阶数（整数） |
| `p` | struct（可选） | 参数结构体；省略或为空时调用 `stm_fig2_parameters()`（默认） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `coeff` | double | 系数 $\varepsilon_{m,n}$，$m=\pm1$ 时为纯虚数 |

**备注：** 依赖 `stm_interval_fourier` 与 `stm_fig2_parameters`。时间因子约定为 $e^{-i m\Omega t}$（与空间因子 $e^{+i n g x}$ 符号相反）。

### `stm_fig2_epsilon`

**文件：** `core/stm_fig2_epsilon.m`

**类型：** 函数

**函数签名：**

```matlab
function epsr = stm_fig2_epsilon(x, t, p)
```

**简介：** 返回 Park & Min 图 2 Floquet 晶体在任意时空点 $(x,t)$ 处的相对介电常数。

**功能：**

输入 `x`、`t` 可为标量或相容数组，函数通过 `x + zeros(size(t))` 与 `t + zeros(size(x))` 把它们广播到同一网格。先取 `xCell = mod(x, Lambda)` 将空间坐标折回一个原胞，再判断是否落在调制区 $[x_{\rm ModStart}, x_{\rm ModEnd})$。静态区取值 $\varepsilon_1$，调制区取 $\varepsilon_c(1+\mathrm{modDepth}\sin(\Omega t))$，其余位置保持 $\varepsilon_1$。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `x` | double | 空间坐标（标量或相容数组） |
| `t` | double | 时间坐标（标量或相容数组） |
| `p` | struct（可选） | 参数结构体；省略或为空时调用 `stm_fig2_parameters()`（默认） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `epsr` | double | 各 $(x,t)$ 处的相对介电常数，尺寸与 `x`/`t` 一致 |

**备注：** 依赖 `stm_fig2_parameters`。调制区为右开区间（`< xModEnd`）。

### `stm_fig2_parameters`

**文件：** `core/stm_fig2_parameters.m`

**类型：** 函数

**函数签名：**

```matlab
function p = stm_fig2_parameters()
```

**简介：** 返回 Park & Min 图 2 复现所用的归一化时空晶体参数。

**功能：**

构造并返回参数结构体 `p`。归一化取 $\Lambda=c_0=1$，从而 $\bar k = k\Lambda/(2\pi)$、$\bar f = \omega\Lambda/(2\pi c_0)$，调制周期为 $T=1/\bar\Omega$。原胞布局为：$[0,\,3\Lambda/4)$ 为静态区 $\varepsilon_1=2$，$[3\Lambda/4,\,\Lambda)$ 为调制区 $\varepsilon_c=6$、调制深度 $\mathrm{modDepth}=0.6$。归一化调制频率 $\bar\Omega=0.20$，对应物理角频率 $\Omega=2\pi c_0/\Lambda\cdot\bar\Omega$ 与周期 $T=2\pi/\Omega$。

**输入参数：**

无

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `p` | struct | 字段：`Lambda`、`c0`、`g`、`eps1`、`epsc`、`modDepth`、`xModStart`、`xModEnd`、`OmegaBar`、`Omega`、`T` |

**备注：** 与 `stm_preset_modulated_slab` 返回相同的数值参数，仅命名/用途不同。

### `stm_fourier_modulated_slab`

**文件：** `core/stm_fourier_modulated_slab.m`

**类型：** 函数

**函数签名：**

```matlab
function coeff = stm_fourier_modulated_slab(m, n, p)
```

**简介：** 解析计算“部分调制平板”原胞介电常数的时空傅里叶系数 $\varepsilon_{m,n}$。

**功能：**

采用约定 $\varepsilon(x,t)=\sum_{m,n}\varepsilon_{m,n}\,e^{i n g x - i m \Omega t}$，返回谐波阶 $(m,n)$ 的系数。空间积分拆成静态区 $[0,x_{\rm ModStart})$ 与调制区 $[x_{\rm ModStart},x_{\rm ModEnd})$ 两段，借助 `stm_interval_fourier` 求解析积分。因为调制是纯正弦的，时间上只有 $m=0,\pm1$ 非零：$m=0$ 为 $\varepsilon_1\,\mathrm{staticPart}+\varepsilon_c\,\mathrm{modPart}$，$m=1$ 为 $i\,(\varepsilon_c\,\mathrm{modDepth}/2)\,\mathrm{modPart}$，$m=-1$ 为 $-i\,(\varepsilon_c\,\mathrm{modDepth}/2)\,\mathrm{modPart}$，其余为 0。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `m` | double | 时间谐波阶数（整数） |
| `n` | double | 空间谐波阶数（整数） |
| `p` | struct（可选） | 参数结构体；省略或为空时调用 `stm_preset_modulated_slab()`（默认） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `coeff` | double | 系数 $\varepsilon_{m,n}$，$m=\pm1$ 时为纯虚数 |

**备注：** 依赖 `stm_interval_fourier` 与 `stm_preset_modulated_slab`。时间因子约定为 $e^{-i m\Omega t}$。

### `stm_interval_fourier`

**文件：** `core/stm_interval_fourier.m`

**类型：** 函数

**函数签名：**

```matlab
function value = stm_interval_fourier(n, xa, xb, Lambda)
```

**简介：** 计算单个空间子区间上的归一化傅里叶积分。

**功能：**

计算 $value = \frac{1}{\Lambda}\int_{x_a}^{x_b} e^{-i n g x}\,dx$，其中 $g=2\pi/\Lambda$。当 $n=0$ 时直接返回 $(x_b-x_a)/\Lambda$；当 $n\neq 0$ 时返回解析结果 $\frac{e^{-i n g x_b}-e^{-i n g x_a}}{-i n g \Lambda}$。该积分是介电常数空间傅里叶系数的基本组成单元，供 `stm_fig2_eps_coeff` 与 `stm_fourier_modulated_slab` 在分段常数介质上调用。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `n` | double | 空间谐波阶数（整数） |
| `xa` | double | 子区间下界 |
| `xb` | double | 子区间上界 |
| `Lambda` | double | 空间周期（归一化通常为 1） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `value` | double | 归一化积分值（$n\neq0$ 时为复数） |

### `stm_permittivity_modulated_slab`

**文件：** `core/stm_permittivity_modulated_slab.m`

**类型：** 函数

**函数签名：**

```matlab
function epsr = stm_permittivity_modulated_slab(x, t, p)
```

**简介：** 返回“部分调制平板”原胞在任意时空点的实空间相对介电常数。

**功能：**

`x`、`t` 可为标量或相容数组，函数把它们广播到同一网格，并用 `mod(x, Lambda)` 将空间坐标折回单个原胞。原胞分段为：$[0,x_{\rm ModStart})$ 静态区取 $\varepsilon_1$，$[x_{\rm ModStart},\Lambda)$ 调制区取 $\varepsilon_c(1+\mathrm{modDepth}\sin(\Omega t))$。返回与输入同尺寸的 `epsr` 数组。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `x` | double | 空间坐标（标量或相容数组） |
| `t` | double | 时间坐标（标量或相容数组） |
| `p` | struct（可选） | 参数结构体；省略或为空时调用 `stm_preset_modulated_slab()`（默认） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `epsr` | double | 各 $(x,t)$ 处的相对介电常数，尺寸与 `x`/`t` 一致 |

**备注：** 依赖 `stm_preset_modulated_slab`。调制区为右开区间（`< xModEnd`）。

### `stm_preset_modulated_slab`

**文件：** `core/stm_preset_modulated_slab.m`

**类型：** 函数

**函数签名：**

```matlab
function p = stm_preset_modulated_slab()
```

**简介：** 返回“部分调制平板”单元胞的参数结构体。

**功能：**

返回一个描述 1D 时空晶体的参数结构体 `p`：单元胞中一部分为静态介质、其余部分随时间正弦调制，是 ST-PWE 研究的常用参考构型。归一化取 $\Lambda=c_0=1$，归一化波矢 $\bar k=k\Lambda/(2\pi)$、归一化频率 $\bar f=\omega\Lambda/(2\pi c_0)$。单元胞布局为 $[0,x_{\rm ModStart})$ 静态 $\varepsilon_1=2$、$[x_{\rm ModStart},\Lambda)$ 调制 $\varepsilon_c(1+\mathrm{modDepth}\sin\Omega t)$，其中 $\varepsilon_c=6$、$\mathrm{modDepth}=0.6$、$x_{\rm ModStart}=3\Lambda/4$，归一化调制频率 $\bar\Omega=0.20$。

**输入参数：**

无

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `p` | struct | 字段：`Lambda`、`c0`、`g`、`eps1`、`epsc`、`modDepth`、`xModStart`、`xModEnd`、`OmegaBar`、`Omega`、`T` |

**备注：** 数值与 `stm_fig2_parameters` 相同；被 `stm_permittivity_modulated_slab` 与 `stm_fourier_modulated_slab` 作为默认参数使用。

### `stm_redblue`

**文件：** `core/stm_redblue.m`

**类型：** 函数

**函数签名：**

```matlab
function cmap = stm_redblue(n)
```

**简介：** 生成蓝—白—红的发散型（diverging）colormap。

**功能：**

生成一个 $n\times3$ 的 RGB 颜色映射矩阵，从蓝色渐变到白色再渐变到红色，用于可视化有正负号的数据（如场的实部、能带或傅里叶系数）。将 $n$ 分成前一半（`n1 = floor(n/2)`）与后一半（`n2 = n - n1`）：前段由蓝色 `[0,0,1]` 线性过渡到白色，后段由白色过渡到红色 `[1,0,0]`，两段拼成完整 colormap。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `n` | double（可选） | 颜色数；省略时默认 256 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `cmap` | double（$n\times3$） | RGB colormap 矩阵，每行一个颜色 |

### `stm_stft`

**文件：** `core/stm_stft.m`

**类型：** 函数

**函数签名：**

```matlab
function [spectrogramValue, f, tCenter] = stm_stft(signal, dt, windowLength, hop, nFFT)
```

**简介：** 不依赖工具箱的短时傅里叶变换（STFT）。

**功能：**

对输入信号做加窗短时傅里叶变换，返回单边（正频率）复数频谱、频率轴与各帧中心时间。信号先被拉成列向量，默认窗长为 128、帧移 `hop = floor(windowLength/8)`、FFT 点数 `nFFT = 2^nextpow2(2*windowLength)`。窗口为 Hann 窗 $0.5-0.5\cos(2\pi j/(L-1))$。以步长 `hop` 取帧，每帧乘窗后做 `nFFT` 点 FFT，只保留前 `floor(nFFT/2)+1` 个正频率分量。频率轴为 `(0:nPositive-1)/(nFFT*dt)`，时间中心为 `((starts-1)+(windowLength-1)/2)*dt`。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `signal` | double | 输入时域信号（任意方向，内部转列向量） |
| `dt` | double | 采样时间间隔 |
| `windowLength` | double（可选） | 窗长（默认 128） |
| `hop` | double（可选） | 帧移采样点数（默认 `floor(windowLength/8)`） |
| `nFFT` | double（可选） | FFT 点数（默认 `2^nextpow2(2*windowLength)`） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `spectrogramValue` | double complex | STFT 复数谱，尺寸 $(\lfloor nFFT/2\rfloor+1)\times N_{\rm frames}$ |
| `f` | double（列向量） | 频率轴（Hz） |
| `tCenter` | double（行向量） | 各帧中心时间 |

**备注：** 仅返回正频率半谱（`0..Nyquist`），不含负频率部分。

### `stpwe_build_system`

**文件：** `core/stpwe_build_system.m`

**类型：** 函数

**函数签名：**

```matlab
function sys = stpwe_build_system(epsCoeff, muCoeff, Nspace, Mtime, g, Omega)
```

**简介：** 构建 1D 时空平面波展开（ST-PWE）所需的卷积矩阵与系统结构体。

**功能：**

根据介电/磁导率的时空傅里叶系数构造 ST-PWE 离散系统的全部矩阵。谐波阶截断为 $n=-N_{\rm space}:N_{\rm space}$、$m=-M_{\rm time}:M_{\rm time}$，总自由度 $S=(2N_{\rm space}+1)(2M_{\rm time}+1)$。卷积矩阵按 Toeplitz 结构填充：$C_{\rm eps}(row,col)=\varepsilon_{dm,dn}$、$C_{\rm mu}(row,col)=\mu_{dm,dn}$，其中 $dm=m_{row}-m_{col}$、$dn=n_{row}-n_{col}$。同时构造对角矩阵 $G=\mathrm{diag}(n g)$、$W=\mathrm{diag}(m\Omega)$，以及一阶形式的 $B_\omega=[0,\,C_{\rm mu};\,C_{\rm eps},\,0]$。`muCoeff` 缺省时为 $\delta_{m,0}\delta_{n,0}$（即 $\mu=1$ 常数）。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `epsCoeff` | function_handle | 返回 $\varepsilon_{m,n}$ 的函数 `epsCoeff(m,n)` |
| `muCoeff` | function_handle（可选） | 返回 $\mu_{m,n}$ 的函数；默认 `@(m,n) double(m==0 && n==0)` |
| `Nspace` | double | 空间谐波截断阶数（保留 $n=-Nspace:Nspace$） |
| `Mtime` | double | 时间谐波截断阶数（保留 $m=-Mtime:Mtime$） |
| `g` | double | 空间倒格矢 $g=2\pi/\Lambda$ |
| `Omega` | double | 调制角频率 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `sys` | struct | 字段：`Nspace`、`Mtime`、`nList`、`mList`、`S`、`g`、`Omega`、`Ceps`、`Cmu`、`G`、`W`、`I`、`Z`、`Bomega` |

**备注：** 系数约定 $p(x,t)=\sum_{m,n}p_{m,n}e^{i n g x - i m \Omega t}$。`sys` 作为后续 `stpwe_solve_k` / `stpwe_solve_omega` / `stpwe_select_mode` 的输入。

### `stpwe_fold_frequency`

**文件：** `core/stpwe_fold_frequency.m`

**类型：** 函数

**函数签名：**

```matlab
function omegaFolded = stpwe_fold_frequency(omega, Omega)
```

**简介：** 将复数频率的实部折叠到第一布里渊区 $[-\Omega/2,\,\Omega/2)$。

**功能：**

对每个复频率 $\omega$，仅折叠其实部：$\mathrm{Re}(\omega_{\rm folded}) = \bmod\big(\mathrm{Re}(\omega)+\Omega/2,\ \Omega\big)-\Omega/2$，虚部保持不变。用于把 ST-PWE 解出的本征频率折叠到调制频率的约化区间，便于绘制折叠能带图。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `omega` | double | 复数频率（标量或数组） |
| `Omega` | double | 调制角频率 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `omegaFolded` | double | 实部折叠后的复数频率 |

### `stpwe_lookup_coefficient`

**文件：** `core/stpwe_lookup_coefficient.m`

**类型：** 函数

**函数签名：**

```matlab
function value = stpwe_lookup_coefficient(table, m, n)
```

**简介：** 从采样系数表中读取谐波阶 $(m,n)$ 的系数，越界时返回 0。

**功能：**

在 `table.mOrders` 与 `table.nOrders` 中查找给定的时间阶 `m` 与空间阶 `n`；若两者都找到则返回 `table.values(im,in)`，否则返回 0。用于以统一的“函数式”接口访问数值采样的傅里叶系数，使采样表可被 `stpwe_build_system` 等期望 `coeff(m,n)` 函数句柄的接口消费。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `table` | struct | 系数表，含字段 `mOrders`、`nOrders`、`values`（见 `stpwe_sample_fourier_coefficients`） |
| `m` | double | 时间谐波阶数 |
| `n` | double | 空间谐波阶数 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `value` | double | 对应系数；阶数不在表中时返回 0 |

### `stpwe_reconstruct_field`

**文件：** `core/stpwe_reconstruct_field.m`

**类型：** 函数

**函数签名：**

```matlab
function field = stpwe_reconstruct_field(mode, x, t, includeGrowth)
```

**简介：** 由单个 Floquet 本征模重建实空间-时间上的实电场。

**功能：**

按完整场表达式 $E(x,t)=\mathrm{Re}\big[e^{i(kx-\omega t)}\sum_{mn} E_{mn} e^{i(n g x - m\Omega t)}\big]$ 重建电场。先 `meshgrid` 生成时空网格，累加周期部分 `periodicPart`，再乘载波 $e^{i(k x-\omega t)}$ 取实部。默认 $\omega$ 只用实部（`includeGrowth=false`）；若 `includeGrowth=true` 则保留复数 $\omega$，从而把增长率纳入场幅度。最后按峰值 `max(|field|)` 归一化。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `mode` | struct | 本征模结构体，含 `Ecoef`、`nList`、`mList`、`g`、`Omega`、`k`、`omega`（见 `stpwe_select_mode`） |
| `x` | double | 空间坐标向量 |
| `t` | double | 时间坐标向量 |
| `includeGrowth` | logical（可选） | 是否保留 $\omega$ 虚部（增长率）；默认 false |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `field` | double | 重建的实电场，尺寸 `numel(t) × numel(x)`，按峰值归一化 |

**备注：** 约定载波 $e^{+i k x}$、时间因子 $e^{-i\omega t}$、空间谐波 $e^{+i n g x}$。

### `stpwe_sample_fourier_coefficients`

**文件：** `core/stpwe_sample_fourier_coefficients.m`

**类型：** 函数

**函数签名：**

```matlab
function table = stpwe_sample_fourier_coefficients(materialFun, mOrders, nOrders, Lambda, T, Nx, Nt)
```

**简介：** 对任意 $p(x,t)$ 用数值采样求其时空傅里叶系数，返回系数表。

**功能：**

对任意材料参数函数 `materialFun(X,Time)`（需返回 `Nt×Nx` 数组）在端点不重复的均匀网格上采样，用离散平均近似傅里叶系数 $p_{mn}=\langle p(x,t)\,e^{-i n g x + i m\Omega t}\rangle$。网格为 $x=(0:N_x-1)\Lambda/N_x$、$t=(0:N_t-1)T/N_t$，$g=2\pi/\Lambda$、$\Omega=2\pi/T$。对每个 $(m,n)$ 构造核 $e^{-i n g X + i m\Omega t}$，取 `mean(values.*kernel,'all')` 作为系数。结果连同阶数、周期、采样点数一起存入结构体 `table`。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `materialFun` | function_handle | 返回材料参数的函数 `materialFun(X,Time)`，需返回 `Nt×Nx` 数组 |
| `mOrders` | double | 需求的时间谐波阶数向量 |
| `nOrders` | double | 需求的空间谐波阶数向量 |
| `Lambda` | double | 空间周期 |
| `T` | double | 时间周期 |
| `Nx` | double（可选） | 空间采样点数（默认 `max(128, 8*numel(nOrders))`） |
| `Nt` | double（可选） | 时间采样点数（默认 `max(128, 8*numel(mOrders))`） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `table` | struct | 字段：`mOrders`、`nOrders`、`values`、`Lambda`、`T`、`Nx`、`Nt` |

**备注：** 系数约定 $p(x,t)=\sum p_{mn}e^{i n g x - i m\Omega t}$（核的虚部符号与此一致）。采样为端点不重复（`0..Nx-1`）网格。

### `stpwe_select_mode`

**文件：** `core/stpwe_select_mode.m`

**类型：** 函数

**函数签名：**

```matlab
function mode = stpwe_select_mode(sys, kTarget, omegaTarget, omegaWindow)
```

**简介：** 在固定波矢 $k$ 的本征解中挑选最接近目标频率的可见 Floquet 模。

**功能：**

先用 `stpwe_solve_omega` 在 `kTarget` 处求解本征值，筛出有限且在 `omegaWindow` 实频窗内的解。然后按加权得分挑选：$\mathrm{score}=\frac{|\mathrm{Re}(\omega)-\omega_{\rm target}|}{scale}+0.15\frac{|\mathrm{Im}(\omega)|}{scale}+0.02(1-m_0Weight)$，其中 $scale=\max(|\omega_{\rm target}|,\Omega)$，得分最小者入选。选出后对本征向量做相位归一（令第一个 $S$ 分量中模最大的元素相位为 0），并拆分电场/磁场系数，打包成 `mode` 结构体。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `sys` | struct | `stpwe_build_system` 产出的系统结构体 |
| `kTarget` | double | 目标波矢 |
| `omegaTarget` | double | 目标（复数）频率 |
| `omegaWindow` | double（可选） | 实频筛选窗 `[min,max]`；默认 `[-inf, inf]` |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `mode` | struct | 字段：`k`、`omega`、`Ecoef`、`Hcoef`、`nList`、`mList`、`g`、`Omega`、`m0Weight` |

**备注：** 依赖 `stpwe_solve_omega`。偏好 $m_0$ 权重高（基模成分大）且虚部小的解；若窗内无有效解则报错。

### `stpwe_solve_k`

**文件：** `core/stpwe_solve_k.m`

**类型：** 函数

**函数签名：**

```matlab
function sol = stpwe_solve_k(sys, omega)
```

**简介：** 求解固定频率 $\omega$ 下的 ST-PWE 本征问题，得到波矢本征值 $k$（Park-Min 式(7)）。

**功能：**

给定 $\omega$，构造 $OW=\omega I+W$ 与矩阵 $A=[-G,\ OW\,C_{\rm mu};\ OW\,C_{\rm eps},\ -G]$，用 `eig` 求其本征值即波矢 $k$。对每个本征向量用左右特征向量内积归一：若 $L_q^\dagger R_q$ 有限且模大于 $10^{-12}$，则 $R_q\leftarrow R_q/(L_q^\dagger R_q)$。返回本征值 `k`、右/左特征向量矩阵 `R`/`L` 与矩阵 `A`。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `sys` | struct | `stpwe_build_system` 产出的系统结构体 |
| `omega` | double | 给定的（复数）频率 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `sol` | struct | 字段：`omega`、`k`、`R`、`L`、`A` |

### `stpwe_solve_omega`

**文件：** `core/stpwe_solve_omega.m`

**类型：** 函数

**函数签名：**

```matlab
function sol = stpwe_solve_omega(sys, k)
```

**简介：** 求解固定波矢 $k$ 下的 ST-PWE 广义本征问题，得到频率本征值 $\omega$（Park-Min 式(6)）。

**功能：**

给定 $k$，构造 $K=k I+G$、$A=[K,\ -W C_{\rm mu};\ -W C_{\rm eps},\ K]$ 与 $B=B_\omega$，求解广义本征问题 $A R=\omega B R$。左右特征向量满足 $L^\dagger A=\omega L^\dagger B$，若非自正交则归一化到 $L_i^\dagger B R_i=1$，否则各自按范数归一。最后计算基模权重 $m_0Weight=\frac{\sum_{m=0}(|E|^2+|H|^2)}{\sum(|E|^2+|H|^2)}$，用于衡量 $m=0$ 分量的占比。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `sys` | struct | `stpwe_build_system` 产出的系统结构体 |
| `k` | double | 给定的波矢 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `sol` | struct | 字段：`k`、`omega`、`R`、`L`、`A`、`B`、`m0Weight` |

**备注：** `R` 前 $S$ 行为电场分量、后 $S$ 行为磁场分量。`m0Weight` 被 `stpwe_select_mode` 用作模式筛选依据。

### `stpwe_static_bands`

**文件：** `core/stpwe_static_bands.m`

**类型：** 函数

**函数签名：**

```matlab
function fBands = stpwe_static_bands(kValues, epsCoeff0, Nspace, g, c0, nBands)
```

**简介：** 在与 ST-PWE 相同的一阶形式下计算静态 1D 光子晶体的能带。

**功能：**

按空间谐波截断 $n=-N_{\rm space}:N_{\rm space}$ 构造倒格矢对角阵 $G=\mathrm{diag}(n g)$，并用 `epsCoeff0`（仅空间，$p(x)=\sum_n p_n e^{i n g x}$）填充卷积矩阵 $C_{\rm eps}$，再组装 $B=[0,\,I;\,C_{\rm eps},\,0]$。对每个波矢 $k$，求解广义本征问题 $[K,0;0,K]\,R=\omega B R$（$K=k I+G$），保留有限、虚部 $<10^{-8}$、实部 $\geq-10^{-9}$ 的实频本征值并升序排列。在 $k=0$ 处 $K$ 奇异会额外产生一个近零的虚假本征值，代码检测到多个近零值时删除全部再补回一个 0。最终归一化频率 $\bar f=\omega/(g c_0)$（因 $\Lambda=2\pi/g$）。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `kValues` | double | 波矢采样点向量 |
| `epsCoeff0` | function_handle | 返回空间系数 $\varepsilon_n$ 的函数 `epsCoeff0(n)` |
| `Nspace` | double | 空间谐波截断阶数 |
| `g` | double | 空间倒格矢 $g=2\pi/\Lambda$ |
| `c0` | double | 归一化光速 |
| `nBands` | double | 需输出的能带条数 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `fBands` | double（`nBands × numel(kValues)`） | 归一化频率 $\bar f=\omega\Lambda/(2\pi c_0)$ 的能带；不足 `nBands` 处为 NaN |

**备注：** 依赖 $g$、$c_0$ 的归一化关系 $\bar f=\omega/(g c_0)$。含 $\Gamma$ 点奇异本征值的去重处理，避免能带索引在 $k=0$ 处跳变。
## tmm/ — 时间传输矩阵法

### `temporal_binary_eps_coeff`

**文件：** `tmm/temporal_binary_eps_coeff.m`

**类型：** 函数

**函数签名：**

```matlab
function coeff = temporal_binary_eps_coeff(m, epsA, epsB, dutyA)
```

**简介：** 计算二值（分段常数）时间晶体的第 m 阶 Fourier 介电系数。

**功能：**

对周期为 T 的二值介电函数 ε(t)（在一个周期内 0 ≤ t < dutyA·T 取 epsA，其余时间取 epsB）按 ε(t) = Σ_m ε_m exp(−i·m·Ω·t) 展开，返回第 m 阶系数 ε_m。零阶项（m = 0）等于时间平均 dutyA·epsA + (1−dutyA)·epsB；非零阶项由解析积分 (epsA−epsB)·(e^{i·2π·m·dutyA}−1)/(i·2π·m) 给出。结果为复数（非零阶）。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `m` | double（整数） | Fourier 阶数，可取任意整数（含负值）。 |
| `epsA` | double | 相位 A 的介电常数。 |
| `epsB` | double | 相位 B 的介电常数。 |
| `dutyA` | double | 相位 A 的占空比，介于 0 与 1 之间。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `coeff` | double | 第 m 阶 Fourier 系数 ε_m。 |

**备注：** 采用时间因子 e^{−i·m·Ω·t} 的约定；本函数只返回单个系数，用于供其它模块构造介电序列。

### `temporal_crystal_bands`

**文件：** `tmm/temporal_crystal_bands.m`

**类型：** 函数

**函数签名：**

```matlab
function bands = temporal_crystal_bands(kValues, epsSequence, ...
    muSequence, durations)
```

**简介：** 由时间单周期单值矩阵求 Floquet 准频率（能带）。

**功能：**

对每个守恒波数 k，调用 `temporal_crystal_monodromy` 得到一周期 D/B 演化矩阵 U，取其本征值 λ，按 ω_F = i·log(λ)/T 反解出 Floquet 准频率。由于主值对数取主分支，Re(ω_F) 落在第一个时间布里渊区 [−π/T, π/T) 内。两个本征值按 (Re ω, Im ω) 排序，与 λ、半迹 trace(U)/2 一起存入结构体返回。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `kValues` | double（向量） | 待扫描的守恒波数 k 序列。 |
| `epsSequence` | double（向量） | 一个周期内各时间层的介电常数序列。 |
| `muSequence` | double（向量/标量） | 各时间层的磁导率序列（标量时自动广播）。 |
| `durations` | double（向量） | 各时间层的持续时间（厚度），周期 T = sum(durations)。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `bands` | struct | 含字段 `k`、`omegaF`（2×Nk 准频率）、`lambda`（2×Nk 本征值）、`halfTrace`、`T`（周期）、`Omega`（Ω = 2π/T）。 |

**备注：** 依赖 `temporal_crystal_monodromy`。归一化单位 eps0=mu0=c0=1。

### `temporal_crystal_monodromy`

**文件：** `tmm/temporal_crystal_monodromy.m`

**类型：** 函数

**函数签名：**

```matlab
function U = temporal_crystal_monodromy(k, epsSequence, ...
    muSequence, durations)
```

**简介：** 计算守恒波数 k 下一个时间周期的 D/B 单值矩阵。

**功能：**

对平面波 exp(i·k·x)，D、B 服从 d/dt [D;B] = −i·k·[0 1/mu; 1/eps 0]·[D;B]。按时间顺序逐层用矩阵指数 expm(generator·duration) 左乘累积，得到整个周期的单值矩阵 U。由于 D、B 在任意时间界面处连续，此状态表示下无需单独的界面矩阵。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `k` | double（标量） | 守恒波数。 |
| `epsSequence` | double（向量） | 一周期内各时间层介电常数序列。 |
| `muSequence` | double（向量/标量） | 各时间层磁导率序列，标量时复制到与 epsSequence 等长。 |
| `durations` | double（向量） | 各时间层持续时间，须与 epsSequence 等长。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `U` | double（2×2 复数） | 一个周期内的 D/B 单值矩阵。 |

**备注：** 采用归一化单位 eps0=mu0=c0=1。epsSequence、muSequence、durations 长度不一致时报错。

### `temporal_db_to_directional`

**文件：** `tmm/temporal_db_to_directional.m`

**类型：** 函数

**函数签名：**

```matlab
function amplitudes = temporal_db_to_directional(state, epsr, mur)
```

**简介：** 将 [D;B] 状态转换为 [E_plus; E_minus] 方向振幅。

**功能：**

这是 `temporal_directional_to_db` 的逆变换。给定 n = sqrt(epsr·mur)，由 D、B 反解出前向/后向电场振幅：E_plus = 0.5·(D/epsr + B/n)，E_minus = 0.5·(D/epsr − B/n)。输入 `state` 为 2×N，逐列处理。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `state` | double（2×N 复数） | D、B 状态列向量（第 1 行 D，第 2 行 B）。 |
| `epsr` | double（标量） | 相对介电常数，须非零。 |
| `mur` | double（标量） | 相对磁导率，须非零。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `amplitudes` | double（2×N 复数） | 方向振幅，第 1 行 E_plus，第 2 行 E_minus。 |

**备注：** 约定 exp(i·k·x−i·ω·t)、k>0，归一化单位。epsr 或 mur 非标量或为零时报错。

### `temporal_delay_matrix`

**文件：** `tmm/temporal_delay_matrix.m`

**类型：** 函数

**函数签名：**

```matlab
function DM = temporal_delay_matrix( ...
    omegaInitial, nInitial, nSlab, duration)
```

**简介：** 计算单个时间层内的相位累积（延迟）矩阵。

**功能：**

由守恒波数得层内频率 ω_slab = (nInitial/nSlab)·ω_initial，相位 phase = ω_slab·duration，返回对角矩阵 diag([e^{i·phase}, e^{−i·phase}])，作用于 [E_plus; E_minus] 列向量，表示前向与后向波在层内传播的相位积累。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `omegaInitial` | double | 初始（背景）介质中的角频率。 |
| `nInitial` | double | 初始介质折射率。 |
| `nSlab` | double | 当前时间层的折射率。 |
| `duration` | double | 该时间层的持续时间。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `DM` | double（2×2 复数对角矩阵） | 层内相位累积矩阵。 |

**备注：** 采用 Ramaccia 等的 exp(+i·ω·t) 约定；本函数不单独使用，通常由 `temporal_multilayer_tmm` 调用。

### `temporal_directional_to_db`

**文件：** `tmm/temporal_directional_to_db.m`

**类型：** 函数

**函数签名：**

```matlab
function state = temporal_directional_to_db(amplitudes, epsr, mur)
```

**简介：** 将 [E_plus; E_minus] 方向振幅转换为 [D;B] 状态。

**功能：**

在约定 exp(i·k·x−i·ω·t)、k>0 及归一化单位下，两个时间频率分支满足 D = epsr·(E_plus+E_minus)、B = n·(E_plus−E_minus)，其中 n = sqrt(epsr·mur)。输入 `amplitudes` 为 2×N，逐列转换为 2×N 的 [D;B] 状态。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `amplitudes` | double（2×N 复数） | 方向振幅列向量（第 1 行 E_plus，第 2 行 E_minus）。 |
| `epsr` | double（标量） | 相对介电常数，须非零。 |
| `mur` | double（标量） | 相对磁导率，须非零。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `state` | double（2×N 复数） | 第 1 行 D，第 2 行 B。 |

**备注：** epsr 或 mur 非标量或为零时报错；是 `temporal_db_to_directional` 的逆。

### `temporal_domain_wall_mode`

**文件：** `tmm/temporal_domain_wall_mode.m`

**类型：** 函数

**函数签名：**

```matlab
function mode = temporal_domain_wall_mode(kValues, ...
    epsLeft, muLeft, durationsLeft, ...
    epsRight, muRight, durationsRight, nLeft, nRight)
```

**简介：** 匹配时间域壁上增长/衰减的 Floquet 态，求界面束缚模式。

**功能：**

左侧时间晶体占据域壁前 NLEFT 个周期、右侧晶体占据域壁后 NRIGHT 个周期。在共同的动量带隙内，局域在域壁的态在左单值矩阵下朝壁增长、在右单值矩阵下离开壁衰减。程序扫描守恒 k，用左、右单值矩阵的本征向量（左取 |λ| 最大、右取 |λ| 最小）构造两个一维 Floquet 本征空间，最小化其行列式失配 det([vGrowth, vDecay])。随后用 fminbnd 在最佳 k 邻域细化，重构界面两侧各周期的 D/B 态并归一化。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `kValues` | double（向量） | 扫描的守恒波数序列。 |
| `epsLeft` | double（向量） | 左晶体的介电常数序列。 |
| `muLeft` | double（向量/标量） | 左晶体的磁导率序列。 |
| `durationsLeft` | double（向量） | 左晶体各层持续时间。 |
| `epsRight` | double（向量） | 右晶体的介电常数序列。 |
| `muRight` | double（向量/标量） | 右晶体的磁导率序列。 |
| `durationsRight` | double（向量） | 右晶体各层持续时间。 |
| `nLeft` | double（正整数） | 域壁左侧的周期数。 |
| `nRight` | double（正整数） | 域壁右侧的周期数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `mode` | struct | 含 `kGrid`、`mismatch`、`growthMagnitude`、`decayMagnitude`、`leftHalfTrace`、`rightHalfTrace`、最佳 `k`、`bestMismatch`、`leftMultiplier`、`rightMultiplier`、`leftEigenvector`、`rightEigenvector`、`cellIndex`、`states`、`stateNorm`、`interfaceId`、`leftMonodromy`、`rightMonodromy`、`T` 等字段。 |

**备注：** 依赖 `temporal_crystal_monodromy`；内含子函数 `wall_mismatch`（对给定 k 计算失配）。左右晶体周期须相同，nLeft/nRight 须为正整数，否则报错。本函数只给出指定传输矩阵的体-界面模式，不自行计算或证明量子化时间 Zak 相位。

### `temporal_finite_crystal_response`

**文件：** `tmm/temporal_finite_crystal_response.m`

**类型：** 函数

**函数签名：**

```matlab
function response = temporal_finite_crystal_response( ...
    kValues, epsSequence, muSequence, durations, nPeriods, ...
    epsBackground, muBackground, inputDirectional)
```

**简介：** 计算有限时间晶体的精确响应（透射/反射方向振幅）。

**功能：**

空间均匀的时间晶体在 NPERIODS 个周期内被调制，嵌入在两段相同的静态时间背景之间，输入与输出均为背景介质中的方向振幅 [E_plus; E_minus]。调制开启与关闭处 D、B 连续。对每个 k，先把输入方向振幅经 `temporal_directional_to_db` 转成 D/B 初态，再用单值矩阵 U 的 nPeriods 次幂演化，最后经 `temporal_db_to_directional` 转回方向振幅，同时记录 U、本征值（乘子）等。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `kValues` | double（向量） | 守恒波数序列。 |
| `epsSequence` | double（向量） | 一周期内介电常数序列。 |
| `muSequence` | double（向量/标量） | 一周期内磁导率序列。 |
| `durations` | double（向量） | 各层持续时间。 |
| `nPeriods` | double（非负整数） | 调制持续的周期数。 |
| `epsBackground` | double | 背景介电常数。 |
| `muBackground` | double | 背景磁导率。 |
| `inputDirectional` | double（2×1 列向量，可选） | 输入方向振幅，默认 [1;0]（纯前向）。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `response` | struct | 含 `k`、`inputDirectional`、`outputDirectional`、`forward`、`backward`、`relativePhase`、`finalState`、`periodMatrix`、`multipliers`、`nPeriods`、`T` 等字段。 |

**备注：** 依赖 `temporal_crystal_monodromy`、`temporal_directional_to_db`、`temporal_db_to_directional`。inputDirectional 须为 2×1 列向量，nPeriods 须为非负整数。

### `temporal_interface_matrix`

**文件：** `tmm/temporal_interface_matrix.m`

**类型：** 函数

**函数签名：**

```matlab
function [MM, tau, rho] = temporal_interface_matrix( ...
    epsBefore, muBefore, epsAfter, muAfter)
```

**简介：** 计算突变时间界面处的 Morgenthaler 矩阵（D/B 连续模型）。

**功能：**

对列向量 [E_plus; E_minus]，由 D、B 连续性给出 after = MM·before，其中 MM = [tau rho; rho tau]。系数 tau、rho 由界面两侧参数决定。这是 Ramaccia 等 APL 118, 101901 (2021) 的式 (5)，实际上是调用 `temporal_interface_matrix_jump` 并指定 'DB' 模型。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `epsBefore` | double | 界面之前的介电常数。 |
| `muBefore` | double | 界面之前的磁导率。 |
| `epsAfter` | double | 界面之后的介电常数。 |
| `muAfter` | double | 界面之后的磁导率。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `MM` | double（2×2） | 时间界面传输矩阵。 |
| `tau` | double | 对角系数（平均项）。 |
| `rho` | double | 非对角系数（反射项）。 |

**备注：** D/B 连续是本构切换模型，并非适用于所有实验时间界面的普适规则；需要保留 E、注入电荷或其它跳变规律时改用 `temporal_interface_matrix_jump`。

### `temporal_interface_matrix_jump`

**文件：** `tmm/temporal_interface_matrix_jump.m`

**类型：** 函数

**函数签名：**

```matlab
function [MM, tau, rho, details] = temporal_interface_matrix_jump( ...
    epsBefore, muBefore, epsAfter, muAfter, model)
```

**简介：** 带显式跳变因子的时间界面传输矩阵。

**功能：**

作用在 [E_plus; E_minus] 上，after = MM·before。跳变因子定义为 D_after = jumpD·D_before、B_after = jumpB·B_before。`model` 可选 'DB'（D、B 连续，标准体时间边界）、'EB'（E、B 连续）、'DH'（D、H 连续）、'EH'（E、H 连续），或 struct('jumpD',qD,'jumpB',qB) 自定义微观模型。由跳变因子计算电、磁比值，进而得 tau = 0.5·(ratioElectric+ratioMagnetic)、rho = 0.5·(ratioElectric−ratioMagnetic)，构造 MM = [tau rho; rho tau]。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `epsBefore` | double | 界面之前介电常数，须非零。 |
| `muBefore` | double | 界面之前磁导率，须非零。 |
| `epsAfter` | double | 界面之后介电常数，须非零。 |
| `muAfter` | double | 界面之后磁导率，须非零。 |
| `model` | char/string/struct（可选） | 跳变模型名或自定义 jumpD/jumpB 结构，默认 'DB'。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `MM` | double（2×2） | 时间界面传输矩阵。 |
| `tau` | double | 对角系数。 |
| `rho` | double | 非对角系数。 |
| `details` | struct | 含 `model`、`jumpD`、`jumpB`、`electricRatio`、`magneticRatio`。 |

**备注：** 命名 E/B/D/H 选项是理想化跳变规律；含源、断连电荷、空间局域切换或时间色散的开关电路一般需要更大状态空间，不能用两个标量跳变因子表示。

### `temporal_multilayer_tmm`

**文件：** `tmm/temporal_multilayer_tmm.m`

**类型：** 函数

**函数签名：**

```matlab
function [TM, details] = temporal_multilayer_tmm(omegaInitial, ...
    epsInitial, muInitial, epsSlabs, muSlabs, durations, ...
    epsFinal, muFinal)
```

**简介：** 计算任意时间层叠的总传输矩阵。

**功能：**

场为列向量，按时间顺序（chronological）映射 TM = MM_final·DM_M·…·MM_2·DM_1·MM_1。对每个时间层依次：用 `temporal_interface_matrix` 得到进入层的匹配矩阵，用 `temporal_delay_matrix` 得到层内延迟矩阵，再左乘下一匹配矩阵累积。若层数 M=0，则退化为单次界面矩阵。中间各匹配矩阵与延迟矩阵均存入 details 供检查。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `omegaInitial` | double | 初始介质中的角频率。 |
| `epsInitial` | double | 初始介质介电常数。 |
| `muInitial` | double | 初始介质磁导率。 |
| `epsSlabs` | double（向量） | 各时间层介电常数序列。 |
| `muSlabs` | double（向量/标量） | 各时间层磁导率序列（标量广播）。 |
| `durations` | double（向量） | 各时间层持续时间，须与 epsSlabs 等长。 |
| `epsFinal` | double | 末层之后介电常数。 |
| `muFinal` | double | 末层之后磁导率。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `TM` | double（2×2 复数） | 整个时间层叠的传输矩阵。 |
| `details` | struct | 含 `tau`、`rho`、`nInitial`、`nSlabs`、`omegaSlabs`、`matchingMatrices`、`delayMatrices`。 |

**备注：** 依赖 `temporal_interface_matrix`、`temporal_delay_matrix`。明确给出左乘顺序以避免左右作用约定的歧义。

### `temporal_tmm_spectrum`

**文件：** `tmm/temporal_tmm_spectrum.m`

**类型：** 函数

**函数签名：**

```matlab
function spectrum = temporal_tmm_spectrum(omegaInitial, ...
    epsInitial, muInitial, epsSlabs, muSlabs, durations, ...
    epsFinal, muFinal)
```

**简介：** 对频率序列计算前向/后向电场透反射系数谱。

**功能：**

对每个输入角频率 ω_initial，调用 `temporal_multilayer_tmm` 得到该频率下的传输矩阵 TM，令输入 [1;0]（纯前向单位振幅），则 output = TM·[1;0]，第一分量 forward、第二分量 backward 即为出射的前向/后向电场系数，同时记录 det(TM)。所有频率结果存入结构体返回。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `omegaInitial` | double（向量） | 输入角频率序列。 |
| `epsInitial` | double | 初始介质介电常数。 |
| `muInitial` | double | 初始介质磁导率。 |
| `epsSlabs` | double（向量） | 各时间层介电常数序列。 |
| `muSlabs` | double（向量/标量） | 各时间层磁导率序列。 |
| `durations` | double（向量） | 各时间层持续时间。 |
| `epsFinal` | double | 末层之后介电常数。 |
| `muFinal` | double | 末层之后磁导率。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `spectrum` | struct | 含 `omegaInitial`、`forward`、`backward`、`detTM` 字段。 |

**备注：** 依赖 `temporal_multilayer_tmm`。
## fdtd/ 与 tests/ — 时域有限差分与测试

### `fdtd1d_db`

**文件：** `fdtd/fdtd1d_db.m`

**类型：** 函数

**函数签名：**

```matlab
function out = fdtd1d_db(cfg)
```

**简介：** 一维 D/B-Yee 时域有限差分（FDTD）求解器，适用于时空变化介质。

**功能：**

本函数在一维 Yee 网格上实现 leapfrog（蛙跳）时间推进。与传统 E/H FDTD 不同，它直接更新电位移场 D 和磁感应场 B，而非电场 E 和磁场 H：因为在时间界面（介电常数/磁导率突变的时刻）处 D 和 B 自动满足连续性条件 D(t⁺)=D(t⁻)、B(t⁺)=B(t⁻)，而每个时间步通过本构关系 E=D/ε(x,t)、H=B/μ(x,t) 恢复 E/H，可自然处理任意时变材料参数。

空间上采用 Yee 交错网格：E（整数格点）与 H（半整数格点，偏移 dx/2）错开半个步长；时间上同样交错半个时间步 dt/2，初始条件 E0 定义在 t=0，Hhalf0 定义在 t=-dt/2。时间推进由法拉第定律 B^{n+1/2}=B^{n-1/2}-(Δt/Δx)·curl(E^n) 与安培定律 D^{n+1}=D^n-(Δt/Δx)·curl(H^{n+1/2}) 交替完成，其中在显式时间界面处 E 用 0.5·[D/ε(t⁻)+D/ε(t⁺)] 的界面平均场计算旋度以消除一阶时间误差。

支持两种边界：sponge（三次方衰减剖面 damp=exp(-strength·s³) 的简化吸收层，默认）与 periodic（通过 circshift 环绕的周期边界）。返回结构体 `out`，包含记录时刻的 E/H/D 历史矩阵（H 已插值到 E 网格）、探针记录、瞬时电磁能量时间序列（0.5·dx·ΣRe[E*·D + H*·B]）以及 CFL 审计得到的最大波速与 Courant 数。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `cfg` | struct | 仿真配置结构体，必填/可选字段如下。 |
| `cfg.x` | double | 均匀空间网格（行或列向量，至少 5 点），定义 E/D 场格点。 |
| `cfg.dt` | double | 时间步长，须满足 CFL 条件 Δt < Δx/v_max。 |
| `cfg.nSteps` | integer | 总推进步数，总仿真时间 = nSteps·dt。 |
| `cfg.epsFun` | function_handle | 相对介电常数函数句柄 `epsFun(x,t)`，输入空间坐标数组 x 与标量时间 t，返回长度与 x 相同的 ε_r 数组。 |
| `cfg.muFun` | function_handle（可选） | 相对磁导率函数句柄 `muFun(x,t)`（默认 `@(x,t) ones(size(x))`）。 |
| `cfg.E0` | double（可选） | 初始电场（t=0，E 网格，长度 Nx，默认全零）。 |
| `cfg.Hhalf0` | double（可选） | 初始磁场（t=-dt/2，H 网格，长度同 xH，默认全零）。 |
| `cfg.boundary` | char（可选） | 边界条件 `'sponge'`（默认）或 `'periodic'`。 |
| `cfg.spongeCells` | integer（可选） | 海绵吸收层宽度（格点数，默认 `min(80, floor(Nx/8))`）。 |
| `cfg.spongeStrength` | double（可选） | 吸收强度（默认 0.12）。 |
| `cfg.recordEvery` | integer（可选） | 场记录间隔（每 N 步记录一次，默认 1）。 |
| `cfg.storeFields` | logical（可选） | 是否存储完整 E/H 历史（默认 true；false 可大幅降低内存占用）。 |
| `cfg.storeD` | logical（可选） | 是否独立存储 E 网格上的 D 历史（默认 false，与 storeFields 相互独立）。 |
| `cfg.temporalInterfaces` | double（可选） | 介电常数突变时刻的严格递增实数向量（默认空；须与 E 整数时间节点 n·dt 对齐且小于 nSteps·dt）。 |
| `cfg.temporalInterfaceTolerance` | double（可选） | temporalInterfaces 的绝对对齐容差（默认由 dt 与总时长的浮点精度确定，须小于 dt/4）。 |
| `cfg.spectralFilterMask` | double/logical（可选） | 周期网格 FFT 掩码（长度 Nx，默认空）；若给出则在每个时间界面更新后同时过滤 D 与 B。 |
| `cfg.probeIndices` | integer（可选） | 探针点索引（E 网格，1..Nx；即使 storeFields=false 也记录）。 |
| `cfg.stabilityTimes` | double（可选） | CFL 审计的采样时刻数组（默认 `linspace(0,nSteps*dt,9)`）。 |
| `cfg.maxWaveSpeed` | double（可选） | 手动指定的最大波速（跳过采样审计）。 |
| `cfg.sourceD` | function_handle（可选） | D 场增量源函数句柄 `sourceD(x,t,n)`，在每步更新 D 后调用（软源注入）。 |
| `cfg.precision` | char（可选） | 数值精度 `'double'`（默认）或 `'single'`。 |
| `cfg.progressBar` | logical（可选） | 是否显示文本进度条（默认 false）。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `out` | struct | 仿真结果结构体，字段如下。 |
| `out.x` | double | E 空间网格。 |
| `out.xH` | double | H 空间网格（相对 E 偏移 dx/2）。 |
| `out.t` | double | 记录时刻数组（列向量）。 |
| `out.E` | complex | 电场历史矩阵 (nRecords × Nx)。 |
| `out.H` | complex | 磁场历史矩阵 (nRecords × Nx)，已插值到 E 网格。 |
| `out.D` | complex | 电位移历史矩阵 (nRecords × Nx)；storeD=false 时为空。 |
| `out.probeIndices` | integer | 探针索引。 |
| `out.probeX` | double | 探针空间位置。 |
| `out.probeE` | complex | 探针处电场记录。 |
| `out.probeH` | complex | 探针处磁场记录。 |
| `out.energy` | double | 瞬时电磁能量时间序列。 |
| `out.finalD` / `out.finalB` / `out.finalE` / `out.finalH` | complex | 最终时刻的 D/B/E/H 场。 |
| `out.dx` / `out.dt` | double | 空间步长 / 时间步长。 |
| `out.boundary` | char | 使用的边界条件类型。 |
| `out.sampledMaxWaveSpeed` | double | CFL 审计得到的最大波速（可能为 NaN）。 |
| `out.sampledCourant` | double | CFL 审计得到的 Courant 数。 |
| `out.precision` | char | 数值精度。 |

**备注：** 依赖局部辅助函数 `yee_h_to_e_grid`（见下）。物理约定：D/B 在时间界面处连续；能量密度对非色散介质为 (1/2)(E·D + H·B)，时变介质中调制度会与场交换能量、总能量一般不守恒。稳定性要求 Courant 数 S = v_max·Δt/Δx < 1（v_max = max[1/√(ε_r·μ_r)]）；对色散/负折射率/有源介质，采样审计返回 NaN 并跳过检查，需通过 `cfg.maxWaveSpeed` 手动指定。sponge 不是严格 CPML，而是简化的乘法衰减。`cfg.spectralFilterMask` 仅对 `boundary='periodic'` 且存在 temporalInterfaces 时有效。

### `yee_h_to_e_grid`

**文件：** `fdtd/fdtd1d_db.m`（局部辅助函数）

**类型：** 函数

**函数签名：**

```matlab
function [HOnE, BOnE] = yee_h_to_e_grid(H, B, boundary, Nx)
```

**简介：** 将 H/B 从 H 网格（半整数格点）插值到 E 网格（整数格点）。

**功能：**

Yee 网格中 E 与 H 在空间上错开 dx/2，本函数把定义在 H 网格上的磁场 H 与磁感应场 B 线性插值到 E 网格的同一空间位置，便于统一绘图和分析。周期边界下所有点都用相邻 H 值的算术平均 `0.5*(H + circshift(H,1))`；吸收边界下内部点用中心平均 `0.5*(H(i-1)+H(i))`，两个端点用最近邻单侧值（零阶外推，只影响 sponge 吸收层内的边界点）。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `H` | double | H 网格上的磁场数组。 |
| `B` | double | H 网格上的磁感应场数组。 |
| `boundary` | char | 边界条件 `'periodic'` 或 `'sponge'`。 |
| `Nx` | integer | E 网格点数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `HOnE` | double | 插值到 E 网格的磁场（长度 Nx）。 |
| `BOnE` | double | 插值到 E 网格的磁感应场（长度 Nx）。 |

——

### `animate_fdtd1d`

**文件：** `fdtd/animate_fdtd1d.m`

**类型：** 函数

**函数签名：**

```matlab
function animate_fdtd1d(out, options)
```

**简介：** 生成 FDTD 1D D/B-Yee 仿真结果的波包时空演化动画。

**功能：**

以 `fdtd1d_db` 的输出结构体 `out` 为输入，播放电场 E(x,t)（或磁场 H、或二者）随时间的演化动画。支持两种视图模式：`line`（一维线图，实时更新 YData）与 `waterfall`（时空瀑布图，用 imagesc 逐行累积显示 E(x,t) 二维分布，并叠加红色进度线 yline 标注当前时刻）。可选的 `fields` 字段控制显示 E、H 或 both；`showEnergy=true` 时叠加一个能量曲线子图并用红色圆点标记当前时刻。

支持将动画保存为 GIF（`rgb2ind` + `imwrite` 逐帧追加）或视频（VideoWriter，.mp4 用 MPEG-4、.avi 用 Motion JPEG AVI）。交互模式下通过 `pause` 按 `speed` 倍率控制播放速度；`skipFrames` 实现跳帧（每 N 帧取一帧）。全部帧基于 `real(fieldData)` 的实部绘制。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `out` | struct | `fdtd1d_db` 的输出结构体（需含 `x`、`t`、`E`/`H`，可选 `energy`）。 |
| `options` | struct（可选） | 可选动画配置结构体，字段如下。 |
| `options.saveAs` | char（可选） | 保存文件名，支持 .gif、.mp4、.avi（默认 `''`，仅显示不保存）。 |
| `options.fields` | char（可选） | 显示的场：`'E'`、`'H'` 或 `'both'`（默认 `'E'`）。 |
| `options.viewMode` | char（可选） | `'line'`（线图）或 `'waterfall'`（时空瀑布图，默认 `'line'`）。 |
| `options.fps` | double（可选） | 保存视频/动画的帧率（默认 30）。 |
| `options.speed` | double（可选） | 播放速度倍率（默认 1）。 |
| `options.xlim` | double（可选） | x 轴范围 `[xmin,xmax]`（默认自动）。 |
| `options.ylim` | double（可选） | 场幅值范围 `[ymin,ymax]`（默认按最大幅值的 ±1.1 倍自动对称）。 |
| `options.showTime` | logical（可选） | 是否显示时间标签（默认 true）。 |
| `options.showEnergy` | logical（可选） | 是否叠加能量曲线子图（默认 false）。 |
| `options.cmap` | char（可选） | waterfall 模式的 colormap（默认 `'parula'`）。 |
| `options.title` | char（可选） | 动画标题（默认自动生成）。 |
| `options.figSize` | double（可选） | 图形窗口大小 `[w,h]` 像素（默认 `[800,500]`）。 |
| `options.skipFrames` | integer（可选） | 跳帧播放，每 N 帧取一帧（默认 1）。 |

**输出参数：** 无（仅绘图/写文件）。若指定 `saveAs` 且扩展名为 .gif，则保存 GIF 动画；若为 .mp4/.avi 则保存视频文件；否则仅在图形窗口中播放。

**备注：** 依赖 `fdtd1d_db` 的输出结构，以及局部辅助函数 `getOpt`（见下）。瀑布图模式只对播放帧逐行累积，无 NaN 间隙，并固定 y 轴范围防止动画过程中跳变。保存模式下图形窗口设为不可见以加速渲染。

### `getOpt`

**文件：** `fdtd/animate_fdtd1d.m`（局部辅助函数）

**类型：** 函数

**函数签名：**

```matlab
function val = getOpt(opts, fieldName, defaultVal)
```

**简介：** 从 options 结构体提取字段，字段不存在或为空时返回默认值。

**功能：**

检查结构体 `opts` 中是否存在字段 `fieldName` 且其值非空，是则返回该字段值，否则返回默认值 `defaultVal`。用于 `animate_fdtd1d` 中解析可选的 `options` 配置字段。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `opts` | struct | 待提取字段的结构体。 |
| `fieldName` | char | 字段名。 |
| `defaultVal` | 任意 | 字段不存在或为空时的默认返回值。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `val` | 任意 | 字段值或默认值。 |

——

### `test_fdtd_optimizations`

**文件：** `fdtd/test_fdtd_optimizations.m`

**类型：** 脚本

**函数签名：** 脚本 (无函数签名)

**简介：** 测试 FDTD 求解器的单精度、进度条与动画功能，并诊断波包传播。

**功能：**

采用归一化单位（c0=1，lambda0=1），构造正向传播的高斯波包（`profile = exp(-((xq-x0)/sigma)^2)·exp(1i·k0·(xq-x0))`，H 与 E 阻抗匹配）作为初始条件，通过 `fdtd1d_db` 分别以双精度与单精度运行真空介质仿真（CFL=0.80，sponge 边界，`progressBar=true`）。测试 1 诊断波包峰值位置并对比实际传播距离与期望值 c0·t；测试 2 对比双精度与单精度结果的相对误差并报告内存节省。

随后调用 `animate_fdtd1d` 分别生成 line 模式与 waterfall 模式的 GIF 动画（`test_wave_line.gif`、`test_wave_waterfall.gif`），最后绘制三个时刻（t=0、t=T/2、t=T）的 E(x) 快照图并保存为 `test_diagnostic_snapshot.png`。

**输入参数：** 无

**输出参数：** 无（仅绘图/写文件）。保存 `test_wave_line.gif`、`test_wave_waterfall.gif` 与 `test_diagnostic_snapshot.png`。

**备注：** 依赖 `fdtd1d_db` 与 `animate_fdtd1d`。正向行波 H 初始条件按 `Hhalf0 = profile(xH + c0·dt/2)` 构造（真空 n=1，故 H=1·E）。

——

### `test_temporal_interface_animation`

**文件：** `fdtd/test_temporal_interface_animation.m`

**类型：** 脚本

**函数签名：** 脚本 (无函数签名)

**简介：** 时间界面（time boundary）处波包演化的动画与诊断演示。

**功能：**

在 t=tSwitch=15 时刻令介电常数突变（折射率 nBefore=1.5 → nAfter=2.5），通过 `fdtd1d_db` 的 `temporalInterfaces` 参数模拟时间折射与时间反射（正向透射波与反向相位共轭反射波）。介电常数由辅助函数 `temporal_eps` 给出分段常数，H 初始条件按折射率 nBefore 介质中的正向波构造（Hhalf0 = nBefore·profile(xH + vBefore·dt/2)）。

仿真后对界面后的场做方向分解 E_± = (E ± H/nAfter)/2，用 L2 范数计算数值透射系数 |tau| 与反射系数 |rho|，并与解析 Morgenthaler 系数 `tauExact = n1(n1+n2)/(2·n2²)`、`rhoExact = n1(n1-n2)/(2·n2²)` 对照。最后调用 `animate_fdtd1d` 生成 line 与 waterfall 两种 GIF（`test_ti_line.gif`、`test_ti_waterfall.gif`），并绘制诊断快照（|E(x,t)| 瀑布图、方向分解、能量随时间变化）保存为 `test_ti_diagnostic.png`。

**输入参数：** 无

**输出参数：** 无（仅绘图/写文件）。保存 `test_ti_line.gif`、`test_ti_waterfall.gif` 与 `test_ti_diagnostic.png`。

**备注：** 依赖 `fdtd1d_db`、`animate_fdtd1d` 与局部辅助函数 `temporal_eps`（见下）。注释指出这些是 E 场系数，D/B 连续会带来相对常见归一化振幅多一个 n1/n2 因子。

### `temporal_eps`

**文件：** `fdtd/test_temporal_interface_animation.m`（局部辅助函数）

**类型：** 函数

**函数签名：**

```matlab
function epsr = temporal_eps(x, t, tSwitch, epsBefore, epsAfter)
```

**简介：** 时间界面两侧的分段常数介电常数。

**功能：**

在 t < tSwitch 时返回 epsBefore 常数数组，在 t ≥ tSwitch 时返回 epsAfter 常数数组，数组长度与输入空间坐标 x 一致。用于构造 `fdtd1d_db` 所需的时变 `epsFun` 函数句柄。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `x` | double | 空间坐标数组（仅用于确定输出尺寸）。 |
| `t` | double | 标量时间。 |
| `tSwitch` | double | 介电常数切换时刻。 |
| `epsBefore` | double | 切换前的介电常数（ε_r = nBefore²）。 |
| `epsAfter` | double | 切换后的介电常数（ε_r = nAfter²）。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `epsr` | double | 与 `x` 等长的相对介电常数数组。 |

——

### `stm_selftest`

**文件：** `tests/stm_selftest.m`

**类型：** 函数

**函数签名：**

```matlab
function stm_selftest()
```

**简介：** STM 工具箱核心引擎、TMM、FDTD 与拓扑例程的快速数值/代数自检。

**功能：**

安装后运行以验证核心引擎、TMM、FDTD 与拓扑例程在容差范围内输出正确结果。首先调用 `stm_init` 初始化，然后依次做五组断言：① 采样与解析的时空 Fourier 系数一致（通过 `stpwe_sample_fourier_coefficients` 与 `stm_fourier_modulated_slab` 对比，误差 < 3e-2）；② 零时长中间介质匹配矩阵等于直接时间界面矩阵（`temporal_interface_matrix` 与 `temporal_multilayer_tmm` 的 Frobenius 范数差 < 1e-12）；③ 无耗时间晶体 monodromy 行列式为 1（`temporal_crystal_monodromy`）；④ Rice-Mele 泵的小网格 FHS Chern 数等于 1（`fhs_chern_number`）；⑤ 短时周期 FDTD 均匀平面波能量漂移有界（< 8%）。全部通过后打印各指标汇总。

**输入参数：** 无

**输出参数：** 无（通过 `assert` 断言，失败抛出错误；通过时打印自检结果）。

**备注：** 依赖 `stm_init` 及 STM 工具箱多个函数（`stm_preset_modulated_slab`、`stpwe_sample_fourier_coefficients`、`stpwe_lookup_coefficient`、`stm_fourier_modulated_slab`、`temporal_interface_matrix`、`temporal_multilayer_tmm`、`temporal_crystal_monodromy`、`fhs_chern_number`、`fdtd1d_db`）。FDTD 段使用周期边界、n=1.5 均匀介质，初始场为单色平面波 `E0=exp(i·k·x)`。

——

### `test_smoke`

**文件：** `tests/test_smoke.m`

**类型：** 函数

**函数签名：**

```matlab
function test_smoke()
```

**简介：** 包级快速数值与代数冒烟测试。

**功能：**

调用 `startup_stm` 后对包的核心例程做一组快速数值/代数检查，内容与 `stm_selftest` 部分重叠但更全面。前两部分（Fourier 系数采样一致性、时间界面矩阵复合性）与 `stm_selftest` 相同；此外还验证：命名跳变律保持其名称所述量（`temporal_interface_matrix_jump` 的 'DB' 保持 D/B、'EB' 保持 E/B，以及双输入相干抵消）；无耗时间晶体 monodromy 行列式为 1；有限晶体响应（`temporal_finite_crystal_response`）与显式矩阵幂一致；反向二进制 PTC 的时间畴壁模式失配最小且在界面处包络峰值（`temporal_domain_wall_mode`）。

随后做更严格的 FDTD 检查：短时周期均匀介质能量漂移 < 8% 且 CFL 审计通过、探针记录正确；独立 D 历史存储（storeFields=false、storeD=true）时 D 尺寸正确且 E/H 为空、初始 D 记录等于 ε·E0；一个二进制时间晶胞的事件感知 FDTD 与精确 D/B monodromy 一致（相对状态误差与增益误差均 < 1e-3）；以及界面时刻谱投影（spectralFilterMask）能移除支撑外的 Fourier 分量且不泄漏（< 1e-12）。全部通过后打印各指标汇总。

**输入参数：** 无

**输出参数：** 无（通过 `assert` 断言，失败抛出错误；通过时打印测试结果）。

**备注：** 依赖 `startup_stm` 及大量 STM 工具箱函数（含 `stm_fig2_parameters`、`stm_fig2_epsilon`、`stm_fig2_eps_coeff`、`temporal_interface_matrix_jump`、`temporal_directional_to_db`、`temporal_finite_crystal_response`、`temporal_domain_wall_mode`、`fhs_chern_number`、`fdtd1d_db` 等）。最后一个 FDTD 检查中 finalB 位于 T-dt/2 半整数时刻，比较前先推进半个步长到整数时刻 T。
## topology/ 与根目录 — 拓扑不变量与顶层脚本

### `fhs_chern_number`

**文件：** `topology/fhs_chern_number.m`

**类型：** 函数

**函数签名：**

```matlab
function [chern, curvature, links] = fhs_chern_number(rightStates, leftStates, metric)
```

**简介：** 在周期性二维动量网格上用 Fukui–Hatsugai–Suzuki（FHS）离散化计算陈数（Chern number）。

**功能：**

在 (k, p) 两个方向的闭合周期网格上，用每个网格点处的占据态构造幺模 link 变量，再由 link 变量的乘积角度得到每个 plaquette 的 Berry 曲率，最后对全部 plaquette 求和除以 2π 得到陈数。每个 link 由重叠矩阵的行列式归一化得到：`det(L'*metric*R)`，取其相位作为 link。对于厄米问题可省略左矢与度规（此时左矢取右矢、度规取单位阵）；对于广义/非厄米问题则需提供双正交左矢 `leftStates` 与重叠度规 `metric`。若某个 link 的模小于 1e-13 则报错，提示该网格点 link 奇异。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `rightStates` | complex double | 右本征态，尺寸 `[Nbasis, Noccupied, Nk, Np]` |
| `leftStates` | complex double（可选） | 双正交左本征态，尺寸须与 `rightStates` 一致（默认取 `rightStates`） |
| `metric` | double（可选） | 重叠度规矩阵，尺寸 `Nbasis x Nbasis`（默认单位阵 `eye(Nbasis)`） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `chern` | double | 陈数，`sum(curvature)/(2π)` |
| `curvature` | double | `Nk x Np` 矩阵，每个 plaquette 的 Berry 曲率（`angle` 取值） |
| `links` | struct | 含字段 `k` 与 `parameter`，均为 `Nk x Np` 的幺模复数 link 变量 |

**备注：** 两个网格方向均按周期闭合（用 `mod(ik,Nk)+1` 回绕），网格不得包含重复的端点。该实现与 `stpwe_bz_sewing_matrix` / `zak_phase_biorthogonal` 独立，自包含无外部依赖。

### `stpwe_bz_sewing_matrix`

**文件：** `topology/stpwe_bz_sewing_matrix.m`

**类型：** 函数

**函数签名：**

```matlab
function sewing = stpwe_bz_sewing_matrix(sys)
```

**简介：** 构造将 k=−g/2 处的态缝接到 k=+g/2 处等价态的 sewing 矩阵。

**功能：**

根据时空平面波展开 `E = exp(i·k·x)·Σ_n u_n·exp(i·n·g·x)`，同一物理态满足 `u_n(k+g) = u_{n+1}(k)`，即在平移一个倒格矢 g 时系数整体右移一格。函数在系统 `sys` 的 `mList`、`nList` 索引表上构建一个 `S x S` 的 shift 矩阵：当 `mList(row) == mList(col)` 且 `nList(col) == nList(row)+1` 时置 1；随后用 `blkdiag(shift, shift)` 复制两份，得到 `2S x 2S` 的 sewing 矩阵。有限傅里叶截断会丢掉最外侧一个系数，因此要求该截断处权重可忽略以保证收敛。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `sys` | struct | 时空平面波系统对象，须含字段 `S`（截断阶数）、`mList`、`nList`（谐波索引表） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `sewing` | double | `2S x 2S` 的 sewing 矩阵（`blkdiag(shift, shift)`） |

**备注：** 依赖 `sys.S`、`sys.mList`、`sys.nList` 字段；物理约定为 `exp(i·k·x)` 平面波展开下的系数平移关系。

### `stpwe_track_band`

**文件：** `topology/stpwe_track_band.m`

**类型：** 函数

**函数签名：**

```matlab
function band = stpwe_track_band(sys, kGrid, omegaSeed, options)
```

**简介：** 通过本征矢连续性追踪一条 Floquet 能带。

**功能：**

沿 `kGrid` 逐点求解系统（调用 `stpwe_solve_omega`），按 `omegaWindow`、`maxImag`、`minM0Weight` 过滤候选模式后，用代价函数挑选与上一点本征态最连续的模式，从而追踪单一能带。首点代价取 `|omega − omegaSeed|/scale + 0.01·(1 − m0Weight)`；其后各点代价取 `1 − min(fidelity,1) + 0.15·frequencyStep + 0.01·(1 − m0Weight)`，其中 fidelity 用左右广义本征矢经 `sys.Bomega` 度规的 forward/backward 重叠之几何平均 `sqrt(|forward·backward|)` 定义，使其在左右矢各自独立缩放时不变。若重叠过低或相邻带隙闭合，提示单带 Zak 相位不可靠，应改用复合带 Wilson loop。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `sys` | struct | 时空平面波系统对象，须含字段 `S`、`Omega`、`Bomega` |
| `kGrid` | double | k 点序列（被展平为行向量） |
| `omegaSeed` | complex double | 初始频率种子，用于首点选模 |
| `options` | struct（可选） | 可选参数，字段见下 |

`options` 可选字段：

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `omegaWindow` | double(1x2) | 频率实部窗口（默认 `[-inf inf]`） |
| `maxImag` | double | 允许的最大 `|Im(omega)|`（默认 `inf`） |
| `minM0Weight` | double | 最小 m=0 参与权重（默认 `0`） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `band` | struct | 追踪结果，含字段 `k`、`omega`、`R`、`L`、`m0Weight`、`selectedId`、`continuity`、`neighborGap` |

`band` 字段说明：

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `k` | double | k 网格（1 x Nk） |
| `omega` | complex double | 追踪到的本征频率（1 x Nk） |
| `R` | complex double | 右本征矢（`2S x Nk`） |
| `L` | complex double | 左本征矢（`2S x Nk`） |
| `m0Weight` | double | 每点 m=0 参与权重（1 x Nk） |
| `selectedId` | double | 每点被选模式的候选索引（1 x Nk） |
| `continuity` | double | 相邻点间的本征态重叠（1 x Nk） |
| `neighborGap` | double | 与最近邻其他模式的频率间距（1 x Nk） |

**备注：** 依赖 `stpwe_solve_omega(sys, k)`。相邻点通过 `exp(−1i·angle(forward))` 做相位规整以保持本征矢连续。

### `zak_phase_biorthogonal`

**文件：** `topology/zak_phase_biorthogonal.m`

**类型：** 函数

**函数签名：**

```matlab
function [zak, info] = zak_phase_biorthogonal(rightStates, leftStates, metric, sewing)
```

**简介：** 计算规范不变的单带 Wilson loop 与 Zak 相位（双正交版本）。

**功能：**

对一列 k 网格点上的左右广义本征矢，构造相邻点间的 link `L_j'*metric*R_{j+1}`，最后一点经 `sewing` 矩阵映射回第一点（对应 k→k+G）。将所有 link 归一化为单位模后连乘得到 Wilson loop，Zak 相位取 `wrap_to_pi(−angle(wilsonLoop))`。该构造对标量规范变换不敏感，适用于非厄米/双正交本征问题。若任一 link 模小于 1e-12 则报错，提示能带可能与其他带接触或傅里叶截断过小。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `rightStates` | complex double | 广义右本征矢，尺寸 `Nbasis x Nk`（每列一个本征态） |
| `leftStates` | complex double | 广义左本征矢，尺寸须与 `rightStates` 一致 |
| `metric` | double | 重叠度规矩阵（`Nbasis x Nbasis`） |
| `sewing` | double（可选） | 将第一点态映射到 k+G 的矩阵（默认单位阵 `eye(Nbasis)`） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `zak` | double | Zak 相位，`wrap_to_pi(−angle(wilsonLoop))` |
| `info` | struct | 含字段 `links`、`unitLinks`、`wilsonLoop`、`minimumLinkMagnitude` |

**备注：** k 网格不得包含重复的右端点。该函数与 `stpwe_bz_sewing_matrix` 配合使用（sewing 由后者构造）。自包含，另含局部辅助函数 `wrap_to_pi`。

### `wrap_to_pi`

**文件：** `topology/zak_phase_biorthogonal.m`

**类型：** 函数（局部子函数）

**函数签名：**

```matlab
function value = wrap_to_pi(value)
```

**简介：** 将角度折返到 (−π, π] 区间。

**功能：**

对输入值做 `mod(value + pi, 2π) − pi`，把任意实数相位折返到主值区间 (−π, π]。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `value` | double | 待折返的角度（弧度） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `value` | double | 折返到 (−π, π] 的角度 |

### `run_all_demos`

**文件：** `run_all_demos.m`

**类型：** 函数

**函数签名：**

```matlab
function run_all_demos(runHeavy)
```

**简介：** 按合理顺序运行软件包内全部演示脚本。

**功能：**

先调用 `startup_stm()` 配置路径，然后依次运行快演示 `demo02_complex_frequency_and_momentum_gaps`、`demo03_ptc_bands_pwe_vs_tmm`、`demo04_temporal_multilayer_tmm`、`demo05_fdtd_temporal_interface`、`demo08_fhs_chern_thouless_pump`、`demo09_finite_ptc_order_and_phase`、`demo10_temporal_domain_wall_mode`、`demo11_coherent_time_interface`。当 `runHeavy` 为真时，额外运行论文级 `demo01_reproduce_fig2_stpwe('paper')`、`demo06_fdtd_spacetime_wavepacket`、`demo07_zak_phase_stpwe`、`demo12_ptc_convergence_audit`；否则打印跳过提示。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `runHeavy` | logical（可选） | 是否运行重计算演示（默认 `false`） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 仅运行各演示脚本并打印进度（各演示自行绘图/写文件） |

**备注：** 依赖 `startup_stm` 与 `demo01`–`demo12` 系列函数；重计算演示较长（FDTD、收敛性审计等）。

### `startup_stm`

**文件：** `startup_stm.m`

**类型：** 函数

**函数签名：**

```matlab
function rootDir = startup_stm()
```

**简介：** 将软件包内所有文件夹加入 MATLAB 搜索路径。

**功能：**

以本文件所在目录为根目录，把根目录及 `core`、`tmm`、`fdtd`、`topology`、`demos`、`tests` 子目录加入 `addpath`；若 `output` 目录不存在则创建之，并打印就绪信息与根目录路径。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无输入 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `rootDir` | char | 软件包根目录的绝对路径 |

**备注：** 会创建 `output/` 目录；与 `stm_init` 的区别在于本函数加入 `demos` 目录而 `stm_init` 加入 `examples` 目录。

### `stm_init`

**文件：** `stm_init.m`

**类型：** 函数

**函数签名：**

```matlab
function rootDir = stm_init()
```

**简介：** Space-Time Media (STM) 工具箱的入口，将各文件夹加入 MATLAB 搜索路径。

**功能：**

以本文件所在目录为根目录，把根目录及 `core`、`tmm`、`fdtd`、`topology`、`examples`、`tests` 子目录加入 `addpath`；若 `output` 目录不存在则创建之，并打印就绪信息与根目录路径。是 STM 工具箱的统一初始化入口。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无输入 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `rootDir` | char | 工具箱根目录的绝对路径 |

**备注：** 会创建 `output/` 目录；与 `startup_stm` 的区别在于本函数加入 `examples` 目录而 `startup_stm` 加入 `demos` 目录。

### `stm_run_examples`

**文件：** `stm_run_examples.m`

**类型：** 函数

**函数签名：**

```matlab
function stm_run_examples(runHeavy)
```

**简介：** 按合理顺序运行 STM 工具箱的全部示例。

**功能：**

先调用 `stm_init()` 配置路径，然后依次运行快示例 `example_complex_gaps`、`example_ptc_pwe_vs_tmm`、`example_tmm_multilayer`、`example_fdtd_interface`、`example_chern_pump`。当 `runHeavy` 为真时，额外运行 `example_floquet_bands_and_fields('paper')`、`example_fdtd_wavepacket`、`example_zak_phase`；否则打印跳过提示。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `runHeavy` | logical（可选） | 是否运行重计算示例（默认 `false`） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 仅运行各示例脚本并打印进度（各示例自行绘图/写文件） |

**备注：** 依赖 `stm_init` 与 `example_*` 系列函数；重计算示例（论文级能带图、FDTD 波包、Zak 相位）较长。
## demos/ — 12 个论文复现演示脚本

### `demo01_reproduce_fig2_stpwe`

**文件：** `demos/demo01_reproduce_fig2_stpwe.m`

**类型：** 函数

**函数签名：**

```matlab
function demo01_reproduce_fig2_stpwe(quality)
```

**简介：** 用时空平面波展开（ST-PWE）复现 Park-Min 论文图 2 的 Floquet 能带结构。

**功能：**

调用 `stpwe_build_system` 按图 2 的时空调制介电常数构造广义本征问题，然后对一系列约化波矢 k 逐点用 `stpwe_solve_omega` 求解复本征频率 ω，并按实部落在 `[0, fMax]`、虚部幅值小于 `imagTolerance`、且 m=0 Floquet 扇区参与度 `m0Weight` 高于 `weightThreshold` 的条件筛选解，得到能带散点。同时用 `stpwe_static_bands` 计算静态参考能带并绘制其 Floquet 平移副本（shift = ±1, ±2），再用 `stpwe_select_mode` 在六个标注点 (c)–(h) 处选取代表模式、用 `stpwe_reconstruct_field` 重构场分布。散点颜色表示 m=0 扇区参与度，作为论文中基于投影的静态带权重的可复现替代。最后绘制空间介电常数曲面、能带图与六个场图并保存 PNG。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `quality` | char（可选） | 计算精度，`'quick'` 或 `'paper'`，默认 `'quick'`。`'paper'` 取 Nspace=20、Nk=181；`'quick'` 取 Nspace=10、Nk=101。其他值报错。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无返回值，仅绘图并保存 `output/demo01_fig2_stpwe_<quality>.png`。 |

**备注：** 依赖 `startup_stm`、`stm_fig2_parameters`、`stm_fig2_eps_coeff`、`stpwe_build_system`、`stpwe_solve_omega`、`stpwe_static_bands`、`stpwe_select_mode`、`stm_fig2_epsilon`、`stpwe_reconstruct_field`、`stm_redblue` 及本文件局部函数 `save_demo_figure`。时间截断取 Mtime=1；频率上限 fMax=0.82，权重阈值 0.035，虚部容差 3e-3。约化坐标 `kBar = k/g`、`fBar = ω/(g·c0)`。

### `save_demo_figure`

**文件：** `demos/demo01_reproduce_fig2_stpwe.m`

**类型：** 函数（局部函数）

**函数签名：**

```matlab
function save_demo_figure(fig, outputFile)
```

**简介：** 将图形以 220 dpi 保存为 PNG 文件，并带打印回退。

**功能：**

优先调用 `exportgraphics` 以 220 dpi 导出图形到 `outputFile`；若导出失败则回退到 `print` 的 `-dpng -r220` 选项。这是一个用于演示脚本统一的图片保存辅助函数。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `fig` | matlab.ui.Figure | 待保存的图形句柄。 |
| `outputFile` | char | 输出 PNG 文件的完整路径。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无返回值，仅写文件。 |

——

### `demo02_complex_frequency_and_momentum_gaps`

**文件：** `demos/demo02_complex_frequency_and_momentum_gaps.m`

**类型：** 函数

**函数签名：**

```matlab
function demo02_complex_frequency_and_momentum_gaps()
```

**简介：** 对比 Park-Min 公式 (6) 与 (7)，展示动量带隙与频率带隙的复能谱。

**功能：**

构建与 demo01 相同的图 2 时空晶体系统。固定 k 扫描时用 `stpwe_solve_omega` 求解复本征频率 ω，筛出实部在 0.075~0.125 范围内的解，得到动量带隙中的复频率（增长/衰减对）；固定 ω 扫描时改用 `stpwe_solve_k` 求解复波矢 k，得到频率带隙中的消逝 Bloch 解。四个子图分别展示：固定 k 的实频率带隙、固定 k 的 Im(ω) 增长/衰减对、固定 ω 的实波矢带隙、固定 ω 的 Im(k) 消逝解。最后保存 PNG。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无输入参数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无返回值，仅绘图并保存 `output/demo02_complex_omega_and_k_gaps.png`。 |

**备注：** 依赖 `startup_stm`、`stm_fig2_parameters`、`stm_fig2_eps_coeff`、`stpwe_build_system`、`stpwe_solve_omega`、`stpwe_solve_k`。系统取 Nspace=10、Mtime=1。固定 k 扫描区间 kBar∈[0.12,0.23]，固定 ω 扫描区间 fBar∈[0.385,0.425]。

——

### `demo03_ptc_bands_pwe_vs_tmm`

**文件：** `demos/demo03_ptc_bands_pwe_vs_tmm.m`

**类型：** 函数

**函数签名：**

```matlab
function demo03_ptc_bands_pwe_vs_tmm()
```

**简介：** 用时间 PWE 与精确 2×2 D/B 单值矩阵（TMM）两种独立方法交叉验证二元光子时间晶体的能带。

**功能：**

设定二值光子时间晶体（ε_A=1.0、ε_B=4.0、占空比 0.5、周期 T=1）。一方面用 `temporal_crystal_bands`（精确 2×2 D/B 单值矩阵，TMM）计算准频率能带作为参考；另一方面用 `stpwe_build_system` 配 `temporal_binary_eps_coeff`（取 Mtime=19 的时间谐波截断）做时间 PWE，再用 `stpwe_fold_frequency` 将频率折叠到第一 Floquet 布里渊区，并通过最小化频率差与 m=0 权重的代价函数将 PWE 解逐支匹配到 TMM 能带。左右两个子图分别比较 Re(ω_F) 与 Im(ω_F)。第二个图绘制二值调制波形 ε(t)。最后保存 PNG。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无输入参数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无返回值，仅绘图并保存 `output/demo03_ptc_pwe_vs_tmm.png`。 |

**备注：** 依赖 `startup_stm`、`temporal_crystal_bands`、`temporal_binary_eps_coeff`、`stpwe_build_system`、`stpwe_solve_omega`、`stpwe_fold_frequency`。注释指出二值时域波形不连续，故 PWE 收敛比正弦调制慢得多，需 Mtime=19。第二段（绘制调制波形的图 2）大部分被注释掉。

——

### `demo04_temporal_multilayer_tmm`

**文件：** `demos/demo04_temporal_multilayer_tmm.m`

**类型：** 函数

**函数签名：**

```matlab
function demo04_temporal_multilayer_tmm()
```

**简介：** 用时间多层传递矩阵复现 Ramaccia 等论文的四类 APL 设计。

**功能：**

设定四组时间多层设计（任意堆叠、在 f₀ 透明、周期透明堆叠、时间放大），每组给出初始折射率 na、末折射率 nb、层折射率向量 n 与各层持续时间 dt，其中透明情形取 `dt = n/(2·na)`（δ=π）、放大情形取 `dt = n/(4·na)`（δ=π/2）。对每个 case 调用 `temporal_tmm_spectrum` 在归一化频率 f/f₀∈[0.5,1.5] 上计算前向 `forward` 与后向 `backward` 透射谱的幅值并绘图。最后保存 PNG。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无输入参数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无返回值，仅绘图并保存 `output/demo04_temporal_multilayer_tmm.png`。 |

**备注：** 依赖 `startup_stm`、`temporal_tmm_spectrum`。频率基准 ω₀=2π·fNorm（T₀=1）。ε 由折射率平方得到（εInitial=na²、εFinal=nb²、εSlabs=n.²），磁导率取 1。

——

### `demo05_fdtd_temporal_interface`

**文件：** `demos/demo05_fdtd_temporal_interface.m`

**类型：** 函数

**函数签名：**

```matlab
function demo05_fdtd_temporal_interface()
```

**简介：** 用 D/B Yee-FDTD 模拟单个时间边界处的波分裂，并与解析 Morgenthaler 系数对比。

**功能：**

以复高斯波包（中心 k₀、σ=3.5）为初值，在 t=tSwitch 处将折射率由 1.5 突变为 2.5（由局部函数 `temporal_eps` 定义）。用 `fdtd1d_db`（D/B 连续的 Yee 格式，海绵吸收边界）推进 Maxwell 方程，取开关前后的电场快照，通过 `E± = 0.5(E ± H/nAfter)` 做前后向分解，数值计算透射 `|τ|` 与反射 `|ρ|`，并与 `temporal_interface_matrix` 给出的解析 Morgenthaler 系数比较。四个子图展示：|E(x,t)| 时空图、开关后的方向分解、能量相对初值随时间变化。最后保存 PNG。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无输入参数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无返回值，仅绘图并保存 `output/demo05_fdtd_temporal_interface.png`。 |

**备注：** 依赖 `startup_stm`、`fdtd1d_db`、`temporal_interface_matrix` 及本文件局部函数 `temporal_eps`。空间步 dx=λ₀/32，时间步 dt=0.80·dx（CFL 稳定），海绵边界 160 格、强度 0.08。数值 τ、ρ 由 ℓ² 范数比值 `norm(E±)/norm(Ebefore)` 计算。

### `temporal_eps`

**文件：** `demos/demo05_fdtd_temporal_interface.m`

**类型：** 函数（局部函数）

**函数签名：**

```matlab
function epsr = temporal_eps(x,t,tSwitch,epsBefore,epsAfter)
```

**简介：** 返回时刻 t 的阶跃式相对介电常数分布（时间边界）。

**功能：**

当 t < tSwitch 时返回全空间等于 `epsBefore` 的向量，否则返回等于 `epsAfter` 的向量，用于模拟折射率在 t=tSwitch 时刻的突变。该函数被 `fdtd1d_db` 的 `epsFun` 回调调用。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `x` | double | 空间坐标向量。 |
| `t` | double | 当前时间（标量）。 |
| `tSwitch` | double | 折射率突变时刻。 |
| `epsBefore` | double | 突变前的相对介电常数。 |
| `epsAfter` | double | 突变后的相对介电常数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `epsr` | double | 与 `x` 同尺寸的相对介电常数向量。 |

——

### `demo06_fdtd_spacetime_wavepacket`

**文件：** `demos/demo06_fdtd_spacetime_wavepacket.m`

**类型：** 函数

**函数签名：**

```matlab
function demo06_fdtd_spacetime_wavepacket()
```

**简介：** 用 D/B Yee-FDTD 模拟高斯波包在图 2 时空晶体中的传播与频率转换。

**功能：**

以图 2 的时空调制介电常数 `stm_fig2_epsilon`（调制贯穿整个模拟域）为介质，用 `fdtd1d_db` 推进一个中心波矢 k₀=0.175g 的简单高斯波包。波包分解为可用的 Floquet-Bloch 模式，从而在探针谱中显现频率转换。取 x=14Λ 处的探针电场，用 `stm_stft` 做短时傅里叶变换得到 Floquet 边带谱。三个子图展示：|E(x,t)| 时空演化、探针谱图（边带）、参数能量交换（瞬时能量/初值）。最后保存 PNG。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无输入参数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无返回值，仅绘图并保存 `output/demo06_fdtd_spacetime_wavepacket.png`。 |

**备注：** 依赖 `startup_stm`、`stm_fig2_parameters`、`stm_fig2_epsilon`、`fdtd1d_db`、`stm_stft`。空间步 dx=Λ/40，时间步 dt=0.55·dx/c₀，模拟至 t=3T；有效折射率 nEff=√(0.75ε₁+0.25ε_c)，海绵边界 120 格。STFT 参数窗长 128、重叠 12、FFT 512。

——

### `demo07_zak_phase_stpwe`

**文件：** `demos/demo07_zak_phase_stpwe.m`

**类型：** 函数

**函数签名：**

```matlab
function demo07_zak_phase_stpwe()
```

**简介：** 计算连续 ST-PWE 能带在空间布里渊区上的 Zak 相位。

**功能：**

分别构造静态参考系统（`modDepth=0`）与驱动（Floquet）系统，用 `stpwe_track_band` 沿空间 Bloch 波矢 k 网格追踪单条能带（给出右/左本征矢 R、L），再用 `zak_phase_biorthogonal` 结合 `stpwe_bz_sewing_matrix` 提供的布里渊区缝合矩阵计算双正交 Wilson 回路，得到 Zak 相位与各段 Wilson 链接幅值。若驱动带的最小链接幅值小于 0.1，则发出警告提示单带 Zak 相位应视为收敛诊断而非稳健拓扑不变量。四个子图展示：追踪带、驱动带稳定性 Im(ω)、m=0 参与度、Wilson 链接幅值。最后保存 PNG。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无输入参数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无返回值，仅绘图并保存 `output/demo07_zak_phase_stpwe.png`。 |

**备注：** 依赖 `startup_stm`、`stm_fig2_parameters`、`stm_fig2_eps_coeff`、`stpwe_build_system`、`stpwe_track_band`、`zak_phase_biorthogonal`、`stpwe_bz_sewing_matrix`。此循环对空间 Bloch 动量 k 积分，不是空间均匀光子时间晶体的时间 Zak 相位（后者 Bloch 坐标为准频率、本征值为 k）。Zak 相位依赖原胞原点，图 2 原胞始于层边界而非反演中心，故未必严格为 0 或 π。

——

### `demo08_fhs_chern_thouless_pump`

**文件：** `demos/demo08_fhs_chern_thouless_pump.m`

**类型：** 函数

**函数签名：**

```matlab
function demo08_fhs_chern_thouless_pump()
```

**简介：** 用 Rice-Mele 泵模型对拓扑不变量算法做基准验证。

**功能：**

在 k×相位 二维网格上构造 Rice-Mele 二带哈密顿量 H=[mass, offDiag; conj(offDiag), -mass]，其中 `t1=t0+δcosφ`、`t2=t0-δcosφ`、`mass=m₀sinφ`、`offDiag=t1+t2·e^{-ik}`，逐点对角化取低能带本征矢并记录带隙。用 `fhs_chern_number` 计算 Fukui-Hatsugai-Suzuki 陈数，同时沿 k 回路求相位累积得到随调制相位 φ 变化的 Zak 相位并解绕（unwrap），展示量子化的极化泵浦。三个子图展示：FHS 曲率、解绕 Zak 相位、带隙。最后保存 PNG。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无输入参数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无返回值，仅绘图并保存 `output/demo08_fhs_chern_thouless_pump.png`。 |

**备注：** 依赖 `startup_stm`、`fhs_chern_number`。此例是拓扑例程的独立量化验证，不是对 demo01–07 连续 Maxwell 时空介质的陈数复现；调制相位 φ 作为合成维度。参数 Nk=Np=61、t0=1、δ=0.6、m₀=1.0。

——

### `demo09_finite_ptc_order_and_phase`

**文件：** `demos/demo09_finite_ptc_order_and_phase.m`

**类型：** 函数

**函数签名：**

```matlab
function demo09_finite_ptc_order_and_phase()
```

**简介：** 比较有限时长 PTC 在 AB 与 BA 两种层序下的可观测相位差。

**功能：**

AB 与 BA 时间原胞具有相同的无限晶体能带（其单值矩阵是循环置换），但有限实验还分辨出时间折射/反射输出的复相位，该相位随时间原点移动而变化。脚本用 `temporal_crystal_bands` 分别计算 AB（[ε_A,ε_B]）与 BA（[ε_B,ε_A]）能带并比较 Floquet 乘子差异，再用 `temporal_finite_crystal_response` 计算有限（nPeriods=8）晶体的前向/后向响应与相对相位，得到增益范数、相对相位及 AB/BA 相位差。四个子图展示：共同能带判别式 Tr(U)/2、有限响应增益、相对相位、序敏感的相位差。最后保存 PNG。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无输入参数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无返回值，仅绘图并保存 `output/demo09_finite_ptc_order_and_phase.png`。 |

**备注：** 依赖 `startup_stm`、`temporal_crystal_bands`、`temporal_finite_crystal_response`。参数 ε_A=3、ε_B=1、T=1、durations=[0.5 0.5]T、背景 ε=2。入射取 [1;0]（仅前向）。相位差用 `angle(exp(i·Δphase))` 包裹到 [-π,π]。

——

### `demo10_temporal_domain_wall_mode`

**文件：** `demos/demo10_temporal_domain_wall_mode.m`

**类型：** 函数

**函数签名：**

```matlab
function demo10_temporal_domain_wall_mode()
```

**简介：** 寻找时间畴壁（AB|BA 界面）处的局域界面模式。

**功能：**

两个层序相反的二元 PTC（AB 与 BA）具有相同的体带谱。在其共同动量带隙内，脚本用 `temporal_domain_wall_mode` 把 AB 的增长 Floquet 本征态与 BA 的衰减本征态匹配起来，包络在时间界面处达到峰值。四个子图展示：增长/衰减分支的 |λ| 共同动量带隙、Floquet 本征空间匹配度 |det(...)|、局域包络 |[D,B]|、以及 AB|BA 的 ε(t) 层序图。最后保存 PNG。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无输入参数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无返回值，仅绘图并保存 `output/demo10_temporal_domain_wall_mode.png`。 |

**备注：** 依赖 `startup_stm`、`temporal_domain_wall_mode`。参数 ε_A=3、ε_B=1、T=1、durations=[0.5 0.5]T、左右各 8 个周期。波矢扫描 kNorm∈[0.45,0.72]（k 以 2π/T 归一化）。

——

### `demo11_coherent_time_interface`

**文件：** `demos/demo11_coherent_time_interface.m`

**类型：** 函数

**函数签名：**

```matlab
function demo11_coherent_time_interface()
```

**简介：** 演示时间边界处双输入相干干涉与两种跃变定律的差异。

**功能：**

两束反向传播的相干波与一个全局时间开关碰撞，其相对相位控制两条出射时间分支。脚本用 `temporal_interface_matrix_jump` 分别以 `'DB'`（D/B 连续）与 `'EB'`（E/B 连续）两种理想化跃变定律构造时间界面矩阵，得到 τ、ρ 系数，并按输入振幅比 |τ/ρ| 扫描相对相位，计算前向/后向输出。四个子图展示：两定律的界面系数对比、DB 与 EB 输出幅值随相位变化（标出相干相消点）、两定律的振幅范数比（调制介导的放大/抑制）。最后保存 PNG。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无输入参数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无返回值，仅绘图并保存 `output/demo11_coherent_time_interface.png`。 |

**备注：** 依赖 `startup_stm`、`temporal_interface_matrix_jump`。材料跳变 εBefore=1.5²、εAfter=2.5²、μ=1。演示要点是：即便材料跳变相同，DB 连续与 EB 连续的理想化跃变定律给出的 τ、ρ 不同，说明微观开关机制必须显式声明。

——

### `demo12_ptc_convergence_audit`

**文件：** `demos/demo12_ptc_convergence_audit.m`

**类型：** 函数

**函数签名：**

```matlab
function demo12_ptc_convergence_audit()
```

**简介：** 以精确 2×2 单值矩阵为参考，审计不连续波形时间 PWE 的收敛性。

**功能：**

以 `temporal_crystal_bands` 给出的精确 2×2 单值矩阵能带为独立参考，对二元 PTC（ε_A=1、ε_B=4、占空比 0.5）用不同时间傅里叶截断 M ∈ {3,5,9,13,19} 运行时间 PWE（`stpwe_build_system` + `stpwe_solve_omega` + `stpwe_fold_frequency`），逐支匹配到参考能带并记录归一化谱误差 |Δω|/Ω 与运行时间。四个子图展示：中位数/最大谱误差随 M 收敛、各 k 处最差支误差热图、精确增长/衰减参考 Im(ω_F)、运行时间。同时把参数、M 序列、误差、运行时间保存为 MAT 文件。最后保存 PNG 与 MAT。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无输入参数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| 无 | — | 无返回值，仅绘图并保存 `output/demo12_ptc_convergence_audit.png` 与 `output/demo12_ptc_convergence_data.mat`（含 `parameters`、`Mvalues`、`errors`、`medianError`、`maximumError`、`runtime`）。 |

**备注：** 依赖 `startup_stm`、`temporal_crystal_bands`、`temporal_binary_eps_coeff`、`stpwe_build_system`、`stpwe_solve_omega`、`stpwe_fold_frequency`。k 扫描 kNorm∈[0.25,1.1]（k 以 Ω=2π/T 归一化）。有效解判定用 `abs(imag(folded)) < 0.45Ω`。该 MAT 文件是“误差、截断、参数与运行时间一并留存的论文配套审计”示例。
## examples/ — 10 个示例脚本（v1）

### `example_band_comparison`

**文件：** `examples/example_band_comparison.m`

**类型：** 函数

**函数签名：**

```matlab
function example_band_comparison(quality)
```

**简介：** 分别绘制三种能带结构（静态能带、空晶格近似、完整 Floquet 能带），便于逐幅比对时间调制的影响。

**功能：**

在同一参数配置下计算并分三张图绘制：图 1 为静态光子能带（无时间调制，仅空间周期性，对应原文黑色实线）；图 2 为空晶格近似（静态能带加上平移 ±Ω、±2Ω 的 umklapp 复制，对应黑色实线加灰色虚线）；图 3 为完整 Floquet 能带（ST-PWE 数值解，彩色散点）。三张图并排打开，用于对比 umklapp 折叠如何引入能带交叉、时间调制如何打开带隙。程序先用 `stpwe_build_system` 构造时空系统矩阵，逐 k 点调用 `stpwe_solve_omega` 求解本征频率，按有限值、正频率、频率上限 `fMax=0.82`、虚部容差 `3e-3`、m=0 参与度阈值 `0.035` 筛选模式；静态参考能带由 `stpwe_static_bands` 计算前 8 条。散点颜色编码 m=0 扇区参与度（品红=强调制，蓝=接近静态）。三张图分别保存为 PNG。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `quality` | char（可选） | 精度档位，`'quick'`（默认）或 `'paper'`；决定 `Nspace`（10 或 20）与 `Nk`（101 或 181）。默认 `'quick'`。 |

**输出参数：**

无（仅绘图/写文件）。保存三张 PNG 到 `output/example_band_comparison_{static_bands|empty_lattice|floquet_bands}_{quality}.png`。

**备注：** 依赖 `stm_init`、`stm_preset_modulated_slab`、`stm_fourier_modulated_slab`、`stpwe_build_system`、`stpwe_solve_omega`、`stpwe_static_bands` 及局部辅助函数 `save_example_figure`。归一化约定：`fBar = ω/(g*c0) = ωΛ/(2πc)`，`kBar = kΛ/(2π)`。

——

### `save_example_figure`（局部辅助函数）

**文件：** `examples/example_band_comparison.m`

**类型：** 函数

**函数签名：**

```matlab
function save_example_figure(fig, outputFile)
```

**简介：** 保存图像窗口到 PNG 文件，优先使用 `exportgraphics`，失败时回退到 `print`。

**功能：**

尝试用 `exportgraphics(fig, outputFile, 'Resolution',220)` 导出图形；若抛出异常则改用 `print(fig, outputFile, '-dpng', '-r220')`。两个后端均以 220 dpi 输出。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `fig` | matlab.ui.Figure | 待保存的图形窗口句柄。 |
| `outputFile` | char | 输出 PNG 的完整路径。 |

**输出参数：**

无（仅写文件）。

——

### `example_chern_pump`

**文件：** `examples/example_chern_pump.m`

**类型：** 函数

**函数签名：**

```matlab
function example_chern_pump()
```

**简介：** 用 Rice-Mele 泵验证 FHS Chern 数算法，计算时间调制二能带模型的量化拓扑不变量。

**功能：**

以 Rice-Mele 循环提供紧凑、独立量化的 Chern 数验证：调制相位作为合成第二维，与波矢 k 构成二维参数空间 (k, phase)，期望 Chern 数等于 1。程序在 61×61 的 (k, phase) 网格上构造 2×2 哈密顿量 `H = [mass, offDiag; conj(offDiag), -mass]`，其中 `t1 = t0 + delta0*cos(phase)`、`t2 = t0 - delta0*cos(phase)`、`mass = mass0*sin(phase)`、`offDiag = t1 + t2*exp(-i*k)`，对角化后取最低能带本征态，调用 `fhs_chern_number` 计算 Chern 数与曲率。另用 Wilson 环乘积相角计算随相位变化的 Zak 相位（极化泵）。绘图展示曲率、解缠绕 Zak 相位和能隙，保存为 PNG。

**输入参数：**

无。

**输出参数：**

无（仅绘图/写文件）。保存 `output/example_chern_pump.png`。命令行打印 FHS Chern 数。

**备注：** 依赖 `stm_init`、`fhs_chern_number`。参数 `t0=1`、`delta0=0.6`、`mass0=1.0`；`kGrid` 取 `-π..π`，`phaseGrid` 取 `0..2π`。

——

### `example_complex_gaps`

**文件：** `examples/example_complex_gaps.m`

**类型：** 函数

**函数签名：**

```matlab
function example_complex_gaps()
```

**简介：** 对比 ST-PWE 的动量禁带与频率禁带，展示两类复本征值问题。

**功能：**

演示 ST-PWE 的两个基本本征值问题：固定 k 求复频率 ω（揭示不稳定的动量禁带，对应式 6），以及固定 ω 求复波矢 k（揭示衰逝的频率禁带，对应式 7）。程序先构造时空系统矩阵，在 `kBar ∈ [0.12,0.23]` 扫描调用 `stpwe_solve_omega` 并按频率、虚部、m=0 权重筛选；再在 `fBar ∈ [0.385,0.425]` 扫描调用 `stpwe_solve_k` 并按 `Re(k)`、`Im(k)` 范围筛选。四张子图分别绘制：动量禁带的 Re(ω)（颜色= |Im(ω)|）、参数增长/衰减对 Im(ω)、频率禁带的 Re(k)（颜色= |Im(k)|）、以及衰逝 Bloch 解 Im(k)。保存为 PNG。

**输入参数：**

无。

**输出参数：**

无（仅绘图/写文件）。保存 `output/example_complex_gaps.png`。

**备注：** 依赖 `stm_init`、`stm_preset_modulated_slab`、`stm_fourier_modulated_slab`、`stpwe_build_system`、`stpwe_solve_omega`、`stpwe_solve_k`。归一化约定：`fBar = ω/(g*c0)`，`kBar = k/g`。

——

### `example_dual_sweep_bands`

**文件：** `examples/example_dual_sweep_bands.m`

**类型：** 函数

**函数签名：**

```matlab
function example_dual_sweep_bands(quality)
```

**简介：** 分别用「扫 k 求 ω」和「扫 ω 求 k」计算能带，并合并对比两种方法。

**功能：**

生成三张能带图：图 1 为固定 k 求 ω 的 2D 频带结构，仅画 Re(ω)，颜色编码 |Im(ω)|（蓝→红 = 带隙来临）；图 2 为固定 ω 求 k 的 3D 复动量能带结构，在 (Re(k), Im(k), ω) 空间中传播模落在 Im(k)=0 平面内、复模向 ±Im(k) 延伸，用于指示频率带隙；图 3 为叠加对比，把两种方法的实部解画在同一 (k, ω) 平面，重合处为传播模、分歧处为禁带。扫 k 使用 `stpwe_solve_omega`，按 fMax=0.82、虚部容差 `5e-2`、m=0 参与度 `0.035` 筛选；扫 ω 使用 `stpwe_solve_k`，按 `Re(k)∈[-0.5,0.5]`、`Im(k)` 容差 `0.12` 筛选。三张图分别保存为 PNG。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `quality` | char（可选） | 精度档位，`'quick'`（默认）或 `'paper'`；决定 `Nspace`、`Nk`、`Nw`。默认 `'quick'`。 |

**输出参数：**

无（仅绘图/写文件）。保存三张 PNG 到 `output/example_dual_sweep_bands_{fix_k_scan|fix_w_scan_3d|overlay}_{quality}.png`。

**备注：** 依赖 `stm_init`、`stm_preset_modulated_slab`、`stm_fourier_modulated_slab`、`stpwe_build_system`、`stpwe_solve_omega`、`stpwe_solve_k` 及局部辅助函数 `save_example_figure`。归一化：`fBar = ω/(g*c0)`，`kBar = k/g`。

——

### `save_example_figure`（局部辅助函数）

**文件：** `examples/example_dual_sweep_bands.m`

**类型：** 函数

**函数签名：**

```matlab
function save_example_figure(fig, outputFile)
```

**简介：** 保存图像窗口到 PNG 文件，优先使用 `exportgraphics`，失败时回退到 `print`。

**功能：**

尝试用 `exportgraphics(fig, outputFile, 'Resolution',220)` 导出图形；异常时改用 `print(fig, outputFile, '-dpng', '-r220')`。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `fig` | matlab.ui.Figure | 待保存的图形窗口句柄。 |
| `outputFile` | char | 输出 PNG 的完整路径。 |

**输出参数：**

无（仅写文件）。

——

### `example_fdtd_interface`

**文件：** `examples/example_fdtd_interface.m`

**类型：** 函数

**函数签名：**

```matlab
function example_fdtd_interface()
```

**简介：** 模拟单个突变时间界面上的波分裂，并与解析 Morgenthaler 系数对比。

**功能：**

一个复高斯波包在均匀介质中传播，在 `t = tSwitch` 处经历介电常数的突变（折射率 1.5 → 2.5），产生时间折射与时间反射。程序用 D/B 形式的 Yee 时域有限差分 `fdtd1d_db` 求解，随后在切换前后取样，通过方向分解 `Eplus = 0.5*(E + H/nAfter)`、`Eminus = 0.5*(E - H/nAfter)` 提取前向/后向散射振幅，与 `temporal_interface_matrix` 给出的解析透射/反射系数对比。绘图展示 |E(x,t)| 时空图、切换后的方向分解、以及瞬时能量变化，保存为 PNG。

**输入参数：**

无。

**输出参数：**

无（仅绘图/写文件）。保存 `output/example_fdtd_interface.png`。命令行打印数值与解析的 |τ|、|ρ| 对比。

**备注：** 依赖 `stm_init`、`fdtd1d_db`、`temporal_interface_matrix` 及局部辅助函数 `temporal_eps`。空间步长 `dx=λ0/32`，时间步长 `dt=0.80*dx`，边界条件为海绵吸收（sponge，160 格，强度 0.08）。

——

### `temporal_eps`（局部辅助函数）

**文件：** `examples/example_fdtd_interface.m`

**类型：** 函数

**函数签名：**

```matlab
function epsr = temporal_eps(x,t,tSwitch,epsBefore,epsAfter)
```

**简介：** 返回时间界面两侧的分段常数相对介电常数。

**功能：**

当 `t < tSwitch` 时返回 `epsBefore*ones(size(x))`，否则返回 `epsAfter*ones(size(x))`，实现介电常数在时间界面上的突变。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `x` | double | 空间网格坐标数组。 |
| `t` | double | 当前时刻标量。 |
| `tSwitch` | double | 切换时刻。 |
| `epsBefore` | double | 切换前的相对介电常数（本示例为 `nBefore^2`）。 |
| `epsAfter` | double | 切换后的相对介电常数（本示例为 `nAfter^2`）。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `epsr` | double | 与 `x` 同尺寸的相对介电常数数组。 |

——

### `example_fdtd_wavepacket`

**文件：** `examples/example_fdtd_wavepacket.m`

**类型：** 函数

**函数签名：**

```matlab
function example_fdtd_wavepacket()
```

**简介：** 对 Park-Min 时空晶体中的波包进行 FDTD 仿真，展示 Floquet 边带。

**功能：**

一个高斯波包穿过 Park-Min（Fig. 2）描述的时空周期介质，调制在整个仿真域内持续存在。程序用 D/B 形式的 Yee 时域有限差分 `fdtd1d_db` 求解（介质由 `stm_permittivity_modulated_slab` 定义），并在 `x=14Λ` 处设置探针，用短时傅里叶变换 `stm_stft` 计算探针谱图以揭示时空调制产生的 Floquet 边带。绘图展示波包演化的 |E(x,t)| 时空图、探针谱图（dB 标度）以及参数能量交换曲线，保存为 PNG。

**输入参数：**

无。

**输出参数：**

无（仅绘图/写文件）。保存 `output/example_fdtd_wavepacket.png`。

**备注：** 依赖 `stm_init`、`stm_preset_modulated_slab`、`stm_permittivity_modulated_slab`、`fdtd1d_db`、`stm_stft`。空间分辨率 40 点/原胞，时间步长 `dt=0.55*dx/c0`，总时长 `3*T`；载波 `k0=0.175*g`；边界为海绵吸收（120 格，强度 0.08）。

——

### `example_floquet_bands_and_fields`

**文件：** `examples/example_floquet_bands_and_fields.m`

**类型：** 函数

**函数签名：**

```matlab
function example_floquet_bands_and_fields(quality)
```

**简介：** 演示完整的 ST-PWE 工作流：Floquet 能带结构与本征模场分布。

**功能：**

在部分调制的原胞上演示完整 ST-PWE 流程：1) 由解析 Fourier 系数 `stm_fourier_modulated_slab` 通过 `stpwe_build_system` 构造系统矩阵；2) 扫 k 求解固定 k 本征问题 `stpwe_solve_omega` 得到 Floquet 能带，按有限值、正频率、fMax=0.82、虚部容差 3e-3、m=0 参与度 0.035 筛选；3) 用 `stpwe_select_mode` 在布里渊区选取 6 个代表性本征模 (c)-(h)；4) 用 `stpwe_reconstruct_field` 重构并绘制每个模的 E(x,t) 场分布。静态参考能带由 `stpwe_static_bands` 计算前 8 条。3×4 布局中绘制 ε(x,t) 曲面、能带图（叠加静态带与 umklapp 复制）及 6 幅场图，保存为 PNG。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `quality` | char（可选） | 精度档位，`'quick'`（默认）或 `'paper'`；决定 `Nspace`（10 或 20）与 `Nk`（101 或 181）。默认 `'quick'`。 |

**输出参数：**

无（仅绘图/写文件）。保存 `output/example_floquet_bands_and_fields_{quality}.png`。

**备注：** 依赖 `stm_init`、`stm_preset_modulated_slab`、`stm_fourier_modulated_slab`、`stm_permittivity_modulated_slab`、`stpwe_build_system`、`stpwe_solve_omega`、`stpwe_static_bands`、`stpwe_select_mode`、`stpwe_reconstruct_field`、`stm_redblue` 及局部辅助函数 `save_example_figure`。时间谐波截断 `Mtime=1`。颜色尺度用 m=0 Floquet 扇区参与度作为模式权重的透明代理。

——

### `save_example_figure`（局部辅助函数）

**文件：** `examples/example_floquet_bands_and_fields.m`

**类型：** 函数

**函数签名：**

```matlab
function save_example_figure(fig, outputFile)
```

**简介：** 保存图像窗口到 PNG 文件，优先使用 `exportgraphics`，失败时回退到 `print`。

**功能：**

尝试用 `exportgraphics(fig, outputFile, 'Resolution',220)` 导出图形；异常时改用 `print(fig, outputFile, '-dpng', '-r220')`。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `fig` | matlab.ui.Figure | 待保存的图形窗口句柄。 |
| `outputFile` | char | 输出 PNG 的完整路径。 |

**输出参数：**

无（仅写文件）。

——

### `example_ptc_pwe_vs_tmm`

**文件：** `examples/example_ptc_pwe_vs_tmm.m`

**类型：** 函数

**函数签名：**

```matlab
function example_ptc_pwe_vs_tmm()
```

**简介：** 对二元光子时间晶体，交叉验证平面波展开（PWE）与精确传输矩阵（TMM）两种方法。

**功能：**

用两种独立方法计算二元光子时间晶体的准频率能带：1) 时间 Fourier 平面波展开（PWE），取时间谐波截断 `Mtime=19`（因二元时间波形不连续，收敛慢于 Park-Min 的正弦调制）；2) 精确 2×2 D/B 单值矩阵（时间传输矩阵 `temporal_crystal_bands`）。对每个 k 点，PWE 解经 `stpwe_fold_frequency` 折叠到主区，并按「与 TMM 目标频率之差 + m=0 权重惩罚」的最小代价匹配两条分支。绘图对比两种方法的 Re(ω_F)（黑线 TMM、散点 PWE，颜色=m=0 权重）与 Im(ω_F)（动量禁带的增长/衰减），保存为 PNG。

**输入参数：**

无。

**输出参数：**

无（仅绘图/写文件）。保存 `output/example_ptc_pwe_vs_tmm.png`。

**备注：** 依赖 `stm_init`、`temporal_crystal_bands`、`temporal_binary_eps_coeff`、`stpwe_build_system`、`stpwe_solve_omega`、`stpwe_fold_frequency`。参数 `epsA=1`、`epsB=4`、`muA=muB=1`、占空比 `dutyA=0.5`、`T=1`；横轴为 `k*c0/Ω`（c0=1）。

——

### `example_tmm_multilayer`

**文件：** `examples/example_tmm_multilayer.m`

**类型：** 函数

**函数签名：**

```matlab
function example_tmm_multilayer()
```

**简介：** 复现 Ramaccia 等人的四类时间多层结构设计。

**功能：**

复现 Ramaccia et al., APL 118, 101901 (2021) 中的四类时间多层设计：1) 任意堆叠（Table I）；2) 在 f0 处透明（等传播距离，δ=π）；3) 周期透明堆叠；4) 时间放大堆叠（四分之一波，δ=π/2）。对每类结构，在 `fNorm ∈ [0.5,1.5]` 上调用 `temporal_tmm_spectrum` 计算前向/后向电场传递函数，并在 2×2 布局中绘制 |E_b^+/E_a^+| 与 |E_b^-/E_a^+| 随归一化频率的变化，保存为 PNG。

**输入参数：**

无。

**输出参数：**

无（仅绘图/写文件）。保存 `output/example_tmm_multilayer.png`。

**备注：** 依赖 `stm_init`、`temporal_tmm_spectrum`。频率设为 `ω0 = 2π*fNorm`（对应 T0=1）；每类结构用 `cases(q)` 结构体存储标题、入射/出射折射率 `na`/`nb`、层折射率 `n` 与层时间厚度 `dt`。

——

### `example_zak_phase`

**文件：** `examples/example_zak_phase.m`

**类型：** 函数

**函数签名：**

```matlab
function example_zak_phase()
```

**简介：** 计算 ST-PWE 连续能带的双正交 Zak 相位（一维 Wilson 环）。

**功能：**

对静态参考能带与 Park-Min 时空晶体的受驱动 Floquet 能带分别计算 Zak 相位（一维 Wilson 环）。静态系统设 `modDepth=0` 且 `Mtime=0`，受驱动系统用 `Mtime=1`；两者均用 `stpwe_track_band` 在 k 网格上追踪能带（返回右/左本征矢量 R、L），再用 `zak_phase_biorthogonal` 结合 `sys.Bomega` 与布里渊区缝合矩阵 `stpwe_bz_sewing_matrix` 计算双正交 Wilson 环。程序检查最小 Wilson 链接模，若低于 0.1 则提示能带近简并、需改用复合能带 Wilson 环。绘图展示追踪能带、驱动带稳定性、m=0 谱权重及 Wilson 链接模，保存为 PNG。

**输入参数：**

无。

**输出参数：**

无（仅绘图/写文件）。保存 `output/example_zak_phase.png`。命令行打印静态与驱动 Zak 相位（单位 π）及最小链接模。

**备注：** 依赖 `stm_init`、`stm_preset_modulated_slab`、`stm_fourier_modulated_slab`、`stpwe_build_system`、`stpwe_track_band`、`zak_phase_biorthogonal`、`stpwe_bz_sewing_matrix`。Zak 相位依赖原胞原点——Park-Min 原胞起始于层边界而非反演中心。种子频率 `omegaSeed=0.344*g*c0`；静态/驱动带宽窗口与 m=0 权重阈值由 `optionsStatic`/`optionsDriven` 设定。
## reproduction/Observation of temporal reflection — otr_* 核心函数

### `otr_abcd_unit_cell`

**文件：** `reproduction/Observation of temporal reflection/otr_abcd_unit_cell.m`

**类型：** 函数

**函数签名：**

```matlab
function [M,info] = otr_abcd_unit_cell(f,state,model,p)
```

**简介：** 计算单个开关传输线单元胞的 ABCD 矩阵（对应补充材料式 S1）。

**功能：**

该函数把单元胞建模为「半段无载传输线 + 并联支路 + 半段传输线」组成的对称 T 网络，逐频率返回 2×2×N 的 ABCD 矩阵。理想模型（`'ideal'`）的归一化导纳为 `yNorm = 1i*bScale*f`（即 b(f)=2.753·f/(100 MHz)）；realistic 模型把负载取为串联 RLC，`Zload = switchRon + Rpar + 1i*omega*Lpar + 1/(1i*omega*Cload)`，从而 `yNorm = Z0./Zload`。`'off'` 状态表示无并联负载，导纳为零。

返回的物理矩阵采用 `[V; I]` 约定（B 项乘 `Z0`、C 项除以 `Z0`），而 `info.normalizedMatrix` 是归一化 `[V; Z0*I]` 形式、严格等于 SI 式 S1。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `f` | double | 频率（Hz），元素必须为正的有限实数，可标量或数组 |
| `state` | char（可选） | 开关状态 `'off'`（无并联负载）或 `'on'`（默认 `'on'`） |
| `model` | char（可选） | 负载模型 `'ideal'` 或 `'realistic'`（默认 `'ideal'`） |
| `p` | struct（可选） | 参数结构，需含 Z0/Cload/switchRon/Rpar/Lpar/thetaScale/bScale 字段（默认 `otr_parameters()`） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `M` | double complex | 物理 ABCD 矩阵 `[V; I]`；标量 f 时为 2×2，否则 2×2×N |
| `info` | struct | 附加信息：`frequency`、`omega`、`theta`、`normalizedAdmittance`、`b`、`normalizedMatrix`（归一化 S1 矩阵）、`state`、`model`、`referenceConvention` |

**备注：** 依赖 `otr_parameters`。理想模型的 b(100 MHz)=2.753 是独立报告的设计导纳，并不等于 `omega*Cload*Z0`（后者属 82 pF 实际电路，在 100 MHz 处为 2.576）。

### `otr_bloch_dispersion`

**文件：** `reproduction/Observation of temporal reflection/otr_bloch_dispersion.m`

**类型：** 函数

**函数签名：**

```matlab
function [beta,ZB,isPass,info] = otr_bloch_dispersion(f,state,model,p)
```

**简介：** 计算周期开关传输线链的 Bloch 波数与 Bloch 阻抗（SI 式 S2–S3）。

**功能：**

该函数对由 `otr_abcd_unit_cell` 单元胞构成的周期链求解色散关系：先取归一化 ABCD 矩阵的迹半 `(A+D)/2`，由 `cos(beta*d) = (A+D)/2` 反解 Bloch 波数 `beta = acos(traceHalf)/d`（单位 rad/m），并按约定选取第一布里渊区内 `Im(beta)<=0`、`Re(beta)>=0` 的分支。Bloch 阻抗按 `ZB = B*Z0/sqrt(A^2-1)` 计算（式 S3），其中 B 为归一化矩阵元素，平方根分支取正实部（无耗通带内）。无耗时 `isPass` 依据 `|(A+D)/2|<=1` 判定通带；有耗单元则以 Bloch 衰减相对相位是否足够小来标记通带。

返回的 `beta` 为复数（rad/m）、`ZB` 为复数（欧姆），均保留输入 `f` 的形状。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `f` | double | 频率（Hz），正有限实数，可标量或数组 |
| `state` | char（可选） | `'off'` 或 `'on'`（默认 `'on'`） |
| `model` | char（可选） | `'ideal'` 或 `'realistic'`（默认 `'ideal'`） |
| `p` | struct（可选） | 参数结构（默认 `otr_parameters()`） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `beta` | double complex | Bloch 波数（rad/m） |
| `ZB` | double complex | Bloch 阻抗（欧姆），无解点取 NaN |
| `isPass` | logical | 是否处于第一通带的标记 |
| `info` | struct | 附加信息：`traceHalf`、`determinant`、`theta`、`b`、`normalizedAdmittance`、`state`、`model`、`firstBrillouinEdge`（=pi/d） |

**备注：** 依赖 `otr_abcd_unit_cell` 与 `otr_parameters`。采用 `exp(+j*omega*t)` / 前向衰减的时间因子约定。无耗时通带判据为 `|(A+D)/2|<=1`（含小容差）。

### `otr_build_mna`

**文件：** `reproduction/Observation of temporal reflection/otr_build_mna.m`

**类型：** 函数

**函数签名：**

```matlab
function model = otr_build_mna(options)
```

**简介：** 构建开关传输线的 Kirchhoff/MNA（改进节点分析）集总离散电路模型。

**功能：**

该函数构造一个离散、集总的电路模型（并非连续介质 FDTD 离散化），未知量为节点电压与物理支路电流/电压，每个矩阵行要么是 KCL 方程、要么是支路本构方程。默认梯形含 30 个串联 RL 单元、31 个电压节点，首节点接 Thevenin 源、末节点接 50 Ω 负载，节点 2:31 各有 30 个可切换的并联负载；单元的 L、C 由 `Z0`、`d` 与微带有效介电常数 `epsEff` 导出。返回的描述方程形如 `C(w)*xdot + G(w)*x = B*vSource(t)`，需配合 `otr_simulate_circuit` 组装 `C(w)`、`G(w)`、施加时间边界条件并积分。

`options.physicsMode` 支持两种模式：`'ideal-capacitance'`（瞬时电容插入/移除，闭合服从 `C_old*v- = C_new*v+`、断开服从 `v+ = v-`，即 SI 式 S14/S15）与 `'rlc'`（固定最大拓扑，每个开关节点带串联 `(Ron+Rpar)-Lpar-Cload` 支路，用连接权重耦合/解耦而不改变状态规模）。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `options` | struct（可选） | 可选字段（SI 单位）：`nCells`（默认 30）、`d`（0.2080 m）、`Z0`（50 Ω）、`epsEff`（8.36）、`c0`、`physicsMode`（`'ideal-capacitance'` 或 `'rlc'`）、`Lseries`、`Ccell`、`Rseries`（0.25）、`lossTangent`（0.0019）、`referenceFrequency`（100 MHz）、`Gshunt`、`sourceResistance`（50）、`loadResistance`（50）、`Cload`（82 pF）、`Ron`（1）、`Rpar`（4）、`Lpar`（8 nH）、`includeParasitics`（true）、`Cnearest`（2 pF）、`CnextNearest`（1 pF）、`switchNodes`（默认节点 2:nNodes） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `model` | struct | 电路模型：`name`、`form`、`physicsMode`、`parameters`、`nCells`、`nNodes`、`nLineBranches`、`nSwitches`、`nState`、`index`（各状态变量索引）、`incidence`、`CnodeBase`、`GnodePhysical`、`Cbase`（即 C 描述矩阵）、`Gbase`（即 G 描述矩阵）、`B`、`switchNodes`、`Cload`、`isMNA`、`notes` |

**备注：** Thevenin 源以精确 Norton 等效表示（`vSource/Rs` 并联 `Rs`）。寄生耦合用节点间真实电容的 Laplacian 戳记 `C*[1 -1; -1 1]` 加入节点电容矩阵。`index` 给出 `voltage`、`lineCurrent`、`switchCurrent`、`switchCapVoltage` 各状态分量的索引区间。

### `otr_circuit_energy`

**文件：** `reproduction/Observation of temporal reflection/otr_circuit_energy.m`

**类型：** 函数

**函数签名：**

```matlab
function energy = otr_circuit_energy(model, result)
```

**简介：** 审计 `otr_simulate_circuit` 结果的储能、损耗、源功与开关跳变功，验证能量平衡。

**功能：**

该函数对 `otr_simulate_circuit` 返回的结果做独立能量审计（SI 单位）。储能直接由物理电容与电感构造：电容能 `0.5*real(v'*Cnode*v)`、电感能 `0.5*sum(Lseries*iLine^2)`（RLC 模式再加 `Lpar`、`Cload` 支路储能）；耗散包含线/并联损耗、源内阻、50 Ω 负载以及适用时的 RLC 开关损耗。瞬时开关事件显式计入跳变功 `W_switch(t_e) = E(t_e^+) - E(t_e^-)`，其中 `E(t_e^-)` 用事件记录的 `stateMinus` 计算。事件之间，梯形 MNA 在功率取算术中点时具有精确的离散二次平衡，因此返回的残差也可作为求解器与边界条件的校验量。

最终给出 `balanceResidual = E - E0 - 连续功累计 - 开关功累计` 及归一化残差。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `model` | struct | 必须由 `otr_build_mna` 返回（含 `isMNA=true`） |
| `result` | struct | 必须由 `otr_simulate_circuit` 返回（含 time/source/weights/voltage/lineCurrent/eventRecords） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `energy` | struct | 能量审计结果：`time`、`stored`、`electric`、`magnetic`、`sourcePowerMid`、`dissipationPowerMid`、`loadPowerMid`、`sourceResistancePowerMid`、`lineLossPowerMid`、`shuntLossPowerMid`、`switchLossPowerMid`、`intervalContinuousWork`、`cumulativeContinuousWork`、`eventEnergyJump`、`cumulativeSwitchWork`、`balanceResidual`、`maxAbsBalanceResidual`、`relativeBalanceResidual`、`balanceEquation` |

**备注：** 负的跳变功对应电荷共享损耗或移除电容带走的能量。RLC 模式要求 `result` 含 `switchCurrent` 与 `switchCapVoltage`，否则报错。

### `otr_fft_spectrum`

**文件：** `reproduction/Observation of temporal reflection/otr_fft_spectrum.m`

**类型：** 函数

**函数签名：**

```matlab
function [f,X,info] = otr_fft_spectrum(t,x,varargin)
```

**简介：** 无工具箱、正确缩放的 FFT 频谱计算。

**功能：**

沿 `x` 的第一个非单例维做 FFT。`t` 需为秒单位的均匀采样向量，返回的频率 `f` 单位为 Hz。实数输入默认输出单边谱（内部对非直流/奈奎斯特分量乘 2 补偿），复数输入默认输出中心化的双边谱以保留负频反射分量。缩放除以窗口相干增益，使落在 bin 中心的单频正弦保持其峰值幅度。返回的 `X` 是复振幅谱（非功率谱），`info` 附有幅度、相位、功率与频率分辨率，便于绘图与提取，无需 Signal Processing Toolbox。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `t` | double | 时间向量（秒），有限实数、严格递增且均匀（至少 2 个样本） |
| `x` | double | 待变换信号数组，`size(x,Dimension)` 须等于 `numel(t)` |
| `varargin` | 名称/值对（可选） | `'Dimension'`（默认第一个非单例维）、`'NFFT'`（默认信号长度）、`'Window'`（`'none'`/`'hann'` 或长度-N 数值向量，默认 `'none'`）、`'Detrend'`（去均值，默认 false）、`'Spectrum'`（`'auto'`/`'onesided'`/`'twosided'`/`'centered'`，默认 `'auto'`） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `f` | double | 频率轴（Hz） |
| `X` | double complex | 复振幅谱（非功率谱） |
| `info` | struct | `magnitude`、`phase`、`power`、`sampleInterval`、`sampleRate`、`frequencyResolution`、`nfft`、`dimension`、`spectrum`、`coherentGain` |

**备注：** 单边缩放仅对实信号定义；复数信号配合 `'onesided'` 会报错。

### `otr_frequency_map`

**文件：** `reproduction/Observation of temporal reflection/otr_frequency_map.m`

**类型：** 函数

**函数签名：**

```matlab
function [fout,betaConserved,info] = otr_frequency_map(fin,fromState,toState,model,p,varargin)
```

**简介：** 动量守恒的 OFF↔ON 频率映射（固定 beta 求频率）。

**功能：**

在 `toState` 的第一通带内求解满足 `beta_TO(fout) = beta_FROM(fin)` 的频率，频率单位为 Hz。`fromState`/`toState` 为 `'off'`/`'on'`，`model` 为 `'ideal'`/`'realistic'`；有耗复数 beta 通过守恒其实部来映射，与 SI 第 S4 节的实验反演一致。实现上先在搜索区间上对 `toState` 色散曲线采样，取出第一通带、排序并去重，再对每个目标 beta 用 `interp1` 给出初值，随后用 `fzero`（有零点）或 `fminbnd`（无符号变化时求最小残差）精确求根。该函数不采用固定阻抗比，而是直接在 S2/S3 色散曲线第一通带内求根，从而复现文章 Fig. 2(e,f) 的宽带频率映射。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `fin` | double | 输入频率（Hz），正有限实数，须全部位于 fromState 第一通带 |
| `fromState` | char（可选） | `'off'` 或 `'on'`（默认 `'off'`） |
| `toState` | char（可选） | `'off'` 或 `'on'`（默认 `'on'`） |
| `model` | char（可选） | `'ideal'` 或 `'realistic'`（默认 `'ideal'`） |
| `p` | struct（可选） | 参数结构（默认 `otr_parameters()`） |
| `varargin` | 名称/值对（可选） | `'SearchRange'`（默认 `[1 kHz, 2.5*fReference]`）、`'GridSize'`（默认 6001） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `fout` | double | 输出频率（Hz），与 `fin` 同形状 |
| `betaConserved` | double complex | 守恒的输入 beta（即 `bIn`） |
| `info` | struct | `inputFrequency`、`outputFrequency`、`ratio`、`fromState`、`toState`、`model`、`targetBeta`、`residualBeta`、`targetFirstPassBand`、`searchRange` |

**备注：** 依赖 `otr_bloch_dispersion` 与 `otr_parameters`。若某输入频率的 beta 在目标第一通带之外会报错（无解）。

### `otr_gate_signal`

**文件：** `reproduction/Observation of temporal reflection/otr_gate_signal.m`

**类型：** 函数

**函数签名：**

```matlab
function gated = otr_gate_signal(time, signal, limits, taperFraction)
```

**简介：** 提取、基线校正并软门控一段波形。

**功能：**

在请求的闭区间 `limits` 内选取样本，去除一条连接两端稳健（截尾）均值的直线以做基线校正，然后在两端各施加余弦锥形窗做软门控。该操作透明、不做隐藏滤波，也无需 Signal Processing Toolbox。`taperFraction` 设 0 时退化为矩形门控。返回结构体含原始、基线、去趋势、窗及最终信号等字段，便于后续 FFT 分析。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `time` | double | 时间向量，须与 `signal` 等长 |
| `signal` | double | 信号向量 |
| `limits` | double | 两元素 `[T0 T1]`，门控区间 |
| `taperFraction` | double（可选） | 两端锥形窗占比，范围 0~0.5（默认 0.10） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `gated` | struct | `time`、`raw`、`baseline`、`detrended`、`window`、`signal`（=detrended.*window）、`limits`、`dt` |

### `otr_parameters`

**文件：** `reproduction/Observation of temporal reflection/otr_parameters.m`

**类型：** 函数

**函数签名：**

```matlab
function p = otr_parameters(varargin)
```

**简介：** 返回开关传输线实验的参数（论文/Methods/SI 值，SI 单位）。

**功能：**

返回文章正文、Methods 与补充材料中报告的参数结构体，全部采用 SI 单位，并派生出后续计算需要的量（相速度、线电导/电感、单元 C/L、`thetaScale`、`bScale`、负载电阻与设计阻抗比等）。可传入 name/value 对覆盖单个字段，或传入一个结构体批量覆盖。覆盖后重新计算全部派生量。无需任何工具箱，支持 MATLAB R2019b 及以后版本。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `varargin` | 名称/值对 或 struct（可选） | 覆盖若干标量字段；或一个标量结构体应用其全部字段（派生字段被忽略并重算） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `p` | struct | 参数结构：原始量 `c0`、`d`、`Z0`、`epsEff`、`Cload`、`switchRon`、`Rpar`、`Lpar`、`nCells`、`Cnearest`、`CnextNearest`、`theta100`、`b100`、`switchingTime`、`fReference`、`alpha100`、`lossTangent`；派生量 `phaseVelocity`、`phaseVelocityMethods`、`ClinePerLength`、`LlinePerLength`、`ClineCell`、`LlineCell`、`sampleLength`、`thetaScale`、`bScale`、`Rload`、`ZonDesign`、`frequencyRatioDesign` |

**备注：** `thetaScale`/`bScale` 供 S1–S3 使用（基于报告的 72° 设计相位），而物理 L'/C' 值使用 `epsEff` 的微带估计（与 MNA 实现一致）。`ZonDesign=Z0/2` 对应 SI 中 50 MHz 附近的设计目标。

### `otr_read_numeric_table`

**文件：** `reproduction/Observation of temporal reflection/otr_read_numeric_table.m`

**类型：** 函数

**函数签名：**

```matlab
function tbl = otr_read_numeric_table(filename, delimiter)
```

**简介：** 稳健读取小型数值文本/CSV 表。

**功能：**

跳过空行与元数据行，直到找到至少含两个数值字段的行作为数值数据起点；紧邻其前、宽度匹配的非数值行被保留为表头。自动去除 UTF-8 BOM（以 Unicode U+FEFF 或三个原始 UTF-8 字节出现）。数据按最长行宽填 NaN 对齐，并删除全空列。该读取器有意避免 `readtable`/`detectImportOptions`（其推断规则随版本变化），也不需要任何导入或信号处理工具箱。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `filename` | char | 输入文件路径 |
| `delimiter` | char（可选） | 分隔符；缺省时自动检测为 tab/comma/semicolon/whitespace |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `tbl` | struct | `data`（数值矩阵，短行补 NaN）、`headers`（列标签 cell）、`metadata`（数据前的原始文本行）、`filename`、`delimiter` |

### `otr_read_scope_csv`

**文件：** `reproduction/Observation of temporal reflection/otr_read_scope_csv.m`

**类型：** 函数

**函数签名：**

```matlab
function scope = otr_read_scope_csv(filename)
```

**简介：** 读取带内嵌元数据的 Tektronix 示波器 CSV 导出。

**功能：**

通过定位 `TIME`/`CH*` 表头（而非假定固定元数据长度）来解析 Tektronix CSV。支持论文 Fig. 2c 源数据这类导出：右侧在空分隔列之后还放置第二块 `TIME`/`MATH`（频谱）数据。用 `textscan` 高效读取（论文导出约 10 万行），去除 BOM，返回时域主块与可选的右侧频谱块。无需任何工具箱函数，兼容 MATLAB R2019b+。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `filename` | char | Tektronix CSV 文件路径 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `scope` | struct | `time`、`channels`、`channelNames`（主时域块）；`frequency`、`spectra`、`spectrumNames`（可选右侧块）；`sampleInterval`、`metadata`、`rawData`、`rawHeaders`、`filename` |

### `otr_save_figure`

**文件：** `reproduction/Observation of temporal reflection/otr_save_figure.m`

**类型：** 函数

**函数签名：**

```matlab
function filename = otr_save_figure(fig, basename, resolution)
```

**简介：** 导出单张复现图，带可移植回退方案。

**功能：**

调用 `otr_setup()` 获取输出目录，把图保存为 `<output>/<basename>.png`。优先用 `exportgraphics` 导出，若失败（如旧版本不支持）则回退到 `print` 的 `-dpng` 方式，并在命令行打印保存路径。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `fig` | figure 句柄 | 要保存的图形对象 |
| `basename` | char | 输出文件名（不含扩展名） |
| `resolution` | double（可选） | 分辨率 DPI（默认 220） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `filename` | char | 保存的完整 PNG 路径 |

**备注：** 依赖 `otr_setup`。

### `otr_setup`

**文件：** `reproduction/Observation of temporal reflection/otr_setup.m`

**类型：** 函数

**函数签名：**

```matlab
function paths = otr_setup()
```

**简介：** 将本复现目录及父工具箱加入 MATLAB 路径，并返回相关目录路径。

**功能：**

可从任意工作目录安全调用。定位脚本所在目录与父工具箱根目录；若父工具箱的核心函数（`fdtd1d_db`、`temporal_crystal_monodromy`）尚不可用且根目录存在 `startup_stm.m`，则加入根目录并执行启动脚本。随后把本复现目录加到路径最前，保证论文专属同名 helper 优先。所有生成的图形与 MAT 文件均保存在本复现目录下，不修改目录之外的源文件。

**输入参数：**

无

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `paths` | struct | `root`（父工具箱根目录）、`reproduction`（本复现目录）、`data`（`data/Time_Interface_Source_Data`）、`output`（输出目录，不存在则创建） |

### `otr_simulate_circuit`

**文件：** `reproduction/Observation of temporal reflection/otr_simulate_circuit.m`

**类型：** 函数

**函数签名：**

```matlab
function result = otr_simulate_circuit(model, t, source, schedule, options)
```

**简介：** 用梯形法则积分离散开关线 MNA 方程。

**功能：**

求解 `C(w)*xdot + G(w)*x = B*source(t)`，其中 `model` 由 `otr_build_mna` 返回、`source` 可为函数句柄/标量/每时间样本一个值、`schedule` 由 `otr_switch_schedule` 返回（也可直接传入规格结构自动转换）。积分约定显式：每个区间用其左端时刻的状态积分，开关事件精确落在网格右端，故 `result.state(:,n)=x(t_n^+)`，事件记录同时保留 `x(t_n^-)`。ideal-capacitance 模式下，权重下降先以电压连续移除电容，权重上升再施加 `C_before*v^- = C_after*v^+`，电感电流始终连续，有限上升沿被离散为一串电荷共享事件；RLC 模式下权重耦合固定的内部 RLC 拓扑，OFF 事件隔离并清零被移除支路。默认在每次边界处验证物理残差，超限报错。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `model` | struct | 必须由 `otr_build_mna` 返回 |
| `t` | double | 严格递增的有限实数时间向量（至少 2 个样本） |
| `source` | function_handle / double | 源激励：函数句柄、标量、或每个时间样本一个值 |
| `schedule` | struct | `otr_switch_schedule` 的返回（或规格结构） |
| `options` | struct（可选） | `initialState`（默认全零）、`method`（仅 `'trapezoidal'`）、`boundaryTolerance`（默认 5e-10）、`verifyBoundary`（默认 true）、`storeState`（默认 true） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `result` | struct | `time`、`source`、`weights`、`voltage`、`lineCurrent`、`inputVoltage`（节点 1）、`outputVoltage`（末节点）、`eventRecords`、`schedule`、`method`、`modelPhysicsMode`、`boundaryConvention`、`switchCurrent`、`switchCapVoltage`、`state`（若 `storeState`） |

**备注：** 依赖 `otr_build_mna` 与 `otr_switch_schedule`。`result.boundaryConvention = 'state(:,n)=x(t_n^+)'` 说明了状态时间约定。

### `otr_switch_schedule`

**文件：** `reproduction/Observation of temporal reflection/otr_switch_schedule.m`

**类型：** 函数

**函数签名：**

```matlab
function schedule = otr_switch_schedule(t, specification)
```

**简介：** 生成开关连接权重与精确事件元数据。

**功能：**

在严格递增时间网格 `t` 上采样空间均匀（或逐开关）的切换律。事件时间必须落在网格上（到舍入精度），否则报错，从而避免在积分步之间静默施加时间边界。`specification` 可为标量结构（字段 `eventTimes`、`states`、`nSwitch`、`riseTime`、`transition`、`snapTolerance`；便捷字段 `initialState`/`eventStates` 可替代 `states`），也可用字符串 `mode`（`'off-on'`、`'on-off'`、`'slab'`、`'inverse-slab'`）。有限 `riseTime` 时，过渡从事件时刻开始、历经 `riseTime` 结束，且各事件过渡不得重叠。返回 `weights`（nSwitch×numel(t)）及 `eventIndices` 等事件元数据。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `t` | double | 严格递增的有限实数时间向量（至少 2 个样本） |
| `specification` | struct（可选） | `eventTimes`（递增向量）、`states`（长度 nEvents+1 的向量或 nSwitch×(nEvents+1) 矩阵）、`nSwitch`（默认 1）、`riseTime`（默认 0）、`transition`（`'linear'`/`'smoothstep'`，默认 `'linear'`）、`snapTolerance`（可选）、`initialState`、`eventStates`、`mode`（字符串） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `schedule` | struct | `time`、`weights`、`states`、`eventTimesRequested`、`eventTimes`、`eventIndices`、`riseTime`、`transition`、`nSwitch`、`isInstantaneous`、`eventType`（每事件方向 ±1） |

### `otr_synthetic_wavepacket`

**文件：** `reproduction/Observation of temporal reflection/otr_synthetic_wavepacket.m`

**类型：** 函数

**函数签名：**

```matlab
function out = otr_synthetic_wavepacket(options)
```

**简介：** 时间界面下波包散射的解析谱解法（SI S5 的理想空间均匀电路极限）。

**功能：**

对理想、空间均匀电路极限做谱解法（对应 SI S5），与离散 Kirchhoff/MNA 电路求解器互补而非替代。每个空间傅里叶分量守恒波数 k：用 FFT 把初始高斯波包分解到 k 域，按界面前后相速度给各分量分配频率 `w = v*|k|`，再用 `otr_temporal_interface` 的透射/反射系数在时间界面（或时间板）处组合，最后 IFFT 回到空间域。支持 `'on'`（OFF→ON）、`'off'`（ON→OFF）、`'slab'`（双界面时间板）三种方向。返回前向、后向及总场。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `options` | struct（可选） | `direction`（默认 `'on'`）、`Zoff`/`Zon`（默认 50/25 Ω）、`fInitial`（默认 50 MHz）、`x`、`t`、`x0`、`sigma`、`tSwitch`（默认 40 ns）、`tau`（默认 25 ns）、`initialVelocity`（默认 1.04e8 m/s） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `out` | struct | `x`、`t`、`field`、`forward`、`backward`、`initial`、`k`、`k0`、`omegaInitial`、`omegaOther`、`fRatio`、`Tfirst`、`Rfirst`、`tSwitch`、`tau`、`direction` |

**备注：** 内部 `otr_temporal_interface_local` 优先调用 `otr_temporal_interface`，若文件不存在则用本地解析公式回退。时间因子约定为 `exp(-1i*w*t)`。

### `otr_temporal_interface`

**文件：** `reproduction/Observation of temporal reflection/otr_temporal_interface.m`

**类型：** 函数

**函数签名：**

```matlab
function [J,T,R,info] = otr_temporal_interface(Zbefore,Zafter,law,varargin)
```

**简介：** 用 SI 微观定律计算时间界面散射。

**功能：**

返回幅度映射 `[V_after+; V_after-] = J*[V_before+; V_before-]`。`law='charge'` 对应连接电容的 `C2*v(0+)=C1*v(0-)`（式 S14），系数用式 S16（'on'）；`law='voltage'` 对应移除电容的 `v(0+)=v(0-)`（式 S15），系数用式 S17（'off'）；`'auto'` 依据 `|Zafter|<|Zbefore|` 自动选 charge，否则选 voltage。别名 `'on'/'S16'` 与 `'off'/'S17'` 均可接受。`T=J(1,1)`、`R=J(2,1)`，J 关于正负频入射幅对称。色散计算中 `Zbefore`/`Zafter` 可为等尺寸数组，J 为 2×2×N。SI 恒等式假设非磁性传输线，故 `C_before/C_after=(Z_after/Z_before)^2`。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `Zbefore` | double complex | 界面之前阻抗（正实部），标量或等尺寸数组 |
| `Zafter` | double complex | 界面之后阻抗（正实部），标量或等尺寸数组 |
| `law` | char（可选） | `'charge'`/`'voltage'`/`'auto'`（默认 `'auto'`）及别名 |
| `varargin` | 名称/值对（可选） | `'CapacitanceRatio'`（Cbefore/Cafter）可覆盖该比值 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `J` | double complex | 时间散射矩阵（2×2 或 2×2×N） |
| `T` | double complex | 透射系数 `J(1,1)` |
| `R` | double complex | 反射系数 `J(2,1)` |
| `info` | struct | `law`、`direction`、`equation`、`Zbefore`、`Zafter`、`capacitanceRatio`、`frequencyRatio`、`boundaryScale`、`difference` |

**备注：** 对于 Zbefore=50、Zafter=25 的 charge 定律给出 T=0.375、R=-0.125；Zbefore=25、Zafter=50 的 voltage 定律给出 T=1.5、R=-0.5。

### `otr_temporal_slab`

**文件：** `reproduction/Observation of temporal reflection/otr_temporal_slab.m`

**类型：** 函数

**函数签名：**

```matlab
function [S,paths] = otr_temporal_slab(Zouter,Zslab,tau,omegaSlab,varargin)
```

**简介：** 两界面时间板及其四条路径散射。

**功能：**

对实验的 OFF-ON-OFF 时间板建模：第一界面用 S14/S16（电荷守恒）、第二界面用 S15/S17（电压守恒）。`omegaSlab` 是板内角频率（rad/s）、`tau` 是板持续时长（秒）。总散射矩阵 `S = Joff*diag(exp(+j*phi),exp(-j*phi))*Jon`，其中 `phi=omegaSlab*tau`。`paths` 给出 SI 式 S17 之后引用的四条基本幅度：`TT=Ton*Toff*exp(+j*phi)`、`RR=Ron*Roff*exp(-j*phi)`、`RT=Ron*Toff*exp(-j*phi)`、`TR=Ton*Roff*exp(+j*phi)`，于是 `S(1,1)=TT+RR`（时间折射）、`S(2,1)=RT+TR`（时间反射）。标量阻抗可与任意等尺寸频率数组搭配。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `Zouter` | double complex | 板外阻抗（标量或等尺寸数组） |
| `Zslab` | double complex | 板内阻抗（标量或等尺寸数组） |
| `tau` | double | 板持续时间（秒，非负标量） |
| `omegaSlab` | double | 板内角频率（rad/s，非负） |
| `varargin` | 名称/值对（可选） | `'FirstLaw'`、`'SecondLaw'` 覆盖默认的 charge/voltage 微观定律（支持倒置 ON-OFF-ON 对照实验） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `S` | double complex | 总散射矩阵（2×2 或 2×2×N） |
| `paths` | struct | `phase`、`TT`、`RR`、`RT`、`TR`、`timeRefracted`、`timeReflected`、`firstMatrix`、`secondMatrix`、`firstInterface`、`secondInterface`、`identityResidualForward`、`identityResidualReflected` |

**备注：** 依赖 `otr_temporal_interface`。返回的 `identityResidual*` 用于核对路径分解与矩阵乘法的一致性。

### `otr_validate`

**文件：** `reproduction/Observation of temporal reflection/otr_validate.m`

**类型：** 函数

**函数签名：**

```matlab
function report = otr_validate()
```

**简介：** 运行 OTR 复现的确定性自动验收测试。

**功能：**

运行一系列无工具箱的确定性检查，覆盖解析传输线模型、时间边界条件、时间板、FFT 工具、官方 source data 以及一个小规模 MNA 电路，包括：SI S1 的精确单元胞矩阵与幺模性、S2/S3 的 Bloch 迹恒等式与 50 MHz 阻抗目标、守恒 beta 频率映射、S14–S17 微观开关定律及四条系数恒等式、时间板四路径分解与首个反射零点、FFT 的实振幅与负频峰值、官方源数据文件可读性、MNA 的 ON/OFF 状态重映射与能量平衡。任何失败的物理/数值不变量都会触发 `assert`，全部通过后把报告写入 `output/validation_report.mat` 与 `output/validation_summary.txt`（位于本复现目录旁）。

**输入参数：**

无

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `report` | struct | 各检查子结构 `s1`、`s2s3`、`frequencyMap`、`temporalInterface`、`temporalSlab`、`fft`、`sourceData`、`mna`，以及 `name`、`timestamp`、`matlabVersion`、`status`、`parameters`、`elapsedSeconds`、`allPassed`、`outputFiles` |

**备注：** 依赖本组几乎全部函数（`otr_setup`、`otr_parameters`、`otr_abcd_unit_cell`、`otr_bloch_dispersion`、`otr_frequency_map`、`otr_temporal_interface`、`otr_temporal_slab`、`otr_fft_spectrum`、`otr_read_numeric_table`、`otr_build_mna`、`otr_switch_schedule`、`otr_simulate_circuit`、`otr_circuit_energy`）。可自任意工作目录直接运行；完整复现还要求 `source_data.zip` 及其解压出的 `data/Time_Interface_Source_Data` 目录存在。
## reproduction/Observation of temporal reflection — fig* 复现脚本

### `demo_circuit_time_interface`

**文件：** `reproduction/Observation of temporal reflection/demo_circuit_time_interface.m`

**类型：** 函数

**函数签名：**

```matlab
function results = demo_circuit_time_interface()
```

**简介：** 用 Kirchhoff/MNA 集总电路模型复现论文的时域反射实验（不调用连续介质 FDTD 求解器）。

**功能：**

用 30 节串联 RL、并联电容（含寄生电容）的传输线梯形网络逼近论文中的弯折线，所有并联负载在单一网格对齐时刻瞬时接入，构成一次时间界面。以两个高斯导数（derivative-of-Gaussian）小波构造“小脉冲先发、大脉冲后发”的非对称宽带入射，利用零直流特性避免缓慢给旁路电容充电。随后做三组独立自检：单节电路的电荷守恒（S14）与电压连续（S15）边界条件单元测试、MNA 中点功率能量平衡残差、以及近输入端探针处的时间反演排序与负极性符号相关检验。最后绘制节点电压时空图、源/探针波形与能量曲线，并按固定 Bloch 波数 k 计算集总梯形色散的频率下移，保存 PNG/FIG/MAT 文件。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无输入参数，全部配置在函数体内硬编码（`nCells=30`、`Z0=50`、`epsEff=8.36`、`Cload=82 pF` 等）。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `results` | struct | 汇总结构体，含 `parameters`（电路参数）、`time`、`source`、`probeNode`、`probe`、`switchTime`、`frequencyRatioLumped`、`signedTRCorrelation`、`energyBalanceResidual`、`onChargeResidual`、`offVoltageJump`、`k`、`frequencyOff`、`frequencyOn`；同时写入 `circuit_time_interface.png`/`.fig`、`circuit_fixed_k_frequency_translation.png`/`.fig` 与 `demo_circuit_time_interface.mat`。 |

**备注：** 依赖 `otr_setup`、`otr_build_mna`、`otr_switch_schedule`、`otr_simulate_circuit`、`otr_circuit_energy`。切换时刻须落在时间网格上（`assert` 校验）；能量平衡用梯形规则，残差按四舍五入量级期望。开关 ON（增大电容）时返回的时间反射具有负极性符号相关。

### `local_peak_near`

**文件：** `reproduction/Observation of temporal reflection/demo_circuit_time_interface.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function peaks = local_peak_near(time, signal, centres, halfWidth)
```

**简介：** 在给定中心时刻附近窗口内取信号绝对值最大的采样值。

**功能：** 对每个中心时刻，在 `[centre-halfWidth, centre+halfWidth]` 窗口内截取信号并取 `abs` 最大点，返回该点的原始（带符号）值，用于提取时间反射峰值。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `time` | double 向量 | 时间轴（秒） |
| `signal` | double 向量 | 与 `time` 等长的信号 |
| `centres` | double 向量 | 各窗口中心时刻 |
| `halfWidth` | double 标量 | 窗口半宽（秒） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `peaks` | double 向量 | 每个窗口内绝对值最大的原始采样值 |

### `fig01_source_data`

**文件：** `reproduction/Observation of temporal reflection/fig01_source_data.m`

**类型：** 函数

**函数签名：**

```matlab
function result = fig01_source_data()
```

**简介：** 用官方源数据复现论文 Fig. 1e/f（色散提取与入射/折射波形）。

**功能：**

读取开关 OFF/ON 的测量与 ADS 仿真相位表、示波器 CSV（`tek0000.csv`）以及 ADS 输入/输出波形，将相位按 $\beta = -\phi(\text{rad})/L$（总长 $L=30\times0.208$ m）转换为传播常数。在 10–70 MHz 频带内用最小二乘斜率估计 OFF/ON 各自的群速度。以 60 MHz 输入做守恒 $\beta$ 检验，得到 ON 态对应的输出频率与频率转换比。将测量时间轴平移使入射峰值与 ADS 对齐后，绘制色散双支曲线与输入/输出端口波形，输出 PNG 与 MAT。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无输入参数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `result` | struct | 含 `description`、`dataDirectory`、`files`、`totalLengthM`、`dispersion`（四组频率/相位/$\beta$）、`scope`、`adsInput`、`adsOutput`、`measuredTimeShiftS`、`fitBandMHz`、`groupVelocityMPerS`、`translationCheck`、`outputFiles`；写入 `fig01_source_data.png` 与 `fig01_source_data.mat`。 |

**备注：** 依赖 `otr_read_numeric_table`、`otr_read_scope_csv`。源数据目录须存在，否则抛 `otr:MissingSourceData` / `otr:DataDirectoryMissing`。测量时间轴仅整体平移以对齐入射峰值，不做重采样。

### `local_phase_to_beta`

**文件：** `reproduction/Observation of temporal reflection/fig01_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function out = local_phase_to_beta(tbl, totalLength)
```

**简介：** 将相位表转换为按频率排序的传播常数 $\beta$。

**功能：** 取数值表前两列（频率、相位度），剔除非法值并排序，按 $\beta=-\phi_{\deg}\pi/180/L$ 计算传播常数，返回含 `frequencyHz`、`phaseDeg`、`beta`、`source` 的结构体。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `tbl` | struct | 数值表结构体，含 `data` 与 `filename` 字段 |
| `totalLength` | double 标量 | 传输线总长度（米） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `out` | struct | `frequencyHz`（Hz）、`phaseDeg`（度）、`beta`（rad/m）、`source` |

### `local_group_velocity`

**文件：** `reproduction/Observation of temporal reflection/fig01_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function velocity = local_group_velocity(dispersion, bandMHz)
```

**简介：** 在指定频带内用最小二乘拟合 $\beta(\omega)$ 的斜率估计群速度。

**功能：** 选取频带内 $\beta>0$ 的采样点，用 $\omega=2\pi f$ 与 $\beta$ 计算线性最小二乘斜率 $v_g=\frac{\beta^\top\omega}{\beta^\top\beta}$（过原点拟合），得到该频带内的平均群速度。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `dispersion` | struct | 含 `frequencyHz` 与 `beta` 字段 |
| `bandMHz` | double 1×2 | 拟合频带 `[f_low, f_high]`（MHz） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `velocity` | double 标量 | 群速度（m/s） |

### `local_interp_unique`

**文件：** `reproduction/Observation of temporal reflection/fig01_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function value = local_interp_unique(x, y, query)
```

**简介：** 对去重排序后的 `x` 做线性插值查询。

**功能：** 排序并去重 `x`，同步重排 `y`，再对查询点做线性插值，用于从色散曲线反查频率。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `x` | double 向量 | 自变量（如 $\beta$） |
| `y` | double 向量 | 因变量（如频率） |
| `query` | double 标量 | 查询点 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `value` | double 标量 | 插值结果 |

### `local_plot_dispersion`

**文件：** `reproduction/Observation of temporal reflection/fig01_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function local_plot_dispersion(ax, d, color, style, measured)
```

**简介：** 在坐标轴上绘制色散曲线及其关于原点的镜像。

**功能：** 选取 0–100 MHz 且 $\beta$ 有限的点，绘制 $\beta$–f 及其镜像 $(-\beta,-f)$；`measured=true` 时用点标记抽稀显示，否则用实线。无返回值，仅绘图。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `ax` | matlab.graphics.axis.Axes | 目标坐标轴 |
| `d` | struct | 含 `frequencyHz`、`beta` 字段 |
| `color` | double 1×3 | 颜色 RGB |
| `style` | char | 线型/标记（如 `'-'`、`'o'`） |
| `measured` | logical | 是否测量点（决定抽稀与线宽） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 仅绘图 |

### `local_quantile`

**文件：** `reproduction/Observation of temporal reflection/fig01_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function value = local_quantile(x, q)
```

**简介：** 计算一维数据的分位数（线性插值）。

**功能：** 排序有限值后按位置公式线性插值得到分位数；空输入返回 NaN。用于提取逻辑高低电平。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `x` | double 向量 | 输入数据 |
| `q` | double 标量 | 分位点（0–1） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `value` | double 标量 | 分位数值 |

### `local_find_data_dir`

**文件：** `reproduction/Observation of temporal reflection/fig01_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function dataDir = local_find_data_dir(scriptDir)
```

**简介：** 定位官方源数据目录 `data/Time_Interface_Source_Data`。

**功能：** 依次在脚本目录、当前目录下的候选路径与 `genpath` 搜索结果中查找叶子名为 `Time_Interface_Source_Data` 的目录，找不到则抛 `otr:DataDirectoryMissing`。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `scriptDir` | char | 脚本所在目录 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `dataDir` | char | 找到的数据目录绝对路径 |

### `local_set_fonts`

**文件：** `reproduction/Observation of temporal reflection/fig01_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function local_set_fonts(fig)
```

**简介：** 统一设置图中所有坐标轴的字体、字号与线宽。

**功能：** 遍历 `fig` 中所有 `Axes`，将 `FontName` 设为 Helvetica、`FontSize` 设为 9、`LineWidth` 设为 0.8。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `fig` | matlab.ui.Figure | 目标图窗 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 仅修改图形属性 |

### `fig02_source_data`

**文件：** `reproduction/Observation of temporal reflection/fig02_source_data.m`

**类型：** 函数

**函数签名：**

```matlab
function result = fig02_source_data()
```

**简介：** 用官方源数据复现论文 Fig. 2c–f（红移/蓝移波包与频率扫描）。

**功能：**

读取红移/蓝移示波器 CSV 与频率扫描 CSV，用显式时间门（`gates`）在原始示波器时间坐标上把 CH1 入射、CH1 反射、CH2 折射三个分量切分出来，经 `otr_gate_signal` 做线性端点基线去除与 10% 余弦锥度（taper）后分别做补零 FFT，用对数抛物线插值求谱峰频率。对频率扫描测量列与作者提供的理论列计算 RMSE。绘制时间波形、归一化频谱及宽带红移/蓝移扫描对比图，输出 PNG 与 MAT。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无输入参数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `result` | struct | 含 `description`、`dataDirectory`、`files`、`gates`（各分量时间门）、`taperFraction`、`redshift`、`blueshift`（各含三分量时间/频谱/峰值）、`redshiftSweep`/`blueshiftSweep`（数据与 RMSE）、`outputFiles`；写入 `fig02_source_data.png` 与 `fig02_source_data.mat`。 |

**备注：** 依赖 `otr_read_scope_csv`、`otr_read_numeric_table`、`otr_gate_signal`。FFT 使用 `nextpow2` 补零到至少 16384 点；谱峰在 5–100 MHz 内搜索并做对数抛物线插值。

### `local_analyse_scope`

**文件：** `reproduction/Observation of temporal reflection/fig02_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function analysed = local_analyse_scope(scope, gateLimits, taperFraction)
```

**简介：** 对一个示波器记录做门控、去趋势与频谱分析。

**功能：** 用 `otr_gate_signal` 对 CH1（入射、反射）与 CH2（折射）加时间门，再对各分量做频谱分析，返回包含三分量时域/频域信息与峰值频率（MHz）的结构体。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `scope` | struct | `otr_read_scope_csv` 输出，含 `time`、`channels` |
| `gateLimits` | struct | 字段 `incident`、`reflected`、`transmitted` 各为 1×2 时间门（秒） |
| `taperFraction` | double 标量 | 锥度比例（如 0.10） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `analysed` | struct | `scope`、`incident`、`reflected`、`transmitted`、`peakMHz` |

### `local_spectrum`

**文件：** `reproduction/Observation of temporal reflection/fig02_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function [frequency, amplitude, peakFrequency] = local_spectrum(gate)
```

**简介：** 对门控信号做单边 FFT 并返回谱峰频率。

**功能：** 对 `gate.signal` 补零到 2 的幂（至少 16384），取单边谱并按 `dt` 归一化幅度，在 5–100 MHz 内寻找最大幅值并做对数抛物线插值得到峰值频率。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `gate` | struct | 含 `signal`、`dt` 字段的门控信号 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `frequency` | double 向量 | 单边谱频率轴（Hz） |
| `amplitude` | double 向量 | 单边谱幅度 |
| `peakFrequency` | double 标量 | 插值后的谱峰频率（Hz） |

### `local_parabolic_peak`

**文件：** `reproduction/Observation of temporal reflection/fig02_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function peak = local_parabolic_peak(frequency, amplitude, index)
```

**简介：** 在幅度对数域做三点抛物线插值求谱峰。

**功能：** 对以 `index` 为中心的三个相邻谱点取对数幅度做抛物线顶点插值，偏移量限制在 ±1 个频率步内，返回更精确的峰值频率。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `frequency` | double 向量 | 频率轴 |
| `amplitude` | double 向量 | 幅度谱 |
| `index` | double 标量 | 最大幅值索引 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `peak` | double 标量 | 插值峰值频率 |

### `local_plot_time`

**文件：** `reproduction/Observation of temporal reflection/fig02_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function local_plot_time(ax, analysed, cIncident, cReflected, cTransmitted)
```

**简介：** 绘制入射、时间反射、折射三个分量的去趋势时域波形。

**功能：** 将 `analysed` 中三个门控分量的 `detrended` 波形按给定颜色绘制到 `ax` 上，标注坐标轴并加图例。无返回值。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `ax` | Axes | 目标坐标轴 |
| `analysed` | struct | `local_analyse_scope` 输出 |
| `cIncident`/`cReflected`/`cTransmitted` | double 1×3 | 各分量颜色 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 仅绘图 |

### `local_plot_spectra`

**文件：** `reproduction/Observation of temporal reflection/fig02_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function local_plot_spectra(ax, analysed, cIncident, cReflected, cTransmitted)
```

**简介：** 绘制三个分量的归一化频谱。

**功能：** 调用 `local_spectrum_line` 绘制入射、反射、折射分量的归一化单边谱，限制在 0–100 MHz、0–1.05 范围内。无返回值。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `ax` | Axes | 目标坐标轴 |
| `analysed` | struct | `local_analyse_scope` 输出 |
| `cIncident`/`cReflected`/`cTransmitted` | double 1×3 | 各分量颜色 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 仅绘图 |

### `local_spectrum_line`

**文件：** `reproduction/Observation of temporal reflection/fig02_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function local_spectrum_line(ax, gate, color, name)
```

**简介：** 绘制单条归一化频谱曲线。

**功能：** 选取 0–100 MHz 内的频谱，按该频带内最大值归一化后绘制。无返回值。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `ax` | Axes | 目标坐标轴 |
| `gate` | struct | 含 `frequencyHz`、`spectrum` 字段 |
| `color` | double 1×3 | 颜色 |
| `name` | char | 图例名 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 仅绘图 |

### `local_rmse`

**文件：** `reproduction/Observation of temporal reflection/fig02_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function value = local_rmse(measured, theory)
```

**简介：** 计算两列数据的均方根误差。

**功能：** 仅对两者均有限的采样点计算 $\sqrt{\frac1N\sum(\text{measured}-\text{theory})^2}$，用于比较测量与理论频率扫描列。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `measured` | double 向量 | 测量值 |
| `theory` | double 向量 | 理论值 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `value` | double 标量 | RMSE |

### `local_find_data_dir`

**文件：** `reproduction/Observation of temporal reflection/fig02_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function dataDir = local_find_data_dir(scriptDir)
```

**简介：** 定位官方源数据目录 `data/Time_Interface_Source_Data`。

**功能：** 同 fig01 中的同名函数：在脚本目录与当前目录的候选路径及 `genpath` 结果中查找目录，找不到则抛 `otr:DataDirectoryMissing`。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `scriptDir` | char | 脚本所在目录 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `dataDir` | char | 数据目录绝对路径 |

### `local_set_fonts`

**文件：** `reproduction/Observation of temporal reflection/fig02_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function local_set_fonts(fig)
```

**简介：** 统一设置图中所有坐标轴的字体、字号与线宽。

**功能：** 将图中所有 `Axes` 设为 Helvetica、9 号字、线宽 0.8。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `fig` | matlab.ui.Figure | 目标图窗 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 仅修改图形属性 |

### `fig03_source_data`

**文件：** `reproduction/Observation of temporal reflection/fig03_source_data.m`

**类型：** 函数

**函数签名：**

```matlab
function result = fig03_source_data()
```

**简介：** 绘制官方 Fig. 3c/d 的原始示波器波形（仅做稳健偏移去除）。

**功能：**

读取 `tau_15/25/35.csv` 三个时间板记录，在原生示波器时间坐标上显示 CH1、CH2 与 CH3 开关监视逻辑。仅对 CH1/CH2 做基于静默段的稳健中位数偏移去除、对 CH3 用 5%/95% 分位数归一化；不推断 Fig. 3d 的谱门（因源数据未给出唯一门或对照采集），输出明确标注为 RAW DATA。同时计算逻辑边沿（上升/下降/高电平持续）与 CH1/CH2 峰峰值摘要，输出 PNG 与 MAT。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无输入参数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `result` | struct | 含 `description`、`limitation`、`dataDirectory`、`files`、`nominalTauNs`、`scope`、`rawDisplayTraces`、`edgeSummary`、`displayWindowNs`、`outputFiles`；写入 `fig03_source_data.png` 与 `fig03_source_data.mat`。 |

**备注：** 依赖 `otr_read_scope_csv`。该脚本明确声明不声称复现 Fig. 3d 的谱（原始数据未提供唯一门与泄漏参考段）。

### `local_raw_trace`

**文件：** `reproduction/Observation of temporal reflection/fig03_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function trace = local_raw_trace(scope, limits)
```

**简介：** 在显示窗口内对示波器记录做稳健偏移去除与逻辑归一化。

**功能：** 截取 `limits` 时间窗内三通道均有限的数据；用窗口前段静默区的中位数估计 CH1/CH2 直流偏移并减去，用 5%/95% 分位数归一化 CH3 逻辑信号；静默样本不足时回退到前 10% 采样。返回含原始与处理后的通道及偏移量的结构体。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `scope` | struct | 含 `time`、`channels` 的示波器记录 |
| `limits` | double 1×2 | 显示时间窗 `[t_low, t_high]`（秒） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `trace` | struct | `time`、`ch1Raw`/`ch2Raw`/`logicRaw`、`ch1OffsetV`/`ch2OffsetV`、`ch1`/`ch2`（去偏移）、`logicLevelsV`、`logicNormalized` |

### `local_edge_summary`

**文件：** `reproduction/Observation of temporal reflection/fig03_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function summary = local_edge_summary(trace)
```

**简介：** 从归一化逻辑信号提取上升/下降/高电平持续与通道峰峰值。

**功能：** 在归一化逻辑信号中定位首次穿过 0.5 的上升沿与随后首个下降沿，用线性插值求精确穿越时刻，计算高电平持续时间，并统计 CH1/CH2 的峰峰值。无有效边沿时返回 NaN。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `trace` | struct | `local_raw_trace` 输出 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `summary` | struct | `riseTimeNs`、`fallTimeNs`、`highDurationNs`、`ch1PeakToPeakV`、`ch2PeakToPeakV` |

### `local_crossing`

**文件：** `reproduction/Observation of temporal reflection/fig03_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function tCross = local_crossing(time, value, rightIndex, level)
```

**简介：** 用线性插值求信号穿过给定电平的时刻。

**功能：** 在 `rightIndex` 与其左邻采样之间按 `level` 做线性插值，返回精确穿越时刻，插值比例截断到 [0,1]。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `time` | double 向量 | 时间轴 |
| `value` | double 向量 | 信号值 |
| `rightIndex` | double 标量 | 穿越点右侧索引 |
| `level` | double 标量 | 目标电平 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `tCross` | double 标量 | 穿越时刻 |

### `local_quantile`

**文件：** `reproduction/Observation of temporal reflection/fig03_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function value = local_quantile(x, q)
```

**简介：** 计算一维数据的分位数（线性插值）。

**功能：** 对排序后的有限值按位置公式线性插值，用于估计逻辑高/低电平。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `x` | double 向量 | 输入数据 |
| `q` | double 标量 | 分位点（0–1） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `value` | double 标量 | 分位数值 |

### `local_find_data_dir`

**文件：** `reproduction/Observation of temporal reflection/fig03_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function dataDir = local_find_data_dir(scriptDir)
```

**简介：** 定位官方源数据目录 `data/Time_Interface_Source_Data`。

**功能：** 同 fig01 的同名函数，查找数据目录，找不到抛 `otr:DataDirectoryMissing`。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `scriptDir` | char | 脚本所在目录 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `dataDir` | char | 数据目录绝对路径 |

### `local_set_fonts`

**文件：** `reproduction/Observation of temporal reflection/fig03_source_data.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function local_set_fonts(fig)
```

**简介：** 统一设置图中所有坐标轴的字体、字号与线宽。

**功能：** 将图中所有 `Axes` 设为 Helvetica、9 号字、线宽 0.8。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `fig` | matlab.ui.Figure | 目标图窗 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 仅修改图形属性 |

### `fig04_boundary_conditions`

**文件：** `reproduction/Observation of temporal reflection/fig04_boundary_conditions.m`

**类型：** 函数

**函数签名：**

```matlab
function results = fig04_boundary_conditions()
```

**简介：** 复现 SI Fig. S3 与式 (S13)–(S17) 的微观时域边界条件。

**功能：**

构造两个微观示例：左列是“并联电容接入”（电荷守恒，$C_2/C_1=4$，频率减半），右列是“并联电容移除”（电压守恒，$C_2/C_1=1/4$，频率加倍）。电压波形在 $t\ge0$ 后按 $v=T\cos(\omega_2 t)+R\cos(-\omega_2 t)$ 用 `local_coeff` 得到的透射/反射系数合成，并分别追踪基电容电荷、新增电容电荷与总电荷（以及移除支路电荷与线路总电荷）。绘制四幅图并做电荷/电压连续性自检（容差 5e-3）。输出 PNG/FIG 与 MAT。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无输入参数，ON/OFF 电容比在函数体内硬编码。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `results` | struct | `t`（归一化时间）、`on`、`off`（各含电压、电荷分量、T、R、w2 等）；写入 `fig04_boundary_conditions.mat`。 |

**备注：** 依赖 `otr_setup`、`otr_save_figure`。时间以切换前周期归一化；电荷连续（接入）与电压连续（移除）是两种不同的微观边界定律，脚本意在使这一区别可见而非假定统一 D/B 跳变。

### `local_coeff`

**文件：** `reproduction/Observation of temporal reflection/fig04_boundary_conditions.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function [M, T, R] = local_coeff(rootRatio, law)
```

**简介：** 按电荷连续或电压连续定律计算时域界面透射/反射系数。

**功能：** `law='charge'`（电荷连续）时 $T=\tfrac12(r^2+r)$、$R=\tfrac12(r^2-r)$；`law='voltage'`（电压连续）时 $T=\tfrac12(1+r)$、$R=\tfrac12(1-r)$，其中 $r$ 为 `rootRatio`。同时返回系数矩阵 $M=[T\ R;\,R\ T]$。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `rootRatio` | double 标量 | 电容（或阻抗）比的平方根 $r$ |
| `law` | char | `'charge'` 或 `'voltage'` |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `M` | double 2×2 | 系数矩阵 `[T R; R T]` |
| `T` | double 标量 | 透射系数 |
| `R` | double 标量 | 反射系数 |

### `fig05_ideal_time_reflection`

**文件：** `reproduction/Observation of temporal reflection/fig05_ideal_time_reflection.m`

**类型：** 函数

**函数签名：**

```matlab
function results = fig05_ideal_time_reflection()
```

**简介：** 展示理想光子时间界面的时间反演与动量（波数）守恒。

**功能：**

调用 `otr_synthetic_wavepacket` 生成空间均匀瞬时切换下的合成波包，绘制 |场| 时空图（标出切换时刻）、切换前后及前向/后向（时折射/时反射）方向分量的一维剖面，以及切换前后经 `fftshift(fft)` 得到的归一化 k 谱（验证 k 守恒并给出频率比 $f_2/f_1$）。输出 PNG/FIG 与 MAT。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无输入参数；切换方向 `on`、`tSwitch=45 ns`、时间轴在函数体内设定。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `results` | struct | 直接为 `otr_synthetic_wavepacket` 的输出（含 `x`、`t`、`field`、`forward`、`backward`、`tSwitch`、`initial`、`k`、`k0`、`fRatio`、`Tfirst`、`Rfirst` 等）；写入 `fig05_ideal_time_reflection.mat`。 |

**备注：** 依赖 `otr_setup`、`otr_synthetic_wavepacket`、`otr_save_figure`。k 谱保留窗为 `|k-k0| <= max(8, 0.8*k0)`。

### `fig06_temporal_slab_theory`

**文件：** `reproduction/Observation of temporal reflection/fig06_temporal_slab_theory.m`

**类型：** 函数

**函数签名：**

```matlab
function results = fig06_temporal_slab_theory()
```

**简介：** 复现时间板（temporal Fabry-Pérot）的反射零点与四条因果路径干涉现象。

**功能：**

用 OFF→ON（电荷连续）与 ON→OFF（电压连续）两组界面系数 $T_{\rm on},R_{\rm on},T_{\rm off},R_{\rm off}$，对三条板持续 $\tau=15/25/35$ ns 计算时间板净反射与净透射幅值 $|R|=|R_{\rm off}T_{\rm on}e^{-i\phi}+T_{\rm off}R_{\rm on}e^{i\phi}|$（$\phi=2\pi f(Z_{\rm on}/Z_{\rm off})\tau$）。绘制反射零点和透射幅值随频率的曲线，并做固定频率 $f=28.66$ MHz 下随板持续时间的连续调谐曲线。输出 PNG/FIG 与 MAT。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无输入参数；阻抗 $Z_{\rm off}=50,\ Z_{\rm on}=25$ 硬编码。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `results` | struct | `f`、`taus`、`R`、`T`、`tauSweep`、`Rtau`、`coefficients`（含 Ton/Ron/Toff/Roff）；写入 `fig06_temporal_slab_theory.mat`。 |

**备注：** 依赖 `otr_setup`、`otr_save_figure`。采用 $e^{-i\phi}$ 时间因子，四条因果路径幅度为 `[Ton*Toff, Ron*Roff, Ron*Toff, Ton*Roff]`，前两条为正、后两条为负（脚本据此做极性自检）。

### `local_interface`

**文件：** `reproduction/Observation of temporal reflection/fig06_temporal_slab_theory.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function [M, T, R] = local_interface(rootRatio, law)
```

**简介：** 按电荷连续或电压连续定律计算时域界面透射/反射系数。

**功能：** 与 fig04 的 `local_coeff` 相同：`charge` 时 $T=\tfrac12(r^2+r)$、$R=\tfrac12(r^2-r)$；`voltage` 时 $T=\tfrac12(1+r)$、$R=\tfrac12(1-r)$；返回矩阵 $M=[T\ R;\,R\ T]$。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `rootRatio` | double 标量 | 阻抗比 $Z_2/Z_1$（或其平方根含义的 $r$） |
| `law` | char | `'charge'` 或 `'voltage'` |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `M` | double 2×2 | 系数矩阵 `[T R; R T]` |
| `T` | double 标量 | 透射系数 |
| `R` | double 标量 | 反射系数 |

### `fig07_finite_switching_and_homogeneity`

**文件：** `reproduction/Observation of temporal reflection/fig07_finite_switching_and_homogeneity.m`

**类型：** 函数

**函数签名：**

```matlab
function results = fig07_finite_switching_and_homogeneity()
```

**简介：** 复现 SI Fig. S6/S7 的有限上升时间与空间同步性趋势。

**功能：**

左侧用一阶 Fourier 重叠估计有限线性频率斜坡对时域反射峰值的抑制 $|\mathrm{sinc}(\Delta\omega\,\tau/2)|$，其中 $\Delta\omega=2\pi f_0(1-\mathrm{ratio})$、$f_0=50$ MHz、$\mathrm{ratio}=0.55$，并把论文 3 ns 约 90% 的测量值作为独立锚点标注。右侧计算 30 节阵列的同步相干 $|\frac1N\sum_n e^{-i2\pi f\,n\,\Delta t}|$，其中控制线延迟设为信号线单元延迟的 1/80。绘制两幅图并输出 PNG/FIG 与 MAT。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无输入参数；$f_0$、ratio、N=30 等在函数体内硬编码。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `results` | struct | `riseNs`、`finiteFactor`、`publishedAnchor`、`frequency`、`coherence`、`controlTotalDelay`；写入 `fig07_finite_switching_and_homogeneity.mat`。 |

**备注：** 依赖 `otr_setup`、`otr_save_figure`。有限上升曲线是显式标注的解析切换时间研究，非未发表的 ADS 数据。

### `fig08_inverted_slab_and_leakage`

**文件：** `reproduction/Observation of temporal reflection/fig08_inverted_slab_and_leakage.m`

**类型：** 函数

**函数签名：**

```matlab
function results = fig08_inverted_slab_and_leakage()
```

**简介：** 生成 Extended Data Fig. 8/9（反相板损耗与视频泄漏抵消）的方程驱动预测。

**功能：**

无公开源数据，故由方程直接生成预测：对反相 ON-OFF-ON 板，在反射幅值中加入随频率的导体损耗（$0.5\log(10)/20$ Np/m，对应 100 MHz 处 0.5 dB/m）与介电损耗项，得到损耗极小值不归零的 |R| 曲线；并构造含指数衰减视频泄漏的仿真信号，演示 $V(t)-V(t-T)$ 泄漏抵消。绘制两幅图并输出 PNG/FIG 与 MAT。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无输入参数；$Z_{\rm off}=50,\ Z_{\rm on}=25$、板持续 $\tau=17/21/25$ ns 等硬编码。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `results` | struct | `f`、`taus`、`R`、`t`、`raw`、`compensated`；写入 `fig08_inverted_slab_and_leakage.mat`。 |

**备注：** 依赖 `otr_setup`、`otr_save_figure`。损耗模型含 $\sqrt{f/100\ \mathrm{MHz}}$ 的趋肤项与介电损耗项；泄漏抵消为时域差分 $V(t)-V(t-T)$。

### `local_interface`

**文件：** `reproduction/Observation of temporal reflection/fig08_inverted_slab_and_leakage.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function [M, T, R] = local_interface(rootRatio, law)
```

**简介：** 按电荷连续或电压连续定律计算时域界面透射/反射系数。

**功能：** 与 fig06 的 `local_interface` 相同：`charge` 时 $T=\tfrac12(r^2+r)$、$R=\tfrac12(r^2-r)$；`voltage` 时 $T=\tfrac12(1+r)$、$R=\tfrac12(1-r)$；返回 $M=[T\ R;\,R\ T]$。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `rootRatio` | double 标量 | 阻抗比 $r$ |
| `law` | char | `'charge'` 或 `'voltage'` |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `M` | double 2×2 | 系数矩阵 `[T R; R T]` |
| `T` | double 标量 | 透射系数 |
| `R` | double 标量 | 反射系数 |

### `fig09_tlm_design`

**文件：** `reproduction/Observation of temporal reflection/fig09_tlm_design.m`

**类型：** 函数

**函数签名：**

```matlab
function results = fig09_tlm_design()
```

**简介：** 复现 SI Fig. S1 的周期加载传输线单元设计分析。

**功能：**

用 `otr_bloch_dispersion` 计算开关 OFF/ON 理想模型下的 Bloch 色散 $\beta d$ 与 Bloch 阻抗 $Z_B$，绘制色散与阻抗曲线。以 $\theta$（$\omega_0$ 处的 Bloch 相位）与截止频率比 $\omega_c/\omega_0$ 为自变量构建设计图：由 T 型电容单元在 $\omega_c$ 处 $A=-1$ 的条件反解 $b_c=2(1+\cos\theta_c)/\sin\theta_c$，在 $\omega/\omega_0=0.5$ 处计算 ON 态归一化阻抗与反射系数 $R_{\rm on}$（S16，归一化 $Z_{\rm off}=1$），并用求根得到等 $\beta$ 频率比映射。输出 PNG/FIG 与 MAT。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无输入参数；`thetaDeg`、`cutoffRatio` 网格在函数体内定义。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `results` | struct | `frequency`、`betaOff`、`betaOn`、`Zoff`、`Zon`、`thetaDeg`、`cutoffRatio`、`reflectionMap`、`frequencyRatioMap`；写入 `fig09_tlm_design.mat`。 |

**备注：** 依赖 `otr_setup`、`otr_parameters`、`otr_bloch_dispersion`、`otr_save_figure`。设计图对应 SI 式 S2/S3 与 S16。

### `fig10_scattering_spectra`

**文件：** `reproduction/Observation of temporal reflection/fig10_scattering_spectra.m`

**类型：** 函数

**函数签名：**

```matlab
function results = fig10_scattering_spectra()
```

**简介：** 复现主图 Fig. 2a/b 的时域散射谱理论内容。

**功能：**

无公开处理后数据，曲线为 SI S1–S3 与两种微观定律 S16/S17 的预测。用 `otr_frequency_map` 把 OFF 态频率映射到 ON 态（固定 k），用 `otr_bloch_dispersion` 求两组 Bloch 阻抗，再分别用 `otr_temporal_interface(...,'charge')` 与 `'voltage'` 计算 OFF→ON（电荷连续）与 ON→OFF（电压连续）的透射/反射系数，绘制幅度谱与相位谱（$|T|,|R|,\angle T,\angle R$ 随 k）。输出 PNG/FIG 与 MAT。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无输入参数。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `results` | struct | `fOff`、`fOn`、`k`、`Zoff`、`Zon`、`Ton`、`Ron`、`Toff`、`Roff`；写入 `fig10_scattering_spectra.mat`。 |

**备注：** 依赖 `otr_setup`、`otr_parameters`、`otr_frequency_map`、`otr_bloch_dispersion`、`otr_temporal_interface`、`otr_save_figure`。常数阻抗目标（$Z_{\rm on}/Z_0=0.5$）对应 $T=0.375,\ R=-0.125$（电荷连续）与 $T=1.5,\ R=-0.5$（电压连续）。

### `fig11_slab_wavepackets`

**文件：** `reproduction/Observation of temporal reflection/fig11_slab_wavepackets.m`

**类型：** 函数

**函数签名：**

```matlab
function results = fig11_slab_wavepackets()
```

**简介：** 复现 SI Fig. S5 与主图 Fig. 3e 的时间板波包预测。

**功能：**

以有效板持续 $\tau_{\rm eff}=1/(4\times21.2\ \mathrm{MHz})$ 使 38 MHz 附近映射到板的首个零点，用 $R=a\,e^{-i2\pi f_{\rm on}\tau_{\rm eff}}+a\,e^{+i2\pi f_{\rm on}\tau_{\rm eff}}$（$a=-0.1875$，对应 $R_{\rm off}T_{\rm on}=T_{\rm off}R_{\rm on}$）构造时间板反射传递函数；对 17.5 MHz-FWHM 高斯包络（中心 30/35/38/41/46 MHz）计算反射波形与频谱，绘制波包与共享 38 MHz 反射零点谱。另做脉冲带宽平均的板持续扫描（Fig. 3e 预测）。输出两张图与 MAT。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无输入参数；`centres`、`tauEff`、`fwhm` 等在函数体内定义。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `results` | struct | `f`、`centres`、`tauEffective`、`R`、`spectra`、`t`、`waveforms`、`tau`、`widths`、`integratedReflection`；写入 `fig11_slab_wavepackets.mat`。 |

**备注：** 依赖 `otr_setup`、`otr_parameters`、`otr_frequency_map`、`otr_save_figure`。反射谱为 $R=2a\cos(2\pi f_{\rm on}\tau_{\rm eff})$ 的对称组合；带宽平均反射为 $\sqrt{\frac{\int|R(f)\cdot S(f)|^2 df}{\int|S(f)|^2 df}}$。

### `run_all_reproductions`

**文件：** `reproduction/Observation of temporal reflection/run_all_reproductions.m`

**类型：** 函数

**函数签名：**

```matlab
function summary = run_all_reproductions(varargin)
```

**简介：** 顺序执行本目录下全部数值复现任务并生成运行清单。

**功能：**

解析可选参数（`'quick'` 模式或名称/值对），按 `local_jobs` 构建任务列表（官方数据图 fig01–03、解析/电路图 fig04–11、MNA 演示 `demo_circuit_time_interface` 与独立验证 `otr_validate`），逐个调用并记录状态、现象与耗时，逐任务捕获异常（按 `ContinueOnError` 决定是否中止），并关闭新建图窗。最后写 `run_manifest.mat` 清单，若有失败则抛 `otr:ReproductionFailed`。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `varargin` | 变长 | 可选。`'quick'`（跳过官方大数据与 MNA 场图）、`'all'`，或名称/值对：`IncludeSourceData`（logical，默认 true）、`IncludeCircuit`（logical，默认 true）、`RunValidation`（logical，默认 true）、`CloseGeneratedFigures`（logical，默认 true）、`ContinueOnError`（logical，默认 false）。 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `summary` | struct | 运行清单，含 `generatedAt`、`matlabVersion`、`root`、`options`、`entries`（各任务名/现象/状态/耗时/消息）、`elapsedSeconds`、`allPassed`；写入 `run_manifest.mat`。 |

**备注：** 依赖 `otr_setup` 及各 fig 脚本与 `otr_validate`。各任务以函数句柄方式组织，图窗通过对比 `groot` 中新建 figure 来关闭。

### `local_jobs`

**文件：** `reproduction/Observation of temporal reflection/run_all_reproductions.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function jobs = local_jobs(options)
```

**简介：** 根据选项构建复现任务列表。

**功能：** 按 `IncludeSourceData` 添加 fig01–03，随后固定添加 fig04–11，再按 `IncludeCircuit`/`RunValidation` 添加 MNA 演示与 `otr_validate`，返回含 `name`/`phenomenon`/`call` 的任务结构体数组。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `options` | struct | 由 `local_options` 生成的选项 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `jobs` | struct 数组 | 每项含 `name`、`phenomenon`、`call`（函数句柄） |

### `local_job`

**文件：** `reproduction/Observation of temporal reflection/run_all_reproductions.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function job = local_job(name, phenomenon, call)
```

**简介：** 构造单个任务结构体。

**功能：** 将名称、现象描述与函数句柄打包成 struct 返回。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `name` | char | 任务名 |
| `phenomenon` | char | 现象描述 |
| `call` | function_handle | 待调用的函数句柄 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `job` | struct | `name`、`phenomenon`、`call` |

### `local_save_manifest`

**文件：** `reproduction/Observation of temporal reflection/run_all_reproductions.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function summary = local_save_manifest(paths, options, entries, elapsedSeconds)
```

**简介：** 保存运行清单到 `run_manifest.mat`。

**功能：** 汇总生成时间、MATLAB 版本、复现根目录、选项、各任务条目与总耗时，计算 `allPassed`，保存并返回 `summary`。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `paths` | struct | `otr_setup` 输出，含 `output`、`reproduction` 路径 |
| `options` | struct | 选项结构体 |
| `entries` | struct 数组 | 各任务状态记录 |
| `elapsedSeconds` | double 标量 | 总耗时（秒） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `summary` | struct | 见 `run_all_reproductions` 输出说明 |

### `local_options`

**文件：** `reproduction/Observation of temporal reflection/run_all_reproductions.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function options = local_options(varargin)
```

**简介：** 解析 `run_all_reproductions` 的可选参数。

**功能：** 初始化默认选项；识别 `'quick'`（关闭源数据与电路）与 `'all'` 模式；否则要求名称/值对，校验选项名存在且值为逻辑标量，返回选项结构体。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `varargin` | 变长 | `'quick'`/`'all'` 或名称/值对 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `options` | struct | 含 `IncludeSourceData`、`IncludeCircuit`、`RunValidation`、`CloseGeneratedFigures`、`ContinueOnError` |

### `local_onoff`

**文件：** `reproduction/Observation of temporal reflection/run_all_reproductions.m`

**类型：** 函数（局部辅助）

**函数签名：**

```matlab
function label = local_onoff(value)
```

**简介：** 将逻辑值转为 'on'/'off' 字符串。

**功能：** 用于打印选项状态：真返回 `'on'`，假返回 `'off'`。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `value` | logical | 逻辑值 |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `label` | char | `'on'` 或 `'off'` |
## reproduction/Topology of photonic crystals — 复现脚本

### `compute_zak_phases`

**文件：** `reproduction/Topology of photonic crystals/compute_zak_phases.m`

**类型：** 函数

**函数签名：**

```matlab
function [zakPhases, info] = compute_zak_phases(varargin)
```

**简介：** 计算二值光子时间晶体 (PTC) 前若干能带的 Zak 相位。

**功能：**

实现对 Lustig, Sharabi 和 Segev, *Optica* 5, 1390 (2018) 中 Eq. (5) 的离散化——闭合 Floquet 频率 Wilson 回路，即 θ_m = ∫_BZ dΩ i⟨u_{m,Ω}|d_Ω u_{m,Ω}⟩，其中内积为 ⟨u|v⟩ = ∫_0^T ε(t) u*(t)v(t) dt。对每个 Ω（在半开时间布里渊区 [-π/T, π/T) 上取中点采样），先由精确二值晶体色散关系反解守恒动量 k_m(Ω)，再通过中心化单值矩阵 `centred_monodromy` 的本征矢重构整个时间周期单元内的位移场 D，并用 ε 加权内积归一化。相邻 Ω 的 Wilson 链由该内积计算，最后一根链使用时间 Bloch 缝合变换 u_{Ω+2π/T}(t) = exp(i·2π·t/T)·u_Ω(t) 闭合回路，保证规范不变。反演对称使 Zak 相位量子化为 0 或 π，函数存储原始 Wilson 相位并在 `real(wilson) ≥ 0` 时量化为 0，否则为 π。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `InversionCenter` | char/string（可选） | 时间反演中心所在的介电常数段，`'low'`（默认，对应论文 Fig. 1 标记）或 `'high'` |
| `NumberOfBands` | double（可选） | 要计算的完整能带数（默认 7） |
| `NOmega` | double（可选） | 半开时间 BZ 上的采样点数（默认 61，最小 15） |
| `NTime` | double（可选） | 一个时间周期内的采样数（默认 301，最小 81） |
| `KMax` | double（可选） | 定位能带所用的最大 k/k0（默认 5） |
| `Nk` | double（可选） | 能带搜索的采样点数（默认 20001，最小 1001） |
| `ComputeBothOrigins` | logical（可选） | 是否独立计算两个反演中心（默认 true） |
| `MakePlots` | logical（可选） | 是否绘制诊断图（默认 `nargout==0`） |
| `SaveOutput` | logical（可选） | 是否保存诊断结果到 output/（默认 `nargout==0`） |
| `Verbose` | logical（可选） | 是否打印数值审计表（默认 true） |

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `zakPhases` | double（1×nBands） | 量子化后的 Zak 相位（0 或 π） |
| `info` | struct | 详细结果：原始 Wilson 相位、能带/带隙边界、Ω 网格、k(Ω)、链重叠诊断、反演中心约定、两个中心的 Zak 序列等 |

**备注：** 依赖本文件内的局部函数 `binary_half_trace`、`logical_segments`、`refine_band_edge`、`solve_band_momentum`、`centred_monodromy`、`segment_propagator`、`centred_displacement_mode`、`advance_state`、`zak_label`、`make_diagnostic_figure`。物理约定：c=1（仓库 TMM 约定）、T=2 fs、k0=2π/T、ε_high=3、ε_low=1；中心化对称单元为 [ε_centre(T/4), ε_other(T/2), ε_centre(T/4)]，采用 D/B 连续约定。Zak 相位依赖时间原点，原点移动 T/2 会使每个能带相位整体加 π。

### `fig1_ptc_bands`

**文件：** `reproduction/Topology of photonic crystals/fig1_ptc_bands.m`

**类型：** 函数

**函数签名：**

```matlab
function fig1_ptc_bands()
```

**简介：** 复现 Lustig et al., *Optica* 5, 1390 (2018) 的 Fig. 1（含 (a) 二元 PTC 示意图与 (b) 扩展布里渊区色散及 Zak 标记）。

**功能：**

用 `temporal_crystal_bands` 计算 k/k0 ∈ [0,5] 的能带结构，通过半迹 `|halfTrace| ≤ 1+5e-11` 判定带内/带隙并分段。在带内按 cos(ΩT)=Tr(U)/2 取 acos 得到扩展区间的 ±ΩT 两支，逐段绘制以避免跨带隙的虚假连线。调用 `compute_zak_phases`（低 ε 中心、7 个能带）计算 Zak 标记，并与论文参考序列 [0 0 π 0 0 π 0]（低 ε 中心）及其高 ε 中心变体（整体加 π）逐一比对，作为验收测试而非绘图数据。左面板 (a) 绘制以整数周期为中心的 ε(t) 台阶（时间反演中心位于低 ε 段），右面板 (b) 绘制灰色带隙阴影、蓝色能带与 Zak 文本标签。

**输入参数：**

无

**输出参数：**

无（仅绘图/写文件）。保存 `output/fig1_ptc_bands.png` 与 `output/fig1_data.mat`（含 `kNorm`、`omegaT`、`halfTrace`、`zakPhases`、`zakInfo` 等）。

**备注：** 依赖 `temporal_crystal_bands`、`compute_zak_phases` 及本文件局部函数 `logical_segments`；运行父目录的 `startup_stm.m`。采用 c=1 TMM 坐标，k0=2π/T，ε_1=3、ε_2=1、T=2。明确指出论文中 ε_1/ε_2 标记在面板 (a) 与正文间不一致。

### `fig1a_ptc_schematic`

**文件：** `reproduction/Topology of photonic crystals/fig1a_ptc_schematic.m`

**类型：** 函数

**函数签名：**

```matlab
function fig1a_ptc_schematic()
```

**简介：** 复现 Lustig et al. (arXiv:1803.08731) 的 Fig. 1(a)——二值光子时间晶体 ε(t) 示意图。

**功能：**

构建交替取值 ε_1=3、ε_2=1、等持续时间 t1=t2=T/2 的介电常数时间分布，其中时间原点 t=0 位于第 1 段中点，以保持时间反演对称。对多个周期（4 个）在精细采样上按 `mod(t+T/2, T)` 折叠到 [0,T) 区间构造 ε(t)，用填充色块与 `stairs` 绘制，并以红色虚线标记 t=0（TR 中心）。标注 ε_1、ε_2 段、周期 T 与 Ω=2π/T。

**输入参数：**

无

**输出参数：**

无（仅绘图/写文件）。保存 `output/fig1a_ptc_schematic.png`。

**备注：** 依赖父目录的 `startup_stm.m`。此处采用 T=2π（Ω=1）的归一化；段时长各为 T/2。

### `fig1b_band_structure`

**文件：** `reproduction/Topology of photonic crystals/fig1b_band_structure.m`

**类型：** 函数

**函数签名：**

```matlab
function fig1b_band_structure()
```

**简介：** 复现 Lustig et al. (arXiv:1803.08731) 的 Fig. 1(b)——二元 PTC 的 Floquet 能带结构并标注各能带 Zak 相位。

**功能：**

用 `temporal_crystal_bands` 在 k ∈ (0, 2.5] 上计算复 Floquet 频率 ω_F 与半迹，判定带内/带隙并分段。对每个能带，用 `temporal_crystal_monodromy` 的单值矩阵左右本征矢（左本征矢由 `computeLeftEigenvector` 经 U^T 求得）构造双正交 Wilson 回路，得 Zak 相位；带内分支通过与前一点本征态的最大重叠来连续追踪。绘图时以灰色阴影标记动量带隙，蓝色绘制两支实频带，红字标注量化 Zak（0 或 π）。

**输入参数：**

无

**输出参数：**

无（仅绘图/写文件）。保存 `output/fig1b_band_structure.png` 与 `output/fig1b_data.mat`（含 `kGrid`、`omegaF`、`halfTr`、`zakPhases`、`bandStarts`、`gapStarts` 等）。

**备注：** 依赖 `temporal_crystal_bands`、`temporal_crystal_monodromy` 及局部函数 `computeLeftEigenvector`、`iif`。采用 T=2π（Ω=1）、ε_1=3、ε_2=1、μ=1。Zak 相位用 `-angle(prod(links))` 后折叠到 [-π, π)。

### `fig1c_fdtd_in_band`

**文件：** `reproduction/Topology of photonic crystals/fig1c_fdtd_in_band.m`

**类型：** 函数

**函数签名：**

```matlab
function fig1c_fdtd_in_band()
```

**简介：** 兼容入口，指向已发表论文的 Fig. 2(a)（带内传播）。

**功能：**

仅发出重命名警告 `fig1c_fdtd_in_band:renamed`，然后调用 `fig2_fdtd_simulations(false)`。旧版独立脚本曾重复实现了一套错误的归一化单位 FDTD 设置，保留此兼容入口是为了防止旧的 arXiv 图名静默再生已废弃的结果。

**输入参数：**

无

**输出参数：**

无（委托 `fig2_fdtd_simulations` 绘图/写文件）。

**备注：** 依赖 `fig2_fdtd_simulations`。发表的面板实际上是 Fig. 2(a)。

### `fig1d_fdtd_in_gap`

**文件：** `reproduction/Topology of photonic crystals/fig1d_fdtd_in_gap.m`

**类型：** 函数

**函数签名：**

```matlab
function fig1d_fdtd_in_gap()
```

**简介：** 兼容入口，指向已发表论文的 Fig. 2(b)（带隙放大）。

**功能：**

仅发出重命名警告 `fig1d_fdtd_in_gap:renamed`，然后调用 `fig2_fdtd_simulations(false)`。与 `fig1c_fdtd_in_band` 同理，保留此入口以防止旧的 arXiv 图名静默再生废弃的归一化单位 FDTD 结果。

**输入参数：**

无

**输出参数：**

无（委托 `fig2_fdtd_simulations` 绘图/写文件）。

**备注：** 依赖 `fig2_fdtd_simulations`。发表的面板实际上是 Fig. 2(b)。

### `fig2_fdtd_simulations`

**文件：** `reproduction/Topology of photonic crystals/fig2_fdtd_simulations.m`

**类型：** 函数

**函数签名：**

```matlab
function fig2_fdtd_simulations(forceRecompute)
```

**简介：** 复现 Lustig et al. (2018) 的 Fig. 2——用精确 k 空间 Maxwell 传播（谱方法，替代 FDTD）展示带内传播与带隙指数放大。

**功能：**

介质空间均匀，故每个空间傅里叶分量独立演化。本实现直接在 k 空间传播精确 Maxwell 态 [D;B]，再用 FFT 重构 D(z,t)，是论文 FDTD 计算的零空间离散误差谱对应物。对两个工况分别运行：带内（λ=1.4 μm）与带隙（λ=0.93 μm，预期放大）。源谱先用 `carrier_connected_zone` 定位包含载波的连通能带/带隙，将高斯脉冲的微小拖尾在该区外置零（避免 60 周期后带外 1e-6 的尾巴盖过带内脉冲）。时间界面处用 `temporal_interface_matrix_jump` 显式拆分折射/时间反射通道（D/B 连续），并构建前若干界面的相干路径树 `build_path_tree` 作诊断。最终以对数幅度 `log10(|D|/max|D(t=0)|)` 成像，并用 `temporal_finite_crystal_response` 的 TMM 结果审计增益（指出论文 60 周期与 20000 倍增益二者自相矛盾）。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `forceRecompute` | logical（可选） | 是否强制重算两个精确 k 工况（默认 false，优先读取有效缓存） |

**输出参数：**

无（仅绘图/写文件）。保存 `output/fig2_fdtd_simulations.png`、`output/fig2_temporal_path_tree.png`、`output/fig2_audit.mat`，以及逐工况缓存 `output/fig2_inband_data.mat`、`output/fig2_ingap_data.mat`（含 `data` 结构体，字段 `D`、`logAmplitude`、`pathTree`、`interfaceAudit` 等）。

**备注：** 依赖 `run_or_load_case`、`propagate_spectrum`、`ptc_epsilon_scalar`、`build_path_tree`、`plot_path_tree`、`plot_case`、`carrier_connected_zone`、`floquet_group_speed`、`window_peaks`（均为本文件局部函数）以及 `temporal_db_to_directional`、`temporal_interface_matrix_jump`、`temporal_finite_crystal_response`、`temporal_crystal_bands`。关键物理参数：c0=0.299792458 μm/fs、ε_1=3、ε_2=1、T=2 fs、tStart=220 fs、tEnd=340 fs（共 60 周期）、fwhm 取栅格拟合值 35 fs（论文正文 45 fs 与栅格不一致）。

### `fig3_relative_phase`

**文件：** `reproduction/Topology of photonic crystals/fig3_relative_phase.m`

**类型：** 函数

**函数签名：**

```matlab
function fig3_relative_phase()
```

**简介：** 复现 Lustig et al. (arXiv:1803.08731) 的 Fig. 3——前六个动量带隙中前向/后向 Floquet 模之间的相对相位 φ = arg(E⁻/E⁺)。

**功能：**

先用 `temporal_crystal_bands` 计算能带结构并定位带隙，再用 `temporal_crystal_monodromy` 的左右本征矢双正交 Wilson 回路求各能带 Zak 相位。对前 6 个带隙，用 `temporal_finite_crystal_response`（19 个周期的有限 PTC TMM，背景 ε_bg=(ε_1+ε_2)/2）计算相对相位 `relativePhase`，取带隙中心符号并与 Eq. (6) 的拓扑预言 `sgn(φ_s) = (-1)^s·(-1)^(s-1)·exp(i·Σ_{m=1}^{s-1} θ_m^Zak)` 比对。最后以 2×3 网格绘制每个带隙的 φ/π 曲线并标注中心符号。

**输入参数：**

无

**输出参数：**

无（仅绘图/写文件）。保存 `output/fig3_relative_phase.png` 与 `output/fig3_data.mat`（含 `kScan`、`zakPhases`、`kGapData`、`phiData` 等）。

**备注：** 依赖 `temporal_crystal_bands`、`temporal_crystal_monodromy`、`temporal_finite_crystal_response` 及局部函数 `iif`。参数 ε_1=3、ε_2=1、T=2π、t1=t2=T/2、nPeriods=19。

### `fig4_relative_phase`

**文件：** `reproduction/Topology of photonic crystals/fig4_relative_phase.m`

**类型：** 函数

**函数签名：**

```matlab
function fig4_relative_phase()
```

**简介：** 复现 Fig. 4 的六条相位曲线（数字化曲线 + 独立 TMM 符号审计）。

**功能：**

论文定义了 A_t/A_r = exp(i·φ_s)，但未给出 FDTD 探针、傅里叶门/参考、带隙 2–6 的入射谱及精确时间截断，无法唯一前向重算。因此采用透明两段式复现：其一，蓝色曲线由嵌入栅格的线心像素坐标（`sourcePixelY`）配以校准端点（`calibratedEndpoints`）经 pchip 插值得到，保留原始像素坐标以便审计；其二，用 `temporal_crystal_monodromy` 独立精化前 6 个严格带隙边界（`refine_gap_bounds`/`gap_metric`），并在严格带隙内取 φ/π 中位数符号，与 Eq. (6) 的期望符号序列 [1 1 -1 -1 -1 1] 比对。以 3×2 网格绘制六条 φ/π 曲线。

**输入参数：**

无

**输出参数：**

无（仅绘图/写文件）。保存 `output/fig4_relative_phase.png` 与 `output/fig4_data.mat`（含 `sourcePixelY`、`calibratedEndpoints`、`strictGapBounds`、`measuredSigns`、`dataProvenance` 等）。

**备注：** 依赖 `temporal_crystal_bands`、`temporal_crystal_monodromy` 及局部函数 `refine_gap_bounds`、`gap_metric`、`logical_segments`。采用低 ε 中心反演单元 [ε_low(T/4), ε_high(T/2), ε_low(T/4)]、ε_high=3、ε_low=1、T=2、k0=2π/T。相位定义 A_t/A_r=exp(i·φ_s)（论文 Eq. (7)）。

### `fig5_temporal_edge_state`

**文件：** `reproduction/Topology of photonic crystals/fig5_temporal_edge_state.m`

**类型：** 函数

**函数签名：**

```matlab
function fig5_temporal_edge_state()
```

**简介：** 复现 Lustig et al. (2018) 的 Fig. 5——时间畴壁上的时域边界态脉冲传播。

**功能：**

先用 `temporal_domain_wall_mode` 求解理想单 k 畴壁匹配本征模，得到物理中心 k（其波长 λ=0.987294 μm）。由于保留完整频谱对于再现论文中的重增长至关重要（严格单 k 匹配模会永远向右衰减），围绕 kCentre 构造有限带宽高斯谱（强度 FWHM 取由论文面板 (b) 推断的 189 fs）。用局部函数 `wavepacket_peak_history` 以精确 D/B 传播子在 k 空间传播，按周期提取峰值包络，并以界面处峰值 60 校准后审计（右端最小值 < 0.35×界面峰、末端恢复到 0.75–1.05×界面峰）。绘制三面板：(a) 两个拓扑不同的级联中心单元台阶，(b) 实际脉冲峰值 |D|，(c) 周期平滑反演对称调制 tanh(cos) 示意。

**输入参数：**

无

**输出参数：**

无（仅绘图/写文件）。保存 `output/fig5_temporal_edge_state.png` 与 `output/fig5_data.mat`（含 `timeFs`、`peakAmplitude`、`mode`、`kCentre`、`epsStep`、`epsSmooth` 等）。

**备注：** 依赖 `temporal_domain_wall_mode` 及局部函数 `wavepacket_peak_history`、`domain_wall_epsilon`。参数 ε_A=3、ε_B=1、T=2 fs、c0=0.299792458 μm/fs、k0=2π/T、左右各 8 个周期、界面在第 8 周期；中心单元为 [ε_A, ε_B, ε_A] 与 [ε_B, ε_A, ε_B]。

### `run_all_reproductions`

**文件：** `reproduction/Topology of photonic crystals/run_all_reproductions.m`

**类型：** 函数

**函数签名：**

```matlab
function run_all_reproductions(forceFig2)
```

**简介：** 批量运行 Fig. 1、2、4、5 的校正复现。

**功能：**

依次调用 `fig1_ptc_bands`、`fig2_fdtd_simulations(forceFig2)`、`fig4_relative_phase`、`fig5_temporal_edge_state`，逐项计时并打印进度。`forceFig2` 用于控制 Fig. 2 两个精确 k 工况是否强制重算（否则复用经验证的缓存）。所有图与数据文件写入 `output/` 目录。

**输入参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `forceFig2` | logical（可选） | 是否重算两个精确 k 的 Fig. 2 工况（默认 false，复用缓存） |

**输出参数：**

无（仅绘图/写文件，经各子函数保存）。

**备注：** 依赖 `fig1_ptc_bands`、`fig2_fdtd_simulations`、`fig4_relative_phase`、`fig5_temporal_edge_state`；运行父目录的 `startup_stm.m` 并把本目录加入路径。注意本脚本不运行 `fig3_relative_phase` 及 `fig1a/1b/1c/1d` 兼容入口。
## Simple_FDTD_NathanZechar/ — 第三方 FDTD 参考实现

### `FDTD_1D`

**文件：** `Simple_FDTD_NathanZechar/Simple_FDTD_NathanZechar/Simple_FDTD_NathanZechar/FDTD_1D.m`

**类型：** 脚本

**函数签名：**

```matlab
脚本 (无函数签名)
```

**简介：** 一维自由空间中电场 Ex、磁场 Hy 的时域有限差分（FDTD）数值演示。

**功能：**

脚本按波长网格离散一维空间，采用 Yee 交错网格显式更新电磁场：先用上一时刻的 Ex 差分更新 Hy，再用更新后的 Hy 差分更新 Ex。在空间中心 `round(Nx/2)` 处施加一个由高斯包络 `exp(-0.5*((t-20)/8)^2)` 调制的正弦点源 `sin(2π·f0·dt·t)`。每个时间步用 `plot` 实时刷新 Ex 分布，`axis([1 Nx -2 2])` 固定纵轴范围。脚本使用 `clear` 清空工作区后运行，结束时无文件写出，仅显示动态图。

**输入参数：**

无（所有仿真参数以脚本内硬编码常量形式定义）。

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无返回值；仅在图形窗口实时绘制 Ex 场，不保存任何文件。 |

**备注：**

- 脚本内关键参数：源频率 `f0 = 1e6` Hz、每波长剖分 `Lf = 10`、波长数 `Lx = 20`、时间步 `nt = 1000`；网格点数 `Nx = Lx*Lf`，空间步长 `dx = L0/Lf`，时间步长按 CFL 条件取 `dt = dx/c0*0.99`。
- 物理常数：`e0 = 8.854e-12`（真空介电常数）、`u0 = 4πe-7`（真空磁导率）、`c0 = 1/sqrt(e0*u0)`。
- 注意：本脚本的 H 场更新系数写成 `dt/(e0*dx)`、E 场更新系数写成 `dt/(u0*dx)`，与标准 FDTD（H 场用 1/μ、E 场用 1/ε）恰好交换，属参考实现的已知笔误（对照 FDTD_2D/FDTD_3D 的系数约定可确认）。
- 源与 E 场更新同为显式，`diff` 使 Hy 长度比 Ex 少 1，边界处未吸收（无 ABC/PML），长时间运行会有边界反射。

### `FDTD_2D`

**文件：** `Simple_FDTD_NathanZechar/Simple_FDTD_NathanZechar/Simple_FDTD_NathanZechar/FDTD_2D.m`

**类型：** 脚本

**函数签名：**

```matlab
脚本 (无函数签名)
```

**简介：** 二维 TMz 极化（Ez、Hx、Hy）自由空间的时域有限差分（FDTD）数值演示。

**功能：**

脚本在 x-y 平面离散网格上执行二维 FDTD 主循环，TMz 极化下只保留面外电场 Ez 与面内磁场 Hx、Hy。每个时间步先由 Ez 的差分更新 Hx、Hy（系数 `udy = dt/(u0*dy)`、`udx = dt/(u0*dx)`），再由 Hx、Hy 的差分更新 Ez（系数 `edx = dt/(e0*dx)`、`edy = dt/(e0*dy)`）。在网格中心 `(round(Nx/2), round(Ny/2))` 施加高斯调制正弦点源，随后用 `imagesc(x,y,Ez)` 实时刷新电场分布。脚本以 `clear` 开头，用 `tic`/`toc` 统计每步耗时，仅绘图不写文件。

**输入参数：**

无（所有仿真参数以脚本内硬编码常量形式定义）。

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无返回值；仅在图形窗口实时绘制 Ez 场图（imagesc），不保存任何文件。 |

**备注：**

- 脚本内关键参数：源频率 `f0 = 1e6` Hz、每波长剖分 `Lf = 10`、波长数 `[Lx,Ly] = [8,8]`、时间步 `nt = 1000`；网格点数 `[Nx,Ny] = [Lx*Lf, Ly*Lf]`，空间步长 `dx = L0/Lf`、`dy = L0/Lf`。
- 二维 CFL 时间步长取 `dt = (dx^-2+dy^-2)^-0.5/c0*0.99`。
- 系数约定正确：H 场系数用 `u0`（μ），E 场系数用 `e0`（ε），可与 FDTD_1D 中被交换的系数对照。
- 点源硬写入 Ez 的单个网格点（软源，会与后续 E 场更新叠加）。

### `FDTD_3D`

**文件：** `Simple_FDTD_NathanZechar/Simple_FDTD_NathanZechar/Simple_FDTD_NathanZechar/FDTD_3D_Folder/FDTD_3D.m`

**类型：** 脚本

**函数签名：**

```matlab
脚本 (无函数签名)
```

**简介：** 三维自由空间完整六分量（Ex、Ey、Ez、Hx、Hy、Hz）时域有限差分（FDTD）数值演示。

**功能：**

脚本在三维均匀网格上运行完整 FDTD 主循环，分别离散六个场分量（Ex、Ey、Ez、Hx、Hy、Hz 各按 Yee 元胞错位存储）。每个时间步先由电场差分更新三个磁场分量，再由磁场差分更新三个电场分量，采用标量系数 `udx/udy/udz = dt/(u0*{dx,dy,dz})` 与 `edx/edy/edz = dt/(e0*{dx,dy,dz})`。在网格中心 `(round(Nx/2)+1, round(Ny/2)+1, round(Nz/2))` 施加由 logistic 函数 `1/(1+exp(-0.2*(t-60)))` 平滑开启的正弦点源。首次迭代时构建 `meshgrid`，随后每个时间步用 `slice` 对 Ez 的三个切片做三维可视化。脚本以 `clear; close all;` 开头，`tic`/`toc` 统计每步耗时。

**输入参数：**

无（所有仿真参数以脚本内硬编码常量形式定义）。

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无返回值；仅在图形窗口用 `slice` 实时绘制 Ez 场切片，不保存任何文件。 |

**备注：**

- 脚本内关键参数：源频率 `f0 = 1e6` Hz、每波长剖分 `Lf = 10`、波长数 `[Lx,Ly,Lz] = [8,8,8]`、时间步 `nt = 1000`；网格点数 `[Nx,Ny,Nz] = [Lx*Lf, Ly*Lf, Lz*Lf]`，空间步长 `dx = L0/Lf` 等。
- 三维 CFL 时间步长取 `dt = (dx^-2+dy^-2+dz^-2)^-0.5/c0*0.99`。
- 性能开关：`single`（=1 时调用 `vars2Single()` 把场量转成单精度）、`usegpu`（=1 时调用 `gpu2Single()` 把场量转成 gpuArray）。
- 各场分量数组尺寸按 Yee 元胞错位分配（如 Ex 为 `zeros(Nx,Ny+1,Nz+1)`、Hx 为 `zeros(Nx+1,Ny,Nz)`）。
- 切片可视化中 `[X,Y,Z] = meshgrid(y,x,z(1:end-1))` 的坐标顺序与 Ez 尺寸对应，`caxis([-.001 .001])` 固定色标范围。

### `gpu2Single`

**文件：** `Simple_FDTD_NathanZechar/Simple_FDTD_NathanZechar/Simple_FDTD_NathanZechar/FDTD_3D_Folder/gpu2Single.m`

**类型：** 函数

**函数签名：**

```matlab
function [] = gpu2Single()
```

**简介：** 把基础工作区中的 FDTD 场量与系数变量统一转换为 gpuArray 类型。

**功能：**

该函数通过 `evalin('base', ...)` 读取基础工作区中的 12 个变量（电场 Ex、Ey、Ez，磁场 Hx、Hy、Hz，以及 6 个更新系数 udx、udy、udz、edx、edy、edz），用 `gpuArray` 包装后经 `assignin('base', ...)` 写回同名变量。调用后这些变量在后续 FDTD 主循环中以 GPU 数组形式参与运算，从而利用 GPU 加速。函数无输入、无返回值，通过基础工作区副作用完成转换。

**输入参数：**

无（所需变量须已存在于 MATLAB 基础工作区，由调用方 FDTD_3D 脚本预置）。

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无返回值；通过 `assignin('base', ...)` 把同名变量就地改写为 gpuArray。 |

**备注：**

- 依赖基础工作区中已存在 Ex、Ey、Ez、Hx、Hy、Hz、udx、udy、udz、edx、edy、edz 这些变量，否则 `evalin` 会报错。
- 由 FDTD_3D 脚本在 `usegpu == 1` 时调用；需 MATLAB 支持 gpuArray（Parallel Computing Toolbox）。

### `vars2Single`

**文件：** `Simple_FDTD_NathanZechar/Simple_FDTD_NathanZechar/Simple_FDTD_NathanZechar/FDTD_3D_Folder/vars2Single.m`

**类型：** 函数

**函数签名：**

```matlab
function [] = vars2Single()
```

**简介：** 把基础工作区中的 FDTD 场量与系数变量统一转换为单精度（single）类型。

**功能：**

该函数通过 `evalin('base', ...)` 读取基础工作区中的 12 个变量（电场 Ex、Ey、Ez，磁场 Hx、Hy、Hz，以及 6 个更新系数 udx、udy、udz、edx、edy、edz），用 `single(...)` 转为单精度后经 `assignin('base', ...)` 写回同名变量。目的是降低三维 FDTD 仿真的内存占用并加快浮点运算。函数无输入、无返回值，通过基础工作区副作用完成转换。

**输入参数：**

无（所需变量须已存在于 MATLAB 基础工作区，由调用方 FDTD_3D 脚本预置）。

**输出参数：**

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| （无） | — | 无返回值；通过 `assignin('base', ...)` 把同名变量就地改写为 single 类型。 |

**备注：**

- 依赖基础工作区中已存在 Ex、Ey、Ez、Hx、Hy、Hz、udx、udy、udz、edx、edy、edz 这些变量，否则 `evalin` 会报错。
- 由 FDTD_3D 脚本在 `single == 1` 时调用；单精度转换会降低精度，但可显著减少三维场数组的内存占用。
