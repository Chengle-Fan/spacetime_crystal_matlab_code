# 07 · 拓扑不变量精讲

> 目标读者：**懂能带、懂 Berry 相位/拓扑的物理图景，但没在代码里见过"拓扑不变量到底怎么算"**的读者。
> 这一章把 `topology/` 目录的四个函数讲透：它们算什么、输入输出是什么、两个 demo 在多大程度上"验证了"物理——以及**在多大程度上没有**。

读完这一章，你应该能做到：

1. 用一句话说清 Zak 相位与 Chern 数的区别，并指出它们在代码里各自的"参数环"是什么；
2. 看懂 Wilson loop 的离散化：逐 k 点链接（link）+ 缝矩阵（sewing）闭合环 + 双正交内积；
3. 看懂 `zak_phase_biorthogonal`、`stpwe_bz_sewing_matrix`、`stpwe_track_band`、`fhs_chern_number` 四个函数的输入输出；
4. **不被 demo 标题带偏**：知道 `demo07` 的 Zak 是"空间"的、`demo08` 的 $C=1$ 是 Rice-Mele 算法基准，而不是时空介质的 Chern 数。

前置知识：MATLAB 语法见 [01-matlab-primer.md](01-matlab-primer.md)；ST-PWE 能带怎么算的见学习路线里的 04 章（[README.md](README.md)）。这一章只用到 `stpwe_build_system` 和 `stpwe_solve_omega` 的几个输出字段：`sol.R`、`sol.L`、`sol.omega`、`sol.m0Weight` 和 `sys.Bomega`。你只需先记住一句话：**`sol.R`/`sol.L` 是每个 k 点解出的左右本征向量，`sys.Bomega` 是广义本征问题里的 $B$ 矩阵**。这些名词在本章 2.3 会用到。

---

## 一、为什么要拓扑：一个"整体相位"变成了整数

### 1.1 体-边对应的一句话直觉

你在能带课上学过：一维晶体的能带可以有**非平凡拓扑**（比如 SSH 模型的绕数 $\nu=1$）。这个"体"（bulk）性质并非摆设——**体-边对应**（bulk-boundary correspondence）说：

> 体拓扑不变量非平凡 ⇔ 这个带隙里必然存在局域在边界/界面的态。

这就是为什么人们花力气去算"一个积分出来的整数"：它预先决定了材料边界上有没有"甩不掉"的态，跟具体边界怎么造无关。时空晶体里同样有这个概念（体带隙内的界面态，见 `demo10`），而这一章只讲"体"那一半：**把那个整数算出来**。

### 1.2 Zak 相位 vs Chern 数：一句话区分

| 不变量 | 一句话 | 参数环 | 代码函数 |
|--------|--------|--------|----------|
| **Zak 相位** | 一条**一维闭环**上的"带内几何相位"（Berry 相位沿布里渊区积一圈） | $k$ 绕布里渊区一周（1D） | `zak_phase_biorthogonal` |
| **Chern 数** | 一个**二维闭合面**上的"带内不变量"（Berry 曲率在环面上积分） | 二维参数环面，如 $(k, \phi)$（2D） | `fhs_chern_number` |

差别本质是积分维数：Zak 是"绕一圈积累了多少相位"，Chern 是"在闭合曲面上积累了多大曲率"，后者除以 $2\pi$ 得到一个整数。

> ⚠️ 注意名字：**Zak 相位**只是"沿布里渊区的几何相位"的一个具体名字，它不等于"时间的拓扑量"。`demo07` 算的是**空间 Bloch 动量** $k$ 上的 Zak（见第五章），别把"Zak"和"时间晶体"自动绑定。

---

## 二、Zak 相位：空间布里渊区上的 Wilson loop

### 2.1 连续图像 → 离散化的 Wilson loop

理想化的定义：令 $|u_n(k)\rangle$ 是第 $n$ 条带的 Bloch 本征向量，则这条带的 Zak 相位是

$$
\gamma_n = -i \int_{-\pi/\Lambda}^{\pi/\Lambda} \langle u_n(k) | \partial_k u_n(k) \rangle \, dk .
$$

你大概也见过它被改写成 Berry 相位、并且（对对称元胞）取量子化值 $0$ 或 $\pi$。

但数值上我们只有**离散的** $k$ 采样点 $k_1, k_2, \dots, k_{N_k}$。离散化的标准做法是 **Wilson loop**：逐点计算相邻两个本征向量之间的"链接相位"，然后乘起来取辐角：

$$
U_j = \frac{\langle \psi(k_j) | \psi(k_{j+1}) \rangle}{\big|\langle \psi(k_j) | \psi(k_{j+1}) \rangle\big|} ,\qquad
\gamma = -\arg \prod_{j=1}^{N_k} U_j .
$$

把 $\langle u|\partial_k u\rangle$ 的积分换成"相邻离散点内积的相位累积"，这在物理上等价，在数值上更稳健。图示：

```
 k_1 ──U_1──> k_2 ──U_2──> ... ──> k_Nk
   ^                               │
   └────────── 最后一个 U 用 sewing ┘
              (缝合矩阵接回 k_1)
```

### 2.2 缝矩阵（sewing matrix）：把环闭合

这里有个物理细节：$k=-\pi/\Lambda$ 和 $k=+\pi/\Lambda$ 是**同一个物理态**，但数值上本征向量 $u_n(k)$ 差一个下标平移。回忆 Floquet 展开

$$
E(x,t) = e^{ikx} \sum_n u_n(k)\, e^{i n g x},
$$

把 $k$ 平移一个倒格矢 $g=2\pi/\Lambda$，等价于把空间谐波下标整体加一：

$$
u_n(k+g) = u_{n+1}(k).
$$

这就是 `stpwe_bz_sewing_matrix` 干的事（[stpwe_bz_sewing_matrix.m](../topology/stpwe_bz_sewing_matrix.m)）：构造一个 $2S\times 2S$ 的置换矩阵 `sewing`，使得最后一个链接写成

$$
\text{link}_{N_k} = L^\dagger(k_{N_k})\, \text{metric}\, \text{sewing}\, R(k_1),
$$

把 $k=+g/2$ 末点的态"缝合"回 $k=-g/2$ 首点的等价态，Wilson 环就真正闭合了。函数的注释里还点明了一个收敛性前提：有限 Fourier 截断只丢掉最外层系数，所以**截断处系数权重必须可忽略**，缝合才可靠。

> 只看一维的读者可能觉得多此一举——但在代码里，如果你漏掉 `sewing`（默认是单位阵），闭合环的最后一步会把两个**物理上不同**的态连起来，Zak 相位立刻失真。`demo07` 里每次调用 `zak_phase_biorthogonal` 都把 `stpwe_bz_sewing_matrix(sys)` 显式传进去，就是这个原因。

### 2.3 为什么必须用双正交内积，而不是普通投影？

这是新手最容易困惑、也最值得停下来理解的一步。

ST-PWE 解的是**广义、非 Hermitian** 本征问题

$$
A\,\mathbf{R} = \omega\, B\,\mathbf{R},
$$

其中 $B=\texttt{sys.Bomega}$，$\omega$ 可以取复值。对这样的问题，普通"右—右"欧氏内积 $\mathbf{R}_i^\dagger \mathbf{R}_j$ **不再是好度量**：非 Hermitian 系统的右本征向量彼此不正交，用 $\mathbf{R}_i^\dagger \mathbf{R}_j$ 做投影会把相邻带"混进来"，算出的相位是错的。

正确的配对是**左右本征向量的双正交**：$\mathbf{L}_i^\dagger B \mathbf{R}_j \propto \delta_{ij}$。`stpwe_solve_omega` 已经按 $\mathbf{L}_i^\dagger B \mathbf{R}_i = 1$ 做了归一化（见 [stpwe_solve_omega.m](../core/stpwe_solve_omega.m)），所以相邻 k 点的"正确重叠"是

$$
\langle \psi(k_j)|\psi(k_{j+1})\rangle_{\text{bi}} = L^\dagger(k_j)\, B\, R(k_{j+1}).
$$

对应到代码里：`zak_phase_biorthogonal` 的第三个参数 `metric` 就是 `sys.Bomega`。[本文章节 2.3](#23-为什么必须用双正交内积而不是普通投影)也明确写了这一点——v2 的带追踪用"前后双正交重叠的几何平均"，**不再用普通右—右欧氏重叠**选择非 Hermitian 带。

> 💡 什么时候 `metric` 是单位阵？当你的本征问题本来是 Hermitian、左右本征向量重合时，$B$ 退化成单位阵，双正交内积就退回普通内积（这正是 `demo08` 的情形）。非 Hermitian 的 ST-PWE 必须显式传 `sys.Bomega`。

### 2.4 读 `zak_phase_biorthogonal` 的输入输出

函数签名（[zak_phase_biorthogonal.m](../topology/zak_phase_biorthogonal.m)）：

```matlab
[zak, info] = zak_phase_biorthogonal(rightStates, leftStates, metric, sewing)
```

| 参数 | 形状 | 含义 |
|------|------|------|
| `rightStates` | `[Nbasis, Nk]` | 单条带的右本征向量，每列一个 k 点（**k 网格不含重复右端点**） |
| `leftStates` | `[Nbasis, Nk]` | 对应的左本征向量 |
| `metric` | `[Nbasis, Nbasis]` | 重叠度量，ST-PWE 场景下传 `sys.Bomega`；默认单位阵 |
| `sewing` | `[Nbasis, Nbasis]` | 闭合 BZ 的缝合矩阵；默认单位阵，但 ST-PWE 应传 `stpwe_bz_sewing_matrix(sys)` |

内部动作只有几行：对 $j=1..N_k-1$ 算 `links(j) = L(:,j)'*metric*R(:,j+1)`，最后一步用 `links(Nk) = L(:,Nk)'*metric*sewing*R(:,1)` 闭合；若任一 `|link| < 1e-12` 直接报错（提示能带可能与其他带相触、或 Fourier 截断太小）；否则 `wilsonLoop = prod(links./abs(links))`，最后

```matlab
zak = wrap_to_pi(-angle(wilsonLoop));   % wrap_to_pi(v) = mod(v+pi, 2*pi) - pi
```

所以 `zak` 被折到 $[-\pi,\pi]$ 区间（代码里按 $\pi$ 的倍数打印，如 `zak/pi`）。返回的 `info` 结构体：

| `info` 字段 | 含义 |
|-------------|------|
| `info.links` | 每个链接的原始复数值（$N_k$ 个） |
| `info.unitLinks` | 归一化到单位模的链接相位 |
| `info.wilsonLoop` | 所有 `unitLinks` 的乘积（一个单位复数） |
| `info.minimumLinkMagnitude` | `min(abs(links))`，**最弱链接的模**，是可靠性诊断 |

**怎么读 `info.minimumLinkMagnitude`？** 健康的单带 Wilson loop 里每个链接模都应接近 1；一旦某处模明显塌陷（`demo07` 用阈值 `0.1`），说明这条带在那里几乎碰到别的带，单带 Zak 相位失去拓扑意义——`demo07` 会打警告，提示"只能当收敛性诊断，不能当鲁棒不变量"。

### 2.5 Zak 依赖单胞原点（demo07 注释的原话精神）

你可能会问：教科书说对称元胞的 Zak 是 0 或 $\pi$，为什么 `demo07` 的数值往往落在中间值？

答案写在 `demo07_zak_phase_stpwe.m` 头注释里（大意）：

> Zak 相位依赖单胞原点。Fig. 2 的单胞起点在**层边界**（layer boundary）而非反演中心（inversion center），所以在这个坐标约定下它不必恰好是 0 或 $\pi$。

也就是说：**Zak 相位不是一个"绝对"的数，它依赖你从哪个起点切单胞**。只有当你把单胞原点选在反演中心（或时间反演对称点）时，它才被强制到 $0/\pi$ 二值。`demo07` 的 Fig. 2 介质（前 3/4 静态、后 1/4 调制）从层边界起切，没有这个对称性，因此得到一般值——这是物理，不是 bug。

---

## 三、缝矩阵与带追踪：Zak 的前置工序

Zak 相位要求你手里有"某一条带"在 $N_k$ 个 k 点上的本征向量。但 `stpwe_solve_omega` 在每个 k 点解出的是 $2S$ 个本征向量——哪一个是"同一条带"？这个问题由 `stpwe_track_band` 回答。

### 3.1 `stpwe_bz_sewing_matrix`（[stpwe_bz_sewing_matrix.m](../topology/stpwe_bz_sewing_matrix.m)）

签名 `sewing = stpwe_bz_sewing_matrix(sys)`。输入一个 `sys` 结构体，返回 $2S\times 2S$ 的缝合矩阵。物理就是 2.2 节的 $u_n(k+g)=u_{n+1}(k)$：构造 $S\times S$ 的置换 `shift`（第 `row` 行第 `col` 列取 1，当且仅当 `mList(row)==mList(col)` 且 `nList(col)==nList(row)+1`），然后 `sewing = blkdiag(shift, shift)`——电场块与磁场块各自平移下标，时间谐波 $m$ 不变。

### 3.2 `stpwe_track_band`：靠"特征向量连续性"追同一条带

签名（[stpwe_track_band.m](../topology/stpwe_track_band.m)）：

```matlab
band = stpwe_track_band(sys, kGrid, omegaSeed, options)
```

它沿 `kGrid` 逐点调用 `stpwe_solve_omega`，从每个 k 点的 $2S$ 个候选模式里挑出"和上一点最连续"的那一个。`options` 是过滤条件：

| `options` 字段 | 默认值 | 作用 |
|----------------|--------|------|
| `options.omegaWindow` | `[-inf inf]` | 只保留 `real(omega)` 在这个区间的模式 |
| `options.maxImag` | `inf` | 只保留 `|Im(omega)|` 不超过它的模式 |
| `options.minM0Weight` | `0` | 只保留 m=0 参与度不低于它的模式 |

挑选的代价函数分两段：

- **首点**（`ik=1`）：按频率最近邻选种子，`cost = abs(omega-omegaSeed)/scale + 0.01*(1-m0Weight)`，`scale = max(abs(omegaSeed), sys.Omega)`；
- **后续点**：以特征向量连续性为主。对每个候选 $q$ 计算

$$
\text{fidelity}(q) = \sqrt{\big|L_{\text{prev}}^\dagger B R_q\big|\cdot\big|L_q^\dagger B R_{\text{prev}}\big|},
$$

即**正反向双正交重叠的几何平均**——注释里特意说明，这个量对左右本征向量的互易重新标度（reciprocal rescaling）不变，所以是一个规范安全的连续性度量。再加上频率差项和 m=0 权重项凑成总代价，取最小者。

选中模式后还会做一步**规范固定**（gauge fixing）：令 `phase = exp(-1i*angle(forward))`，把相邻 k 点的重叠拨成正实数，消除逐点的 U(1) 相位歧义——Wilson loop 计算要求规范固定，否则相位累积没有意义。

返回的 `band` 结构体字段：

| 字段 | 含义 |
|------|------|
| `band.k` | k 网格 |
| `band.omega` | 追踪到的频率（可能复） |
| `band.R` / `band.L` | 逐 k 点的左右本征向量（`[nState, Nk]`） |
| `band.m0Weight` | m=0 参与度 |
| `band.selectedId` | 每个 k 点选中的候选下标 |
| `band.continuity` | 每点与上一点的重叠 fidelity，标识追踪质量 |
| `band.neighborGap` | **选中模式与最近其他候选模式的频率间隔** |

`neighborGap` 是 [04 章 §五](04-st-pwe.md#五收敛与截断怎么知道算得够不够) 明确点名的重要诊断：**间隔过小说明带隙闭合**，单带 Zak 相位不可靠，应改用复合子空间 Wilson loop（见第五章 5.4）。

---

## 四、FHS Chern 数：二维环面上的格点规范

### 4.1 从连续 Berry 曲率到格点规范

Chern 数定义在**二维闭合参数面**（环面）上：把 Berry 曲率 $F$ 在面上积分，除以 $2\pi$。数值上如何从一个离散网格算出这个整数，是 Fukui–Hatsugai–Suzuki（FHS）1985 年那篇文章解决的"格点规范"问题：每个小格子（plaquette）贡献一个曲率角，全部求和除以 $2\pi$。

### 4.2 `fhs_chern_number` 的输入

签名（[fhs_chern_number.m](../topology/fhs_chern_number.m)）：

```matlab
[chern, curvature, links] = fhs_chern_number(rightStates, leftStates, metric)
```

| 参数 | 形状 | 含义 |
|------|------|------|
| `rightStates` | `[Nbasis, Noccupied, Nk, Np]` | 右本征态在二维周期网格 $(k, \phi)$ 上的排列 |
| `leftStates` | 同 `rightStates` | 双正交左本征态；**Hermitian 问题可省略**（默认等于 `rightStates`） |
| `metric` | `[Nbasis, Nbasis]` | 重叠度量；Hermitian 问题可省略（默认单位阵），非 Hermitian 场景传 `sys.Bomega` |

关键约定：**两个方向都周期闭合，端点不重复**（代码里 `ikNext = mod(ik,Nk)+1`，`ip` 方向同理）。

### 4.3 算法逐层读

**第 1 层：每个格点算两个方向上的 link 相位。** 对每个 $(ik,ip)$，构造占据子空间的重叠矩阵并取行列式（当 $N_{\text{occupied}}>1$ 时，行列式把多带子空间压成一个 U(1) 相位）：

$$
U_k(k_i,\phi_j) = \frac{\det\!\big(L^\dagger(k_i,\phi_j)\,M\,R(k_{i+1},\phi_j)\big)}
{\big|\det\!\big(L^\dagger(k_i,\phi_j)\,M\,R(k_{i+1},\phi_j)\big)\big|},
$$

$U_\phi$ 类似地把 $k$ 换到 $(k_i,\phi_{j+1})$。代码里若行列式模小于 `1e-13` 会报错"link 奇异"——对应带交叉或网格太粗。

**第 2 层：plaquette 曲率。** 每个小格子由四条 link 围成（`linkK(ik,ip)*linkP(ikNext,ip)/(linkK(ik,ipNext)*linkP(ik,ip))`），取辐角得到该格点的 Berry 曲率：

```matlab
plaquette = linkK(ik,ip)*linkP(ikNext,ip) / (linkK(ik,ipNext)*linkP(ik,ip));
curvature(ik,ip) = angle(plaquette);
```

```
        (ik,   ip) ──U_k──> (ik+1, ip)
           |                  ↑
          U_p              U_p
           ↓                  |
        (ik,   ip+1) ──U_k──> (ik+1, ip+1)
```

**第 3 层：Chern 数。** 所有格点曲率角求和除以 $2\pi$：

```matlab
chern = sum(curvature(:))/(2*pi);
```

每个曲率角都在 $(-\pi,\pi]$，网格足够细时求和结果是精确的整数。

返回的第三个输出 `links` 是结构体：`links.k` 和 `links.parameter` 分别是两个方向的归一化 link 相位矩阵，方便你复检。

### 4.4 两个默认参数何时是"单位阵"

`fhs_chern_number` 的两个可选参数 `leftStates`、`metric` 默认都取平凡值（`leftStates=rightStates`、`metric=eye`），这对应 **Hermitian 问题**——正是 `demo08` 的 Rice-Mele 模型：哈密顿量是 Hermitian 的，`eig` 给的本征向量直接可用，普通内积即可。**只有当你把非 Hermitian 的 ST-PWE 本征态送进来时**，才需要同时传 `leftStates` 和 `metric=sys.Bomega`，让每个 link 形如 `det(L'*metric*R)`。这一点在第五章的澄清里要再次出现。

---

## 五、重要澄清：demo 的拓扑量 vs 它验证了什么

这是本章最重要的部分。两个 demo 的**函数名**可能会让你以为"工具箱已经算出了时空晶体的 Chern 数 / 时间 Zak 相位"——**并没有**。下面用代码和文档原话逐条澄清。

### 5.1 `demo08`：$C=1$ 是 Rice-Mele 泵浦基准，不是时空介质的 Chern 数

`demo08_fhs_chern_thouless_pump.m` 头注释的原话精神：

> A Rice-Mele cycle provides a compact, independently quantized validation of the topology routine. The modulation phase is a synthetic dimension. **It is not a Chern-number reproduction** for demos 01–07's continuum Maxwell space-time medium.

翻译过来：

- `demo08` 用的哈密顿量是 **Rice-Mele 两带模型**，一个 2×2 的紧束缚矩阵，跟 Maxwell 方程、跟 `demo01`–`demo07` 的连续时空介质**毫无关系**；
- 它的"调制相位" $\phi$ 是**合成维度（synthetic dimension）**——用来凑出一个二维环面 $(k,\phi)$，而不是时空晶体里的时间；
- 它验证的是 **`fhs_chern_number` 这个算法实现本身正确**：格点规范取相正确、能给出精确量子化整数（打印 `Chern number = 1.000000000000`）。

[本文章节 五](#五重要澄清demo-的拓扑量-vs-它验证了什么) 把这点说得更硬：

> `demos/demo08` 的 $C=1$ 来自 Rice–Mele 泵，是对通用 FHS 实现的独立测试。**除非把目标 Maxwell/电路模型的本征态真正送入同一算法，否则这个整数不代表目标时空介质。**

### 5.2 `demo07`：算的是空间 Bloch 动量的 Zak，不是"时间 Zak 相位"

`demo07_zak_phase_stpwe.m` 头注释的原话精神：

> IMPORTANT: this loop integrates over **spatial Bloch momentum k**. It is **not the temporal Zak phase** of a spatially uniform photonic time crystal, whose Bloch coordinate is quasifrequency and whose eigenvalue is k.

也就是说，`demo07` 的闭合环路是**空间布里渊区里的 $k$**（2.1 节那个环），对应空间 Bloch 带。而光子时间晶体文献里的"时间 Zak 相位"以**准频率**为时间倒空间坐标，需要固定准频率求 $k$ 带、保持时间单胞与反演中心一致，再做 Wilson 环或反演本征值判据。[本文章节 5.2](#52-demo07算的是空间-bloch-动量的-zak不是时间-zak-相位) 明确说：**当前 v2 尚未把这一量冒充为已完成**（"当前 v2 尚未把这一量冒充为已完成"）。

### 5.3 对照表（防止误读）

| demo | 它计算的量 | 参数环 / 几何对象 | 这个数验证了什么 | 它**没有**声明什么 |
|------|-----------|-------------------|------------------|---------------------|
| `demo08` | FHS Chern 数 $C\approx 1$ | $(k,\phi)$ 环面：$k$ 是晶格动量，$\phi$ 是 Rice-Mele 调制相位（**合成维度**） | FHS 算法实现正确：格点规范取相、量子化整数 | 时空连续介质（`demo01`–`demo07` 的 Maxwell 介质）的 Chern 数——它根本没把介质的本征态送进 `fhs_chern_number` |
| `demo07` | 空间 BZ 的 Zak 相位 | 空间 Bloch 动量 $k$ 绕布里渊区一周 | 空间 Bloch 带的单带 Wilson loop 工作流（静态参照 + 驱动）与 `info.minimumLinkMagnitude` 诊断 | "时间 Zak 相位"；且当 `minimumLinkMagnitude<0.1` 时，连"鲁棒不变量"都算不上，只能当收敛性诊断 |

### 5.4 什么时候单带不变量会失效

[本文章节 五](#五重要澄清demo-的拓扑量-vs-它验证了什么) 末尾与 `stpwe_track_band` 注释共同强调：

> 遇到**异常点（exceptional point）、带简并或非孤立带**时，单带 Berry/Zak 相位通常失效，应改用**复合子空间 Wilson loop**，并报告最小能隙（`neighborGap`）和链接奇异值。

判定口诀：**看 `info.minimumLinkMagnitude` 和 `band.neighborGap`**。链接模接近 1、`neighborGap` 大 → 单带不变量可信；链接模塌陷或 `neighborGap` 趋零 → 带隙闭合，数值不可信。`demo07` 把这两个诊断全部打印出来，正是为了让你学会"先验尸、再下结论"。

---

## 六、可运行例子

下面两段分别复刻 `demo07` 的静态参照路径和 `demo08` 的完整调用序列，是**可以直接复制运行**的最小版本。

### 6.1 Zak 相位：状态怎么造、怎么调用

状态（`band.R`、`band.L`）不是手工造的，而是用 `stpwe_track_band` 从 `stpwe_solve_omega` 的结果里"追"出来的：

```matlab
startup_stm;                        % 把 core/topology 等加入路径
p = stm_fig2_parameters();          % Fig.2 介质参数（Lambda = c0 = 1）

Nk = 161;                           % 布里渊区采样点数
kGrid = -p.g/2 + (0:Nk-1)*p.g/Nk;   % k 从 -g/2 到 +g/2，不含重复右端点
omegaSeed = 0.344*p.g*p.c0;         % 种子频率，用于选第一条带

pStatic = p;  pStatic.modDepth = 0; % 静态参照：调制深度清零
sys = stpwe_build_system( ...       % 构造 ST-PWE 系统矩阵
    @(m,n) stm_fig2_eps_coeff(m,n,pStatic), [], ...
    12, 0, p.g, p.Omega);           % Nspace=12, Mtime=0

opt.omegaWindow = [0.25 0.70]*p.g*p.c0;  % 只保留该频率窗口内的模式
opt.maxImag     = 1e-7*p.g*p.c0;         % 只保留近实数模式
opt.minM0Weight = 0;                     % 不限制 m=0 参与度
band = stpwe_track_band(sys, kGrid, omegaSeed, opt);  % 逐 k 追踪一条带

sewing = stpwe_bz_sewing_matrix(sys);    % BZ 首尾缝合矩阵
[zak, info] = zak_phase_biorthogonal( ...
    band.R, band.L, sys.Bomega, sewing); % 双正交 Zak 相位

fprintf('Zak = %.4f*pi, min link = %.3e\n', ...
    zak/pi, info.minimumLinkMagnitude);
```

运行完整版 `demo07_zak_phase_stpwe`（约 1–3 分钟）会额外画出四张子图：追踪到的能带、驱动带的 `Im(omega)` 稳定性、m=0 参与度、以及 `|Wilson link|` 随 k 的分布——最后一张图就是"链接有没有塌陷"的直观体检。

### 6.2 FHS Chern 数：4 维 `states` 怎么造、怎么调用

`states` 的形状是 `[Nbasis, Noccupied, Nk, Np]`：前两维是每个 $(k,\phi)$ 点的本征向量及其占据带，后两维是二维周期网格：

```matlab
Nk = 61;  Np = 61;                 % (k, 调制相位) 两个方向各 61 点
kGrid = -pi + (0:Nk-1)*2*pi/Nk;    % k 周期环，不含重复右端点
phaseGrid = (0:Np-1)*2*pi/Np;      % 调制相位，同样是一个周期环

t0 = 1;  delta0 = 0.6;  mass0 = 1.0;
states = complex(zeros(2, 1, Nk, Np));  % [Nbasis, Nocc, Nk, Np]

for ik = 1:Nk
    for ip = 1:Np
        phase = phaseGrid(ip);
        t1 = t0 + delta0*cos(phase);     % 最近邻跳跃随相位变化
        t2 = t0 - delta0*cos(phase);
        mass = mass0*sin(phase);         % 交错势也随相位变化
        offDiag = t1 + t2*exp(-1i*kGrid(ik));
        H = [mass, offDiag; conj(offDiag), -mass];  % Rice-Mele 哈密顿量
        [V, D] = eig(H, 'vector');
        [~, order] = sort(real(D));      % 按能量排序
        states(:, 1, ik, ip) = V(:, order(1));  % 取基态（占据带）
    end
end

chern = fhs_chern_number(states);   % Hermitian 问题：leftStates 和 metric 都省略
fprintf('Chern number = %.12f\n', chern);  % 应精确等于 1
```

`demo08` 还会画三张图：FHS 曲率热图（标题直接打印 $C$）、展开后的 Zak 相位/2π 随相位（一条从 0 到 1 的直线，即泵浦量）、能隙热图。这个 demo 秒级跑完，是验证你环境就绪后最快能看到"整数"的一个。

### 6.3 怎么把这些变成你自己的研究

一个值得做的练习：把 6.2 的 `states` 换成**真正的 ST-PWE 本征向量**，构成 $(k, \phi)$ 环面（$\phi$ 是调制相位），然后调用 `fhs_chern_number(bandStates.R, bandStates.L, sys.Bomega)`。如果你能稳定地得到整数，那才真正算出了"时空介质的 Chern 数"——这正是 [本文章节 5.1](#51-demo08c1-是-rice-mele-泵浦基准不是时空介质的-chern-数) 说的"把目标模型的本征态真正送入同一算法"。在此之前，任何整数都只属于算法基准。

---

## 七、与 docs 的关系

本章只做导读，函数签名速查见 [docs/tool-reference.md](../docs/tool-reference.md)：

| 主题 | 权威出处 |
|------|----------|
| 完整的 API 用法与可照抄代码 | [本文章节 六](#六可运行例子) |
| 三种"拓扑量"的区分与适用域 | [本文章节 五](#五重要澄清demo-的拓扑量-vs-它验证了什么) |
| PWE 研究级检查（收敛、双正交、`neighborGap`） | [04 章 §五](04-st-pwe.md#五收敛与截断怎么知道算得够不够) |
| 每个函数的输入/输出/依赖速查 | [docs/tool-reference.md](../docs/tool-reference.md) 的 `topology/` 一节 |
| demo 逐个导读 | 后续 08 章（demo 导读） |

记住本章的四句话，就可以放心去跑 `demo07` 和 `demo08` 了：**Zak 是空间一维环，Chern 是二维环面；ST-PWE 的重叠要用双正交内积 `sys.Bomega`；`minimumLinkMagnitude` 和 `neighborGap` 是数值可靠性的体检报告；demo 的整数只验证算法，不自动等于时空介质的拓扑。**
