# 从这里开始：把代码包用成科研工具

这套 v2 代码不是 20 篇论文各写一个“看起来相似”的脚本，而是把它们按物理模型分开。库中 20 份 PDF 对应 17 项独立工作：其中有一篇预印本和两个相同的 Optica 版本，以及一组 2025 年预印本/正式发表版本。

先记住四条边界：

1. 空间均匀的光子时间晶体用时间 TMM/单周期矩阵最准确。
2. 空间和时间同时周期变化的连续 Maxwell 介质用 ST-PWE。
3. 传输线、超表面、合成光子网格和量子紧束缚模型都有额外自由度，不能直接替换成一个标量 $\epsilon(t)$。
4. 时间界面的连续量取决于微观切换机制。`D/B` 连续是一个常用模型，不是所有实验装置的普适定律。

## 1. 第一次运行

```matlab
startup_stm
test_smoke
demo03_ptc_bands_pwe_vs_tmm
demo09_finite_ptc_order_and_phase
demo10_temporal_domain_wall_mode
demo11_coherent_time_interface
```

这五步分别检查：

- 软件与基本代数是否正常；
- 时间 PWE 与精确单周期矩阵是否得到同一复能带；
- 无限晶体能带相同而有限时间响应相位不同；
- 共同动量带隙中是否存在增长—衰减匹配的时间畴壁模；
- 两束相干输入和不同界面跃迁律如何改变输出。

然后运行较重的收敛审计：

```matlab
demo12_ptc_convergence_audit
```

不要先从 `demo01` 的大矩阵彩色能带图入手。先把 2×2 时间单周期矩阵吃透，能够手算一个时间界面和一个周期，再进入 ST-PWE 会更稳。

## 2. 建议的学习顺序

### 阶段 A：时间界面

阅读并运行：

- `tmm/temporal_interface_matrix_jump.m`
- `demos/demo05_fdtd_temporal_interface.m`
- `demos/demo11_coherent_time_interface.m`

应当能回答：

- 为什么空间均匀的时间切换保持 $k$，而通常不保持 $\omega$？
- 对指定微观切换，究竟是 $D/B$、$E/B$ 还是更一般的状态跃迁？
- 两个输入的相对相位为什么能消除一个输出通道？

### 阶段 B：无限与有限光子时间晶体

阅读并运行：

- `tmm/temporal_crystal_monodromy.m`
- `tmm/temporal_crystal_bands.m`
- `tmm/temporal_finite_crystal_response.m`
- `demos/demo03_ptc_bands_pwe_vs_tmm.m`
- `demos/demo09_finite_ptc_order_and_phase.m`

应当能区分：

- $|\operatorname{Tr}U/2|\leq1$ 的传播区；
- $|\operatorname{Tr}U/2|>1$ 的动量带隙；
- $|\lambda|>1$ 与 $|\lambda|<1$ 的增长/衰减 Floquet 模；
- 无限晶体的谱与有限开关实验的复散射相位。

### 阶段 C：时间拓扑

先运行：

- `demos/demo10_temporal_domain_wall_mode.m`

它实现的是共同动量带隙内的 Floquet 本征态匹配。它给出一个可验证的时间局域态，但不把这个匹配量冒充 Zak 相。

当前包里的：

- `demo07_zak_phase_stpwe` 是沿空间 Bloch 波数 $k$ 的空间 Zak 相；
- `demo08_fhs_chern_thouless_pump` 是 Rice–Mele 有效模型上的算法基准。

它们都不是 Lustig 等人光子时间晶体“时间 Zak 相”的直接复现。若选择时间拓扑作为课题，第一项正式扩展应当是固定准频率求 $k$、重构一周期时间 Bloch 模，并沿时间 Brillouin 区计算对称性约束的 Wilson 环。见 `RESEARCH_ROADMAP_CN.md`。

### 阶段 D：一般时空周期介质

最后进入：

- `demo01_reproduce_fig2_stpwe`
- `demo02_complex_frequency_and_momentum_gaps`
- `demo06_fdtd_spacetime_wavepacket`
- `demo07_zak_phase_stpwe`

此时要同时检查复本征值、Fourier 截断、带追踪、独立求解器和 FDTD 网格收敛，而不只是看图形是否接近论文。

## 3. 每次形成科研结果时必须留下的证据

至少保留以下五项：

1. 参数、单位、场约定、时间因子和单胞原点；
2. 截断阶数或网格、运行时间和输出的原始 `.mat` 数据；
3. 与解析式、精确单周期矩阵或另一种算法的交叉验证；
4. 改变截断/步长后的误差曲线；
5. 结论成立所需的能隙、对称性、孤立带和界面跃迁律。

`demo12_ptc_convergence_audit` 是最小模板：它同时保存图片和误差数据，而不是只保存一张最终图。

## 4. 最适合从本包启动的课题

结合你已有的拓扑光子学基础，最推荐的主线是：

> 二元或平滑光子时间晶体的时间拓扑：严格时间 Zak 相 → 有限周期散射相位 → 时间畴壁态 → 对时间色散和非理想切换的稳健性。

这条主线可以复用当前 TMM、PWE、FDTD 和双正交拓扑工具，又有明确的“尚未完成”部分，不会停留在重复论文图片。

