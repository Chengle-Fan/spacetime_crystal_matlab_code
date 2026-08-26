# 09 · 动手练习：从改参数到设计时空晶体

> 目标读者：已经把前几章读完、看懂了 demo、但**还没自己动手改过代码**的读者。
> 这一章不教新函数——它让你把前面学过的函数真正用起来。练习分三级：**改参数**（热身，只动数字不动结构）→ **改材料**（把材料的物理模型换掉）→ **设计**（从零设计一个时空晶体，并用三个引擎互相验证）。
> 每一道练习都给了三样东西：**目标**、**起点代码**、**数值化的成功判据**。判据是你自己"做完了没有"的标尺——先想清楚"什么算对"，再动手。

---

## 一、使用说明：练习怎么练最有效

### 1.1 先跑起来，再谈优化

动手之前，请先确认两件事（如果还没做过，回 [01-matlab-primer.md](01-matlab-primer.md) 第 10 节）：

```matlab
startup_stm          % 把 core/tmm/fdtd/topology/demos/tests 加入路径
test_smoke           % 自检测试，全过会打印 All smoke tests passed.
```

记住本项目的一切都在归一化单位下运行：$\varepsilon_0=\mu_0=c_0=1$、$\Lambda=1$，时间因子约定为 $e^{ikx-i\omega t}$，介电常数按 $\varepsilon(x,t)=\sum_{m,n}\varepsilon_{m,n}e^{i n g x-i m\Omega t}$ 展开。动手前翻一遍 [02 章](02-physics-background.md)，那是最容易踩坑的地方。

### 1.2 三道题的"科学方法"：预测 → 运行 → 对照

每一道练习都建议按这个顺序做，别跳步：

```
① 读目标，写下一句预测（用文字，不写代码）
        ↓
② 跑起点代码，得到数值结果
        ↓
③ 对照"成功判据"：达标 → 写两句解释；不达标 → 回到②排查
```

**"写预测"不是形式**。改参数之前先预测"带隙会变宽还是变窄、往哪边移动"，运行后你的直觉才会被校准。这也是"什么时候可以说复现"的习惯养成。

### 1.3 改代码的三个纪律

1. **不要改原 demo 文件**。在 MATLAB 编辑器里把要改的 demo "另存为"一个自己的脚本（比如把 `demo05_fdtd_temporal_interface.m` 存成 `my_e2.m`），在原文件上改来改去，下次对照论文就找不到"原始版本"了。
2. **一次只改一处**。一次改三个参数，报错时你不知道是哪个引起的。
3. **报错不慌，按套路来**。`dbstop if error` 让 MATLAB 在出错处暂停；读第一行错误信息；不确定某个函数怎么用就 `help 函数名`。

### 1.4 三道练习的地图

```
第 3 级 · 设计      E7 造晶体 → E8 三方验证 → E9 看拓扑      （自己说了算）
第 2 级 · 改材料    E4 换 FDTD 材料 → E5 自写系数 → E6 换跳变定律
第 1 级 · 改参数    E1 截断 → E2 折射率 → E3 占空比         （只动数字）
```

| 练习 | 名称 | 难度 | 主要引擎/工具 |
|------|------|------|--------------|
| E1 | 观察 ST-PWE 的空间截断收敛 | ★ | `stpwe_build_system` / `stpwe_solve_omega` |
| E2 | 时间界面的 Morgenthaler 系数 | ★ | `fdtd1d_db` / `temporal_interface_matrix` |
| E3 | 二元 PTC 的带隙随 contrast/duty 变化 | ★ | `temporal_crystal_bands` |
| E4 | 给 FDTD 自定义材料 | ★★ | `fdtd1d_db` / `out.energy` |
| E5 | 自写满足实函数条件的 `epsCoeff` | ★★ | `stpwe_build_system` / `stpwe_solve_omega` |
| E6 | DB 与 EB 跳跃律的差别 | ★★ | `temporal_interface_matrix_jump` |
| E7 | 设计一个 PTC 并做收敛审计 | ★★★ | ST-PWE + TMM 参照 |
| E8 | 三个引擎验证动量带隙 + FDTD 增长 | ★★★ | TMM + `fdtd1d_db` |
| E9 | Zak 相位对元胞原点的依赖 | ★★★ | `stpwe_track_band` / `zak_phase_biorthogonal` |

> ⏱️ 预计总时间 3–6 小时。每一级做完可以休息，不用一次做完。

---

## 二、第 1 级 · 改参数（热身）

> 这一级只改**数字**，不改结构。目标是让你对"每个参数在控制什么"有手感，并建立"数值结果要对得上解析/参照"的直觉。

### E1 · 观察 ST-PWE 的空间截断收敛（demo01）

**目标**

`demo01_reproduce_fig2_stpwe` 用 ST-PWE 复现论文 Fig.2 的能带。它的精度由两个截断数控制：`Nspace`（空间谐波 $n=-N_{space}\dots N_{space}$）和 `Mtime`（时间谐波 $m=-M_{time}\dots M_{time}$）。谐波总数 $S=(2N_{space}+1)(2M_{time}+1)$，广义本征问题矩阵是 $2S\times 2S$。你的任务是：改变 `Nspace`，观察能带是否稳定，并感受**运行时间随截断的爆炸式增长**。

**提示**

- Fig.2 介质的调制是纯正弦 $\sin(\Omega t)$，Fourier 系数只有 $m=0,\pm1$ 非零，所以 **`Mtime=1` 在时间维已经是精确的**——这是 demo01 固定 `Mtime=1` 的原因。真正需要收敛测试的维度是 `Nspace`。
- `Nspace=10` 时 $S=63$、矩阵 $126\times126$；`Nspace=20` 时 $S=123$、矩阵 $246\times246$。`Nspace` 翻倍，矩阵尺寸和求解时间涨得很快。

**起点代码**

下面这段自己建系统、扫能带、把四个 `Nspace` 的结果叠加在一张图上（筛选条件与 demo01 完全相同）：

```matlab
%% E1：观察空间谐波截断 Nspace 的收敛
p  = stm_fig2_parameters();              % Lambda=1, c0=1, g=2π, Omega=2π*0.20
kBar = linspace(-0.5, 0.5, 101);         % 第一布里渊区
kValues = p.g*kBar;                      % 物理波矢
fMax = 0.82; weightThreshold = 0.035; imagTolerance = 3e-3;   % 与 demo01 相同
epsCoeff = @(m,n) stm_fig2_eps_coeff(m,n,p);

Nlist = [4 6 10 20];                     % 尝试的空间截断
bandK = cell(1,numel(Nlist));
bandF = cell(1,numel(Nlist));
for iN = 1:numel(Nlist)
    sys = stpwe_build_system(epsCoeff, [], Nlist(iN), 1, p.g, p.Omega);
    allK = cell(numel(kBar),1); allF = cell(numel(kBar),1);
    for ik = 1:numel(kBar)
        sol = stpwe_solve_omega(sys, kValues(ik));
        fBar = sol.omega/(p.g*p.c0);
        keep = isfinite(fBar) & real(fBar) >= 0 & real(fBar) <= fMax ...
             & abs(imag(fBar)) <= imagTolerance ...
             & sol.m0Weight >= weightThreshold;
        allK{ik} = repmat(kBar(ik), sum(keep), 1);
        allF{ik} = real(fBar(keep));
    end
    bandK{iN} = vertcat(allK{:});
    bandF{iN} = vertcat(allF{:});
end

figure('Color','w'); hold on;
for iN = 1:numel(Nlist)
    scatter(bandK{iN}, bandF{iN}, 6, 'filled', ...
        'DisplayName', sprintf('Nspace=%d', Nlist(iN)));
end
xlabel('k\Lambda/(2\pi)'); ylabel('\omega\Lambda/(2\pi c)');
legend('show','Location','best'); grid on; box on;
```

> 💡 这段代码演示了"用 `cell` 攒数据、最后 `vertcat` 拼起来"的习惯——每个 `k` 点筛出来的能带点数不一样，没法用一个矩阵预先装下，所以先用 `cell`，最后统一合并。

**数值化的成功判据**

- 四组曲线叠加后，**`Nspace=10` 与 `Nspace=20` 的 `Re(omega)` 散点几乎完全重合**；`Nspace=4`、`6` 与它们相比有可见的偏离（尤其能带拐弯处）。
- 把 `kBar = linspace(-0.5, 0.5, 101)` 改成 `31` 个点，体会运行时间随 `Nspace` 的暴涨（`Nspace=20` 最明显）。若太慢，先用 31 个点快速验证趋势，再回来跑完整版。
- 追加思考：把上面循环里的 `1`（即 `Mtime`）改成 `2` 重跑，能带图应该**几乎完全一样**——因为调制只有 $m=0,\pm1$。这解释了为什么 `Mtime` 不需要收敛测试（对比 [04 章 §五](04-st-pwe.md#五收敛与截断怎么知道算得够不够)：方波调制才需要高时间截断）。

### E2 · 时间界面的 Morgenthaler 系数随折射率跳变（demo05）

**目标**

`demo05_fdtd_temporal_interface` 用 FDTD 模拟一个时间界面：折射率在 $t_{switch}$ 时刻从 $n_{before}$ 跳到 $n_{after}$，入射波包分裂为前向（时间折射）与后向（时间反射）两支，数值系数与解析的 Morgenthaler 系数对照。你的任务是：**改 `nBefore` / `nAfter`**，验证 FDTD 的数值 `tau` / `rho` 仍然和解析公式一致。

**提示**

- 改 demo05 顶部三行：`nBefore = 1.5;` 和 `nAfter = 2.5;`。建议试 `nBefore=2.5, nAfter=4.0`，以及 `nBefore=1.0, nAfter=4.0`（反差更大）。
- 物理约定：时间界面处**波矢 $k$ 守恒、角频率改变**（$\omega_{after}=(n_{before}/n_{after})\omega_{before}$）。解析系数（Morgenthaler，DB 模型）为
$$
\tau=\frac12\Bigl(\frac{\varepsilon_1}{\varepsilon_2}+\frac{n_1}{n_2}\Bigr),\qquad
\rho=\frac12\Bigl(\frac{\varepsilon_1}{\varepsilon_2}-\frac{n_1}{n_2}\Bigr),
$$
其中 $\varepsilon=n^2$（因为 $\varepsilon_0=\mu_0=c_0=1$）。注意 $\rho$ 可为负（相位反转），比较时用 $|\rho|$。
- 改动后，`demo05` 里的 `temporal_eps` 子函数用 `nBefore^2`、`nAfter^2` 作为介电常数，`Hhalf0` 用 `nBefore` 构造初始行波——这些都会自动跟着变，你只需要改两个 `n`。

**起点代码**

先复制 demo05 为 `my_e2.m`、改好 `nBefore`/`nAfter` 并运行，得到 `out`、`nBefore`、`nAfter`、`tSwitch`。然后运行下面这段，从 `out` 里把数值系数取出来算相对误差：

```matlab
%% E2：从 out 里取数值 tau/rho，并对照解析 Morgenthaler 公式
[~, beforeId] = min(abs(out.t - (tSwitch-2*out.dt)));
[~, afterId ] = min(abs(out.t - (tSwitch+2*out.dt)));
Ebefore = out.E(beforeId,:);
Eafter  = out.E(afterId,:);
Hafter  = out.H(afterId,:);
Eplus  = 0.5*(Eafter + Hafter/nAfter);   % 前向分量（时间折射）
Eminus = 0.5*(Eafter - Hafter/nAfter);   % 后向分量（时间反射）
tauNumeric = norm(Eplus)/norm(Ebefore);
rhoNumeric = norm(Eminus)/norm(Ebefore);

[~, tauExact, rhoExact] = temporal_interface_matrix( ...
    nBefore^2, 1, nAfter^2, 1);          % 解析 Morgenthaler 系数

relErrTau = abs(tauNumeric - abs(tauExact))/abs(tauExact);
relErrRho = abs(rhoNumeric - abs(rhoExact))/abs(rhoExact);
fprintf('relErr(tau)=%.2e, relErr(rho)=%.2e\n', relErrTau, relErrRho);
```

> 💡 `Eplus = 0.5*(E + H/nAfter)` 是把界面后的场按波阻抗 $n_{after}$ 分解为前向/后向——这正是 [05 章时间传输矩阵] 里"方向态 $[E_+;E_-]$"的实空间版本。`norm(...)/norm(...)` 是在求两种波的整体幅度之比。

**数值化的成功判据**

- `relErrTau` 和 `relErrRho` 都 **< 1e-2**（量级大约在 1e-3 或更小）。
- 若反差调得太大（比如 1 → 4）误差变大，先检查是不是波包在 `tEnd=30` 之前撞上了 sponge 吸收层：把 `tEnd` 调小一点，或把波包起点 `x0` 提前。
- 追加思考：记下 $|\tau|$ 和 $|\rho|$ 随 $(n_{after}/n_{before})$ 的变化趋势——时间界面不保模、可以给场注入或抽走能量，这正是动量带隙里指数增长的"种子"（呼应 [02 物理背景] 与 demo06）。

### E3 · 二元 PTC 的动量带隙随 contrast 和 duty 变化（demo03）

**目标**

`demo03_ptc_bands_pwe_vs_tmm` 研究一个**二元光子时间晶体**：$\varepsilon(t)$ 在 $\varepsilon_A$ 与 $\varepsilon_B$ 之间按占空比 `dutyA` 方波跳变。带隙的位置与宽度由两个量决定：**contrast**（$\varepsilon_B/\varepsilon_A$，反差）和 **duty**（`dutyA`）。你的任务是改 `epsA` / `epsB` / `dutyA`，记录动量带隙（$|\mathrm{Im}\,\omega_F|>0$ 的区域）如何移动，并**用解析估计解释它**。

**提示**

- demo03 里带隙的"权威判据"是单值矩阵的 halfTrace：$|\mathrm{Tr}(U)/2|>1 \iff$ 在动量带隙内（见 [05 章时间传输矩阵] 与 `temporal_crystal_bands` 返回的 `halfTrace` 字段）。
- 弱调制近似下，第一个动量带隙中心位于 $\bar{k}_c = k_c c_0/\Omega \approx \bar n/2$，其中 $\bar n=\sqrt{\mathrm{dutyA}\cdot\varepsilon_A+(1-\mathrm{dutyA})\cdot\varepsilon_B}$ 是时间平均折射率。**先写预测再运行**：把 `epsB` 变大（contrast 增大）→ $\bar n$ 变大 → 带隙中心往**哪个方向**移？带隙是变宽还是变窄？把 `dutyA` 从 0.5 改成 0.3 呢？
- 细节：`dutyA=0.5` 时方波只有奇次谐波，偶次谐波严格为零；改 `dutyA` 会让偶次谐波出现，带隙结构会发生定性变化（不只是平移）。

**起点代码**

下面这段独立算带隙区间（不依赖 demo03 的绘图）：

```matlab
%% E3：定位二元 PTC 的动量带隙区间（|halfTrace|>1 判据）
epsA = 1.0; epsB = 4.0; dutyA = 0.5;     % ← 改这里
T = 1; Omega = 2*pi/T;
durations = [dutyA, 1-dutyA]*T;
kNorm = linspace(0.02, 1.45, 181);
bands = temporal_crystal_bands(kNorm*Omega, [epsA epsB], [1 1], durations);

inGap = abs(bands.halfTrace) > 1;        % |Tr(U)/2|>1 ⇔ 动量带隙
kGap = kNorm(inGap);
if isempty(kGap)
    fprintf('本 k 范围内没有带隙，把 kNorm 范围拉大再试。\n');
else
    fprintf('gap: kc/Omega in [%.3f, %.3f], width %.3f\n', ...
        kGap(1), kGap(end), kGap(end)-kGap(1));
end
```

**数值化的成功判据**

- 你在运行前写下的预测（带隙中心移动方向 + 宽度变化），与上面的 `gap` 输出一致；并且你能用 $\bar n=\sqrt{\mathrm{dutyA}\,\varepsilon_A+(1-\mathrm{dutyA})\varepsilon_B}$ 这个公式解释位移方向。
- 具体检查建议：固定 `dutyA=0.5`，取 `epsB` = 4、9、16 三档，打印三行的 `gap` 区间——中心应单调右移、宽度应单调变宽。
- 追加思考：为什么 `dutyA=0.5` 时带宽对 `epsB` 的变化特别敏感？提示：此时方波 Fourier 谱只有奇次谐波、幅度 $\propto 1/|m|$，带隙的"深度"由基频分量的强度决定（对照 [demo12 的收敛审计](../demos/demo12_ptc_convergence_audit.m)）。

---

## 三、第 2 级 · 改材料

> 这一级不再只动数字，而是**换掉材料的物理模型**：给 FDTD 自定义介电函数、自写 Fourier 系数、换时间界面的跳变定律。做完你会明白"材料"在代码里其实就是一个函数或一组数。

### E4 · 给 FDTD 自定义材料（demo05 基础上）

**目标**

`fdtd1d_db` 的材料完全由 `cfg.epsFun` 这个函数句柄决定。demo05 里它是一个时间阶跃函数；现在你把它换成**只含空间依赖**的介电函数（一维光子晶体），检验能量守恒——再（可选）加上微弱的时间调制，观察守恒如何被打破。

**提示**

- 无时间调制时，总电磁能量守恒，`out.energy` 应几乎不变。sponge 吸收层会吸掉撞到边界的波，所以让波包远离边界、在能量测量完成前别碰到吸收层。
- 空间周期介质的平均折射率用 $\bar n=\sqrt{\overline{\varepsilon}}$ 估计，`Hhalf0` 用它构造初始行波即可（不必精确）。
- "含三次谐波"的意思：给空间剖面再加一个 $n=3$ 的空间谐波，例如 $\varepsilon(x)=1+0.3\cos(2\pi x/\Lambda)+0.1\cos(6\pi x/\Lambda)$。对应到 ST-PWE 里就是"需要 `Nspace>=3`"（想想为什么，见 E7）。

**起点代码**

```matlab
%% E4：给 FDTD 自定义材料——静态空间光栅，检验能量守恒
p = stm_fig2_parameters();
epsSpatial = @(xq) 1 + 0.3*cos(2*pi*xq/p.Lambda);       % 空间周期，时间无关
nEff = sqrt(mean(epsSpatial(linspace(0,p.Lambda,401))));

dx = p.Lambda/40;
x  = 0:dx:60*p.Lambda;
dt = 0.55*dx/p.c0;
k0 = 0.30*p.g;  x0 = 15*p.Lambda; sigma = 2.5*p.Lambda;
profile = @(xq) exp(-((xq-x0)/sigma).^2).*exp(1i*k0*(xq-x0));
E0 = profile(x);
xH = x(1:end-1)+dx/2;                     % sponge 边界：H 网格少一个点
Hhalf0 = nEff*profile(xH + dt/2);

cfg.x = x; cfg.dt = dt; cfg.nSteps = ceil(25/dt);
cfg.epsFun = @(xq,tq) epsSpatial(xq);     % 材料：只有空间依赖
cfg.muFun  = @(xq,tq) ones(size(xq));
cfg.E0 = E0; cfg.Hhalf0 = Hhalf0;
cfg.boundary = 'sponge'; cfg.spongeCells = 120; cfg.spongeStrength = 0.08;
cfg.recordEvery = 2;
out = fdtd1d_db(cfg);

fprintf('energy(end)/energy(1) = %.4f\n', out.energy(end)/out.energy(1));
```

**数值化的成功判据**

- **静态情况**：`out.energy(end)/out.energy(1)` 与 1 的偏差 **< 5%**（偏离来自数值耗散与吸收层的轻微吸收，波包别碰到吸收层即可压到很小）。
- 可选加分：把 `epsSpatial` 改成含三次空间谐波 `+ 0.1*cos(6*pi*xq/p.Lambda)`，重跑，能量判据仍应成立；再改成时间调制 `+ 0.05*sin(p.Omega*tq)`，画出 `out.energy/out.energy(1)` 随 `out.t` 的曲线——这时能量会振荡甚至净变化，**这不是 bug**：时间调制与场交换能量，总能量不再守恒（[02 章 §5.1](02-physics-background.md#51-能量不守恒但物理没被打破)）。你的判据变成"你能用'调制做功'解释这条曲线"。
- 追问：用 `stm_stft` 对 `out.E` 的某个空间点做短时傅里叶变换，应该能看到主峰两侧距离 $\Omega$ 的 Floquet 边带（对照 [demo06](../demos/demo06_fdtd_spacetime_wavepacket.m)）。

### E5 · 自写满足实函数条件的 epsCoeff 并检查能带对称性

**目标**

ST-PWE 的材料是一个系数函数 `epsCoeff(m,n)`，返回 $\varepsilon_{m,n}$。物理上 $\varepsilon(x,t)$ 是**实函数**，这强制 Fourier 系数满足共轭关系
$$
\varepsilon_{-m,-n}=\mathrm{conj}(\varepsilon_{m,n}).
$$
你的任务：自己写一组满足这个关系的系数，验证关系本身，再检查由此得到的能带关于 $k=0$ 的对称性。

**提示**

- 最简单的"实时空介质"是行波调制 $\varepsilon(x,t)=1+\delta\cos(gx-\Omega t)$。按约定 $\varepsilon=\sum_{m,n}\varepsilon_{m,n}e^{i n g x-i m\Omega t}$ 展开，它只有两支系数：$(\varepsilon_{1,1},\varepsilon_{-1,-1})$，且两者相等（都是 $\delta/2$）。请验证它满足共轭关系。
- 物理背景：材料为实 ⇒ 若 $\omega$ 是波矢 $k$ 处的本征值，则 $\mathrm{conj}(\omega)$ 是波矢 $-k$ 处的本征值（把 Maxwell 方程整体取复共轭即可看出）。所以 $k$ 与 $-k$ 处的本征频率**实部完全相同、虚部相反**。

**起点代码**

```matlab
%% E5：自写"实"时空介质的 epsCoeff，检查能带在 ±k 的对称性
p = stm_fig2_parameters();
% 行波调制 eps(x,t) = 1 + 0.4*cos(g*x - Omega*t)
epsCoeff = @(m,n) double(m==0 && n==0) ...
    + 0.2*double((m==1 && n==1) | (m==-1 && n==-1));

% (a) 自检共轭关系：eps_{-m,-n} = conj(eps_{m,n})
err = 0;
for im = -4:4, for in = -4:4
    err = max(err, abs(epsCoeff(-im,-in) - conj(epsCoeff(im,in))));
end, end
fprintf('max conj-relation error = %.2e\n', err);      % 应为 0

Nspace = 10; Mtime = 1;
sys = stpwe_build_system(epsCoeff, [], Nspace, Mtime, p.g, p.Omega);

% (b) 比较 k 与 -k 处本征频率的实部（应完全一样）
kTest = 0.2*p.g;
fPlus  = sort(real(stpwe_solve_omega(sys, +kTest).omega));
fMinus = sort(real(stpwe_solve_omega(sys, -kTest).omega));
fprintf('max |Re(omega(+k)) - Re(omega(-k))| = %.2e\n', max(abs(fPlus-fMinus)));
```

> 💡 `sort(real(...))` 这一行很关键：`eig` 返回的本征值顺序是乱的，`+k` 和 `-k` 两组解要对齐，就得先把实部排序再逐项比。

**数值化的成功判据**

- (a) 的 `max conj-relation error` **等于 0**（或 ~1e-16）。
- (b) 的 `max |Re(omega(+k))-Re(omega(-k))|` **在 1e-9 量级或更小**（远小于能带间距）。这说明 $k$ 与 $-k$ 的实频能带一致。
- 追加实验：把 `epsCoeff` 改成**不满足**共轭关系的复数系数（例如把 `(m==1 && n==1)` 那支系数改成 `0.2i`），重新跑 (b)，对称性应被破坏——这就是"实材料"这条假设在数值上的分量。

### E6 · 比较 DB 与 EB 跳跃律的 tau / rho（呼应 demo11）

**目标**

时间界面的物理取决于"界面时刻到底什么量连续"。`temporal_interface_matrix_jump` 提供了多种理想跳变定律；`demo11_coherent_time_interface` 已经展示 DB 与 EB 两种定律给出不同的系数和不同的相干相消点。你的任务是亲手算两种定律的 `tau` / `rho`，并**给出物理解释**。

**提示**

- 通用公式（见 [05 章时间传输矩阵] 与 `temporal_interface_matrix_jump` 源码）：
$$
r_E=\mathrm{jumpD}\cdot\frac{\varepsilon_b}{\varepsilon_a},\qquad
r_M=\mathrm{jumpB}\cdot\sqrt{\frac{\varepsilon_b\mu_b}{\varepsilon_a\mu_a}},
\qquad \tau=\tfrac12(r_E+r_M),\quad \rho=\tfrac12(r_E-r_M).
$$
- `'DB'` 模型：`jumpD=1, jumpB=1`（D、B 连续，本构切换模型）；`'EB'` 模型：`jumpD=epsAfter/epsBefore, jumpB=1`（E、B 连续）。注意对 `'EB'`，$r_E = (\varepsilon_a/\varepsilon_b)\cdot(\varepsilon_b/\varepsilon_a)=1$，电场跳变被 `jumpD` "补偿"掉了。
- 对 $\varepsilon_b=1.5^2$、$\varepsilon_a=2.5^2$（demo11 的材料），预期数值：DB 给 $\tau=0.48,\ \rho=-0.12$；EB 给 $\tau=0.80,\ \rho=0.20$。跑出来对一下。

**起点代码**

```matlab
%% E6：同一材料跳变，不同跳变定律给出不同的 tau/rho
epsBefore = 1.5^2; epsAfter = 2.5^2;     % 与 demo11 相同的材料
[~, tauDB, rhoDB] = temporal_interface_matrix_jump(epsBefore,1,epsAfter,1,'DB');
[~, tauEB, rhoEB] = temporal_interface_matrix_jump(epsBefore,1,epsAfter,1,'EB');
fprintf('DB: tau=%.4f  rho=%.4f\n', tauDB, rhoDB);
fprintf('EB: tau=%.4f  rho=%.4f\n', tauEB, rhoEB);
```

**数值化的成功判据**

- 打印值与上面的预期值一致。
- **给出物理解释**：DB 模型认为开关瞬间无自由电荷/电流注入，$D$、$B$ 这两个"积分量"守恒，$E=D/\varepsilon$ 随 $\varepsilon$ 的跳变而跳变（$\varepsilon$ 变大 $\Rightarrow E$ 变小），所以透射系数偏小；EB 模型假设开关电路把 $E$ 也钉住，等价于向 $D$ 注入 $\mathrm{jumpD}$ 倍的电位移，"吸收"了 $E$ 的跳变，所以 $\tau,\rho$ 都更大。这正解释了 demo11 里 DB 与 EB 的相干相消点不同——**报实验必须声明跳变定律**（[05 章 §2.3](05-temporal-tmm.md#23-通用版本-temporal_interface_matrix_jump跳变律是第-5-个参数)）。
- 追问：`rho` 在 DB 下是负的、EB 下是正的，这个符号差对应什么？（提示：时间反射波的相位）。

---

## 四、第 3 级 · 设计

> 这一级开始"自己说了算"：设计一个晶体，然后用三种独立引擎（ST-PWE、TMM、FDTD）验证它。这正是研究工作的缩影。为了让三引擎能互相验证，我们选择**空间均匀、时间二元**的设计——它是时空晶体中唯一能被 TMM"精确参照"的一类（TMM 要求空间均匀、$k$ 守恒）。

### E7 · 设计一个光子时间晶体，并做收敛审计

**目标**

选定一组材料参数（$\varepsilon_A,\varepsilon_B,\mathrm{dutyA}$）和调制周期 $T$，用 ST-PWE 画出准频率能带、找出动量带隙，并**用 demo12 的收敛审计思想**说明你的时间截断 `Mtime` 够不够。

**提示**

- 二元方波在时间上不连续，Fourier 系数只按 $\sim 1/|m|$ 衰减，收敛比正弦调制慢得多——所以 `Mtime` 要取大（demo03 用 19）。这正是 [demo12](../demos/demo12_ptc_convergence_audit.m) 要量化的事情。
- 收敛审计的标准做法：把"空间均匀、时间二元"介质的 ST-PWE 能带，与**精确的 TMM 能带**（`temporal_crystal_bands`，无截断）逐分支对比，误差应随 `Mtime` 单调下降。
- 先抄 demo03 的参数（$\varepsilon_A=1,\varepsilon_B=4,\mathrm{dutyA}=0.5$）跑通，再改成你自己的设计。**每改一组参数，收敛审计都必须重新做一遍**——这就是研究级检查（[04 章 §五](04-st-pwe.md#五收敛与截断怎么知道算得够不够)）。

**起点代码**

```matlab
%% E7：设计空间均匀、时间二元的 PTC，并做收敛审计（demo12 的套路）
epsA = 1.0; epsB = 4.0; dutyA = 0.5;     % ★ 你的设计
T = 1; Omega = 2*pi/T;
kNorm = linspace(0.2, 1.4, 61);
kValues = kNorm*Omega;

% 精确参照：TMM（无截断）
reference = temporal_crystal_bands(kValues, [epsA epsB], [1 1], ...
    [dutyA 1-dutyA]*T);

% 收敛审计：误差应随 M 单调下降
Mvalues = [5 9 13 19];
for iM = 1:numel(Mvalues)
    Mtime = Mvalues(iM);
    epsCoeff = @(m,n) double(n==0)*temporal_binary_eps_coeff(m, epsA, epsB, dutyA);
    muCoeff  = @(m,n) double(m==0 && n==0);
    sys = stpwe_build_system(epsCoeff, muCoeff, 0, Mtime, 1, Omega);
    err = nan(2,numel(kValues));
    for ik = 1:numel(kValues)
        sol = stpwe_solve_omega(sys, kValues(ik));
        folded = stpwe_fold_frequency(sol.omega, Omega);
        ids = find(isfinite(folded) & abs(imag(folded)) < 0.45*Omega);
        used = false(size(ids));
        for br = 1:2
            target = reference.omegaF(br,ik);
            cost = abs(folded(ids)-target);  cost(used) = inf;
            [bestCost, loc] = min(cost);
            if isfinite(bestCost), err(br,ik) = bestCost/Omega; used(loc) = true; end
        end
    end
    fprintf('M=%2d: median |dw|/Omega=%.3e  max=%.3e\n', ...
        Mtime, median(err(:),'omitnan'), max(err(:),[],'omitnan'));
end

% 画能带（精确参照），肉眼看动量带隙的位置和宽度
figure('Color','w'); hold on;
plot(kNorm, real(reference.omegaF.')/Omega, 'k-', 'LineWidth',1.5);
plot(kNorm, imag(reference.omegaF.')/Omega, 'r-', 'LineWidth',1.2);
xlabel('kc/\Omega'); legend('Re(\omega_F)/\Omega','Im(\omega_F)/\Omega');
grid on; box on;
```

> 💡 这段是 demo12 的"精简版"。核心技巧：`stpwe_fold_frequency` 把 PWE 的高次 Floquet 副本折回第一时间布里渊区，再用"距离 TMM 目标最近"逐分支配对，`cost(used)=inf` 保证两条带一对一。

**数值化的成功判据**

- 打印的 `median` / `max` 误差随 `M=5,9,13,19` **单调下降**；在 `M=19` 时 `max` 误差 **< ~1e-2**（尤其注意带隙边缘处的误差——那里最难收敛）。这就是"截断够了"的证据。
- 你还能用一句话说明：为什么方波需要大 `Mtime`，而 [E1](#e1--观察-st-pwe-的空间截断收敛demo01) 的正弦调制 `Mtime=1` 就够了。

### E8 · 用 TMM 验证动量带隙，用 FDTD 观察指数增长

**目标**

E7 已经用 ST-PWE 找到了带隙。现在用另外两个引擎独立验证同一个带隙：**TMM** 给出精确的带隙区间（`halfTrace` 判据），**FDTD** 打一个中心波矢落在带隙内的波包，观察它在时间上的**指数增长**。判据：三个引擎给出的带隙区间一致。

**提示**

- TMM 的带隙判据 $|\mathrm{Tr}(U)/2|>1$ 是无截断的精确判据；把它当"答案"。
- FDTD 的"信号"：中心波矢落在动量带隙内的波包，其场峰值 $\max_x|E|$ 应随时间**指数增长**，增长速率 $\gamma=|\mathrm{Im}\,\omega_F|$（对应带隙中心）。半对数坐标下指数增长是一条直线，直线斜率应和 $\gamma$ 量级一致。
- 用 `boundary='periodic'`：波包永远留在域内，不会撞上吸收层；但也正因如此，任何带隙内的数值噪声也会被放大（[02 章 §4.3](02-physics-background.md#43-动量带隙场随时间指数增长衰减)），所以跑一段时间看到增长即可，不必跑满。

**起点代码**

```matlab
%% E8：三个引擎互相验证动量带隙（沿用 E7 的设计）
epsA = 1.0; epsB = 4.0; dutyA = 0.5;     % 与 E7 一致
T = 1; Omega = 2*pi/T;
kNorm = linspace(0.2, 1.4, 61);
reference = temporal_crystal_bands(kNorm*Omega, [epsA epsB], [1 1], ...
    [dutyA 1-dutyA]*T);

% --- E8a: TMM 精确定位每个连续带隙区间 ---
inGap = abs(reference.halfTrace) > 1;
d = diff([false, inGap, false]);         % 找连续区间的起止
starts = find(d == 1);  ends = find(d == -1) - 1;
kCenter = [];
for q = 1:numel(starts)
    fprintf('gap %d: kc/Omega in [%.3f, %.3f]\n', q, ...
        kNorm(starts(q)), kNorm(ends(q)));
    kCenter(q) = mean([kNorm(starts(q)), kNorm(ends(q))]);
end
kGapCenter = kCenter(1);                 % 选第一个带隙的中心

% --- E8b: FDTD 打波包，观察带隙内的指数增长 ---
dx = 0.05;  x = 0:dx:300;                % 周期边界，域内 6001 格点
dt = 0.4*dx;                             % CFL = 0.4，安全
k0 = kGapCenter*Omega;                   % 带隙中心波矢
x0 = 100; sigma = 10;
profile = @(xq) exp(-((xq-x0)/sigma).^2).*exp(1i*k0*(xq-x0));
E0 = profile(x);
xH = x + dx/2;                           % 周期边界：H 网格与 E 等长
nEff = sqrt(dutyA*epsA + (1-dutyA)*epsB);
Hhalf0 = nEff*profile(xH + dt/2);

cfg.x = x; cfg.dt = dt; cfg.nSteps = 3000;
cfg.epsFun = @(xq,tq) epsA + (epsB-epsA)*(mod(tq,T) >= dutyA*T);  % 时间方波
cfg.muFun  = @(xq,tq) ones(size(xq));
cfg.E0 = E0; cfg.Hhalf0 = Hhalf0;
cfg.boundary = 'periodic';
cfg.recordEvery = 20;
out = fdtd1d_db(cfg);

peakE = max(abs(out.E),[],2);
figure; semilogy(out.t, peakE, 'LineWidth',1.4); grid on;
xlabel('t'); ylabel('max_x |E|');
```

**数值化的成功判据**

- E8a 打印的带隙区间，与 [E7](#e7--设计一个光子时间晶体并做收敛审计) 能带图里红色 `Im(omegaF)` 非零的区域**重合**——TMM 与 ST-PWE 一致。
- E8b 的半对数图上 `max|E|` 的后半段呈**直线上升**（指数增长），直线斜率换算出的增长率与 `Im(omegaF)/Omega` **量级一致**。把波包中心 `k0` 移到带隙外再跑一次，应**不增长**——这是干净的对照组。
- 若场增长太快导致数值溢出（`Inf`/`NaN`），把 `cfg.nSteps` 减小（比如 1500）再试；你只需要看到若干倍的指数增长就够了，不必跑满。

### E9 · 计算 Zak 相位，改变元胞原点看它如何变（仿 demo07）

**目标**

[demo07](../demos/demo07_zak_phase_stpwe.m) 对 Fig.2 介质的能带算空间 Zak 相位，并特别注明：**Zak 相位依赖于元胞原点的选取**。你的任务是把静态参照介质整体平移 `x0`（等价于移动元胞原点），观察 Zak 相位如何变化，并解释为什么。

**提示**

- 把介质平移 `x0`，就是把它的 Fourier 系数乘上相位因子：$\varepsilon_{m,n}\to\varepsilon_{m,n}e^{i n g x_0}$（想想为什么：平移 $x\to x-x_0$ 给 $e^{i n g x}$ 一个相位 $e^{i n g x_0}$）。
- 物理背景：平移原点会改变 Bloch 波函数的规范（乘上一个 $k$ 依赖的相位），这给每条 Wilson link 都乘上 $e^{-i\Delta k\,x_0}$，于是 Zak 相位整体平移 $\approx -g\,x_0 = -2\pi(x_0/\Lambda)$（模 $2\pi$）。**平移一个整周期 $\Lambda$ 后应回到原值**。
- demo07 的注释提醒：只有当元胞原点落在反演中心时，Zak 相位才取 0 或 $\pi$ 这两个量化值；Fig.2 元胞从层边界开始，所以不必是 0 或 $\pi$。跑出来看 `info.minimumLinkMagnitude`，太小（<0.1）就说明单带 Zak 不可靠。

**起点代码**

```matlab
%% E9：元胞原点平移如何改变 Zak 相位（仿 demo07 的静态参照）
p = stm_fig2_parameters();
Nk = 161;
kGrid = -p.g/2 + (0:Nk-1)*p.g/Nk;       % 不含右端点，闭合 BZ
omegaSeed = 0.344*p.g*p.c0;

pStatic = p; pStatic.modDepth = 0;      % 静态参照介质（空间二元）
for x0 = [0 0.25 0.5 0.75 1.0]          % 原点平移（单位：Lambda）
    epsCoeff0 = @(n) stm_fig2_eps_coeff(0,n,pStatic).*exp(1i*n*p.g*x0);
    sysS = stpwe_build_system(@(m,n) double(m==0)*epsCoeff0(n), ...
        [], 12, 0, p.g, p.Omega);
    opt.omegaWindow = [0.25 0.70]*p.g*p.c0;
    opt.maxImag = 1e-7*p.g*p.c0;
    opt.minM0Weight = 0;
    bandS = stpwe_track_band(sysS, kGrid, omegaSeed, opt);
    [zak, info] = zak_phase_biorthogonal(bandS.R, bandS.L, sysS.Bomega, ...
        stpwe_bz_sewing_matrix(sysS));
    fprintf('x0=%.2f: Zak/pi = %+.4f  minLink=%.2e\n', ...
        x0, zak/pi, info.minimumLinkMagnitude);
end
```

> 💡 这里用 `stpwe_track_band` 沿 $k$ 追踪同一条带（双正交重叠判连续），再用 `zak_phase_biorthogonal` 做闭合 Wilson loop。`stpwe_bz_sewing_matrix(sysS)` 负责把 $k=-g/2$ 与 $k=+g/2$ 的等价态缝起来。这一整套就是 [07 拓扑] 章的内容。

**数值化的成功判据**

- 打印的 `Zak/pi` 随 `x0` 近似**线性平移**：每平移 $0.25\Lambda$，`Zak/pi` 大约改变 $0.5$（以 2 为模，因为最终要 `wrap_to_pi` 到 $(-\pi,\pi]$）；`x0=1.0` 时回到 `x0=0` 的值。
- 你能用一句话解释：**Zak 相位不是任意原点下的不变量**——它随原点平移 $-\mathrm{g}\,x_0$（模 $2\pi$），只有在原点落在反演中心（或平移整数个晶格）时才取量化值。这正好呼应 demo07 注释里的警告。
- 追加思考：对照 [demo08](../demos/demo08_fhs_chern_thouless_pump.m) 的 Chern 数——为什么 Chern 数（二维环面上的积分）没有这种原点问题？提示：二维的"闭环"让规范相位完全抵消，一维的开区间端点则留下边界项。

---

## 五、成功判据汇总表

| 练习 | 数值化的成功判据（一句话） | 主要工具/函数 | 对应章节 |
|------|--------------------------|--------------|---------|
| E1 | `Nspace>=10` 后 `Re(omega)` 曲线几乎不变；`Mtime` 从 1 加到 2 结果不变 | `stpwe_build_system`, `stpwe_solve_omega` | 04 ST-PWE、08 demo |
| E2 | 改 `nBefore`/`nAfter` 后数值 `tau`/`rho` 相对解析公式误差 < 1e-2 | `fdtd1d_db`, `temporal_interface_matrix`, `out.E`/`out.H`/`out.t` | 06 FDTD、02 时间界面 |
| E3 | 带隙中心/宽度随 contrast、duty 的变化与"$\bar n/2$"估计一致，并能解释 | `temporal_crystal_bands`, `bands.halfTrace` | 05 TMM、02 动量带隙 |
| E4 | 静态空间光栅 `out.energy(end)/out.energy(1)` 偏离 1 < 5%；时间调制下能解释能量为何不守恒 | `fdtd1d_db`, `out.energy` | 06 FDTD、conventions 第 6 条 |
| E5 | 共轭关系自检为 0；`Re(omega(+k))=Re(omega(-k))`（~1e-9） | `stpwe_build_system`, `stpwe_solve_omega` | 04 ST-PWE、conventions 第 5 节 |
| E6 | DB 与 EB 的 `tau`/`rho` 符合预期值，并能给出物理解释 | `temporal_interface_matrix_jump` | 05 时间界面、02 |
| E7 | 收敛审计误差随 `M` 单调下降，`M=19` 时 max 误差 < ~1e-2 | `stpwe_build_system`, `stpwe_fold_frequency`, `temporal_crystal_bands` | 04、05、08 demo12 |
| E8 | TMM 带隙区间与 ST-PWE 图一致；FDTD 带隙内波包指数增长、带隙外不增长 | `temporal_crystal_bands`, `fdtd1d_db` | 04、05、06 |
| E9 | `Zak/pi` 随原点平移约 $-2(x_0/\Lambda)$（模 2），`x0=1` 回到原值，能解释原因 | `stpwe_track_band`, `zak_phase_biorthogonal`, `stpwe_bz_sewing_matrix` | 07 拓扑、08 demo07 |

---

## 六、解答提示（简短，不贴完整答案）

**E1**：正弦调制只有 $m=0,\pm1$，所以 `Mtime=1` 时间维已精确；要收敛测试的是 `Nspace`。`Nspace=10` 与 `20` 曲线几乎重合即收敛。

**E2**：相对误差应稳定在 1e-2 以下。若反差过大误差变大，检查波包是否在 `tEnd` 前撞到吸收层（把 `x0` 提前或 `tEnd` 调小）。分解时记住必须用新的 `nAfter`。

**E3**：带隙中心约在 $\bar k_c\approx\bar n/2$（$\bar n=\sqrt{\mathrm{dutyA}\,\varepsilon_A+(1-\mathrm{dutyA})\varepsilon_B}$）。`epsB` 增大 → 中心右移、宽度变宽；`dutyA` 从 0.5 改走，偶次谐波出现，带隙结构定性改变。

**E4**：无时间调制 → 能量守恒，判据 <5%；加时间调制 → 调制做功，能量振荡/净变，属正常（conventions 第 6 条）。

**E5**：实函数 ⇒ 共轭关系 ⇒ 若 $\omega$ 是 $+k$ 本征值则 $\mathrm{conj}(\omega)$ 是 $-k$ 本征值 ⇒ `Re(omega)` 在 $\pm k$ 相同。行波调制是最简单的实时空介质。

**E6**：DB 下 $r_E=\varepsilon_b/\varepsilon_a<1$（E 随 ε 跳变、变弱），EB 下 $r_E=1$（`jumpD` 补偿了 E 的跳变），所以 EB 的 `tau`/`rho` 更大；`rho` 的符号差对应时间反射波的相位。

**E7**：方波系数 $\sim1/|m|$ 衰减 → 收敛慢；以 TMM 为参照，误差随 `M` 单调下降，带隙边缘最难收敛。`M=19` 时 max 误差 ~1e-2 量级即"够"。

**E8**：带隙内的半对数 `max|E|` 是直线，斜率对应 `Im(omegaF)`；把 `k0` 移出带隙作对照组应不增长。增长太快就减 `nSteps`。

**E9**：原点平移给 Bloch 波函数乘 $k$ 依赖相位 → 每条 Wilson link 乘 $e^{-i\Delta k\,x_0}$ → Zak 平移 $-g\,x_0$（模 $2\pi$）；平移整周期回到原值；原点落反演中心才取 0/π（呼应 demo07 注释）。

---

## 下一步

- 想把自己的设计变成可发表的图？先跑一遍 [08 章的 `run_all_demos`](../README.md)，再对照 [docs/validation.md](../docs/validation.md) 看"可信度"要补哪些验证。
- 想找研究方向？读 [docs/roadmap.md](../docs/roadmap.md)。动手之前，用 [03 章 §二](03-methods-and-map.md#二选型决策树) 检查你的模型属于三类中的哪一类。
- 遇到报错卡住了？翻 [10-faq-troubleshooting.md](10-faq-troubleshooting.md)（随用随查），或回到 [01-matlab-primer.md](01-matlab-primer.md) 的调试一节。
- 做完所有练习后，跑一遍 `test_smoke`，看你能不能读懂每个断言在验证什么——那就是你以后自己写"判据"的模板。

---

*练习维护提示：所有练习中的函数名、字段名、参数顺序都以 [docs/tool-reference.md](../docs/tool-reference.md) 的函数签名为准；改过任何函数签名后，请回来同步本章的起点代码。*
