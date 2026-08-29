# V3（3.0.0）科研模板重大重构计划

## 目标与交付物

- 在仓库根目录保留本文件，记录完整重构任务。
- 本次只创建计划文档，不执行目录删除、代码重写、版本更新或 Git 标签操作。
- V3 定位为破坏性更新：不保留旧接口、兼容入口、样例、测试或路径初始化。
- 新版统一采用一维、标量、无色散、归一化单位 $c_0=1$，默认 $\mu_r=1$。
- 所有 MATLAB 文件平铺在根目录，无需 `addpath`、`startup` 或切换目录。

目标结构：

```text
README.md
VERSION.txt
.gitignore
V3_REFACTOR_PLAN.md
run_pwe.m
pwe_fourier.m
pwe_bands.m
run_tmm.m
tmm_bands.m
run_fdtd.m
fdtd1d.m
fdtd_fft_bands.m
```

## 核心重构任务

### 1. 清理旧版

- 重构前建立 `legacy-pre-v3-refactor` Git 标签，确保旧代码和论文资料可恢复。
- 从新版工作树移除 `core/`、`tmm/`、`fdtd/`、`topology/`、`demos/`、`examples/`、`tests/`、`learn/`、`docs/`、`reproduction/`、`Simple_FDTD_NathanZechar/` 和全部生成结果。
- 删除 `startup_stm.m`、`stm_init.m`、所有批量运行器、动画、缓存、测试媒体及旧接口文档。
- 不建立 `archive/`，避免新版目录再次膨胀。

### 2. PWE：一个入口、两个函数

- `run_pwe.m` 固定按“定义 $\varepsilon_r(x,t)$ → 周期和截断参数 → 计算 → 作图”组织。
- `pwe_fourier.m` 在无重复端点的时空原胞上采样材料，按

  $$
  \varepsilon(x,t)=\sum_{m,n}\varepsilon_{mn}e^{i n g x-i m\Omega t}
  $$

  生成覆盖 $[-2M,2M]\times[-2N,2N]$ 的 Fourier 系数和卷积数据。
- `pwe_bands.m` 统一支持固定 $k\rightarrow\omega$ 和固定 $\omega\rightarrow k$，保留原始复谱、第一 Floquet 布里渊区折叠结果、中心谐波权重及可选本征矢。
- 删除材料预设、静态参考带、模式选择、PWE 场重构及拓扑调用；场分布统一交给 FDTD。

公开调用：

```matlab
fourier = pwe_fourier(epsFun, muFun, pweCfg);
pweOmega = pwe_bands(fourier, pweCfg, 'omega');
pweK = pwe_bands(fourier, pweCfg, 'k');
```

### 3. TMM：一个入口、一个函数

- `run_tmm.m` 定义一个时间周期内的 `epsLayers`、`muLayers`、`durations` 和波矢扫描，随后计算并绘制复准频率与带隙。
- `tmm_bands.m` 使用连续状态 $[D;B]$ 构造单周期矩阵，计算 Floquet 乘子、复准频率、半迹和带隙标记。
- 明确 TMM 仅适用于空间均匀、分段常数、瞬时切换的 $\varepsilon(t)$，并采用时间界面上 $D/B$ 连续的模型。
- 不保留方向振幅多层链、有限周期散射、畴壁搜索和界面定制函数。

公开调用：

```matlab
tmmResult = tmm_bands(kValues, epsLayers, muLayers, durations);
```

输出包含 `U`、`lambda`、`omega`、`halfTrace`、`gapMask` 和完整输入参数。

### 4. FDTD 与 FFT：一个入口、两个函数

- `run_fdtd.m` 共用一个材料定义，但执行两次目的明确的计算：
  - 波包加吸收边界，用于场分布；
  - 宽带局域源加周期边界，用于 FFT 能带。
- `fdtd1d.m` 仅保留 D/B-Yee 推进、周期/海绵边界、源、记录间隔、CFL 检查和突变时间界面的中心修正。
- 删除探针工具、动画、缓存、频谱过滤、单精度模式、进度条和论文专用审计。
- `fdtd_fft_bands.m` 对完整 $D(x,t)$ 做空间 FFT、时间 FFT、窗函数、功率归一化以及 Floquet 折叠：
  - 有空间周期时同时折叠 $k$ 和 $\omega$；
  - 空间均匀时间晶体可关闭 $k$ 折叠；
  - 输出谱强度、频率/波矢网格、脊线和真实 FFT 分辨率。
- FFT 函数不得暗中调用 TMM/PWE；理论曲线只能由入口脚本显式叠加。
- 文档明确 FFT 直接给出实频谱脊线，不能单独恢复带隙中的 $\operatorname{Im}\omega$。

公开调用：

```matlab
fieldResult = fdtd1d(fieldCfg);
bandField = fdtd1d(bandCfg);
fftResult = fdtd_fft_bands( ...
    bandField.D, bandField.x, bandField.t, fftCfg);
```

## 文档、版本与科研使用方式

- 将 `VERSION.txt` 更新为 `3.0.0`。
- 重写 `README.md`，只保留三个执行入口、统一归一化方式、相位约定、方法适用范围和最短运行说明。
- 三个入口脚本均在顶部集中放置材料和计算参数，计算结果保留在工作区，默认不创建缓存或输出目录。
- 图像和 MAT 数据仅通过入口末尾清晰可见的可选保存块导出。
- 最低环境定为 MATLAB R2020a、仅 Base MATLAB，不支持 2D/3D、各向异性或色散材料。

## 实施顺序

1. 固化旧版 Git 标签，记录重构前提交。
2. 编写并验证 PWE 三个文件。
3. 编写并验证 TMM 两个文件。
4. 编写并验证 FDTD/FFT 三个文件。
5. 用新入口完成三种方法的交叉验收。
6. 重写 `README.md`、更新 `VERSION.txt` 和 `.gitignore`。
7. 删除旧目录、旧入口及生成文件，确认目标结构。

## 一次性验收

验收只在重构实施期间执行，不在新版保留 `tests/`、冒烟测试或收敛性测试文件。

- 三个入口从仓库根目录直接运行，且不存在 `addpath` 或旧函数依赖。
- 均匀介质下 PWE/TMM 恢复解析色散。
- 二元时间晶体中 PWE 与 TMM 的折叠带位置、增长/衰减支和带隙位置一致。
- FDTD 的 CFL 严格小于 1，突变时刻与时间网格对齐，D/B 连续修正生效。
- 时间晶体基准中，FDTD-FFT 与 TMM 通带脊线满足中位误差 $|\Delta\omega T|\leq 0.10$、最大误差不超过 `0.30`。
- 一般 $\varepsilon(x,t)$ 算例能够输出第一时空布里渊区内的双重折叠谱。
- 最终确认工作树仅包含目标文件，旧资料可从 `legacy-pre-v3-refactor` 恢复。

## 实施约束

- PWE Fourier 采样必须覆盖卷积矩阵所需的差阶，不能对超出采样范围的系数静默补零。
- PWE 必须保留未筛选原始本征谱；绘图筛选不得覆盖原始结果。
- TMM 不计算高次单周期矩阵幂，避免强增益情况下的溢出和病态问题。
- FDTD 始终推进 D/B；不得改成会破坏理想时间界面连续条件的直接 E/H 跃迁。
- 突变时间必须落在整数时间节点，界面更新保留两侧本构关系的中心平均。
- FFT 使用整数个调制周期和无重复端点的周期网格；零填充只细化绘图网格，不宣称提高真实分辨率。
- 海绵边界需明确标注为简化吸收层而非 PML。
- 不添加隐藏路径操作、自动缓存、隐式理论曲线或论文专用参数。
