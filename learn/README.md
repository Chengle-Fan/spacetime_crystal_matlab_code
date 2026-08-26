# Learn：从零上手时空介质工具箱

> 面向**有物理基础、但缺乏科学计算与编程经验**的读者。假设你懂麦克斯韦方程、能带、Fourier 级数这些物理数学概念，但**没怎么用过 MATLAB**，也**没写过仿真代码**。
> 本系列教程会把本项目（`spacetime_crystal_code`，一个用于研究**光子时间晶体 / 时空晶体**的 MATLAB 工具箱）从"能跑起来"讲到"能自己动手做研究"。

---

## 一、这个项目是干什么的？

一句话：**用数值方法研究"介电常数既随空间、又随时间变化"的材料中电磁波的行为。**

- 如果介电常数只在**空间**上周期变化，就是普通的光子晶体（能带、带隙）；
- 如果介电常数只在**时间**上周期变化，就是**光子时间晶体（PTC）**，它有不寻常的"动量带隙"（场在时间上指数放大/衰减）；
- 如果**空间和时间同时**周期变化，就是**时空晶体（spacetime crystal）**，本项目的核心研究对象。

为了研究这些介质，工具箱提供了三套数值"引擎"，它们各有所长、互相验证：

| 引擎 | 目录 | 一句话原理 | 擅长什么 |
|------|------|-----------|---------|
| **ST-PWE**（时空平面波展开） | `core/` | 把场和介质都展开成 Fourier 级数，把 Maxwell 方程变成**本征值问题** | 计算能带、复频率/复波数、场型 |
| **时间 TMM**（时间传输矩阵） | `tmm/` | 把时间分成一段段"台阶"，用 2×2 矩阵精确演化和/或匹配 | 时间界面、时间多层、精确能带、畴壁态 |
| **FDTD**（时域有限差分） | `fdtd/` | 在空间和时间网格上直接"蛙跳"推进 Maxwell 方程 | 波包演化、时域验证、任意时空变化材料 |

再加上 `topology/` 目录的**拓扑不变量**（Zak 相位、Chern 数），以及 12 个**论文复现 demo** 和两个**论文复现工程**。

> 这三套引擎为什么要并存、各自在什么时候用，详见 [03-methods-and-map.md](03-methods-and-map.md)。

---

## 二、学习路线图

整套教程按"**先会跑 → 再懂物理 → 再懂代码 → 再动手**"组织。建议**从上到下、一章一章**来。

| 阶段 | 教程文件 | 学完你能做到 | 预计时间 |
|------|---------|-------------|---------|
| **0. 准备** | [01-matlab-primer.md](01-matlab-primer.md) MATLAB 与科学计算入门 | 看懂 MATLAB 语法、结构体、函数句柄、画图；把工具箱路径设好；跑通自检测试 | 2–4 小时 |
| **1. 物理** | [02-physics-background.md](02-physics-background.md) 时变介质物理背景 | 用直觉理解时间界面、时间反射/折射、PTC、时空晶体、动量带隙 | 1–2 小时 |
| **2. 方法** | [03-methods-and-map.md](03-methods-and-map.md) 三种数值方法与工具箱地图 | 明白三个引擎各自怎么选、代码放在哪、v1/v2 命名怎么分辨 | 1 小时 |
| **3. 引擎 1** | [04-st-pwe.md](04-st-pwe.md) ST-PWE 精讲 | 从 Maxwell 方程一路看懂 `stpwe_build_system` → `stpwe_solve_omega` 的完整管线 | 2–3 小时 |
| **4. 引擎 2** | [05-temporal-tmm.md](05-temporal-tmm.md) 时间传输矩阵精讲 | 看懂时间界面匹配、单值矩阵、能带、多层、畴壁 | 2–3 小时 |
| **5. 引擎 3** | [06-fdtd.md](06-fdtd.md) FDTD 时域仿真精讲 | 看懂 `fdtd1d_db` 的每个配置字段，会自己搭一个 FDTD 仿真 | 2–3 小时 |
| **6. 拓扑** | [07-topology.md](07-topology.md) 拓扑不变量精讲 | 看懂 Zak 相位、Chern 数在代码里怎么算 | 1–2 小时 |
| **7. 实战** | [08-demos-walkthrough.md](08-demos-walkthrough.md) 演示脚本导读 | 12 个论文复现 demo 逐个能跑、能读懂图 | 2–3 小时 |
| **8. 练习** | [09-exercises.md](09-exercises.md) 动手练习 | 独立完成改参数→改材料→设计新晶体的系列练习 | 3–6 小时 |
| **9. 排错** | [10-faq-troubleshooting.md](10-faq-troubleshooting.md) 常见问题与排错 | 遇到报错不慌，能自己定位 | 随用随查 |

> ⏱️ 时间只是粗略估计，取决于你对 MATLAB 的熟悉程度。**所有章节都建议实际动手运行代码**，只看不跑效果会大打折扣。

---

## 三、与 docs/ 参考文档的分工

`docs/` 目录剩下的 4 个文件是**参考书**（随查随用），本 `learn/` 系列是**教材**（从头讲到尾）：

- **物理约定与完整用法示例的权威来源现在是本系列各章**：时间因子、归一化单位、界面连续量、模型选择边界见 [02 章](02-physics-background.md)（物理背景）与 [03 章](03-methods-and-map.md)（选型决策树）；三个引擎的完整用法见 [04](04-st-pwe.md)、[05](05-temporal-tmm.md)、[06](06-fdtd.md) 章；拓扑量的区分见 [07 章](07-topology.md)。原来单独的 conventions / manual 文档已并入本系列，不再单列。
- `docs/` 剩下的文件负责本系列**没有复制**的独特内容：

| 现有文档 | 它的定位 | 什么时候用 |
|---------|---------|-----------|
| [docs/tool-reference.md](../docs/tool-reference.md) | **函数字典**：全部 113 个函数的输入/输出/依赖（只给签名，无用法示例） | 遇到陌生函数，先回对应章节的"代码导读"，再查这里对签名 |
| [docs/paper-map.md](../docs/paper-map.md) | 文献覆盖矩阵 | 想知道"这套代码能复现哪些论文"时 |
| [docs/roadmap.md](../docs/roadmap.md) | 研究路线图 | 学完后想找研究方向时 |
| [docs/validation.md](../docs/validation.md) | 验证记录 | 想确认"这套代码可信吗"时 |

> 用法示例回本系列对应章节，`docs/tool-reference.md` 只给函数签名；物理约定以本系列各章为准。

> 阅读顺序建议：本 `learn/` 系列是"教材"（从头讲到尾），`docs/` 是"参考书"（随查随用）。两者用相对路径互相链接。

---

## 四、项目目录速览

```
spacetime_crystal_code/
├── README.md                 # 项目主页（英文+中文摘要）
├── startup_stm.m             # v2 入口：把 core/tmm/fdtd/topology/demos/tests 加入路径
├── stm_init.m                # v1 入口：把 core/tmm/fdtd/topology/examples/tests 加入路径
├── run_all_demos.m           # 一键运行 12 个 demos（false=快速，true=含重计算）
├── stm_run_examples.m        # 一键运行 v1 examples
├── core/                     # 引擎 1：ST-PWE（时空平面波展开）
├── tmm/                      # 引擎 2：时间传输矩阵
├── fdtd/                     # 引擎 3：一维 D/B-Yee FDTD
├── topology/                 # 拓扑不变量（Zak 相位、Chern 数）
├── demos/                    # 12 个论文复现演示脚本（v2，重点学习对象）
├── examples/                 # 10 个示例脚本（v1，旧版）
├── tests/                    # 自检测试（test_smoke 等）
├── docs/                     # 参考文档（见上表）
├── reproduction/             # 两个论文复现工程（进阶）
│   ├── Observation of temporal reflection/      # Nature Physics 2023 时间反射
│   └── Topology of photonic crystals/           # Optica 2018 PTC 拓扑
├── output/                   # 所有脚本保存图片/数据的目录（自动创建）
├── Simple_FDTD_NathanZechar/ # 第三方教学 FDTD（1D/2D/3D），可选学习
└── learn/                    # 本教程
```

> ⚠️ **两个入口的区别**（新手最容易踩的坑）：
> - `startup_stm` 把 **demos/** 加入路径（v2 演示脚本需要它）；
> - `stm_init` 把 **examples/** 加入路径（v1 示例脚本需要它）。
> 想跑哪个系列就用对应的入口。详见 [03-methods-and-map.md](03-methods-and-map.md#五v1-vs-v2同一物理量两套名字)。

---

## 五、快速开始（先跑起来再说）

在 MATLAB 中，把当前目录切到本项目根目录，然后：

```matlab
% 1) 设置路径（v2 系列，推荐）
startup_stm

% 2) 运行自检测试（会做十多个数值断言，全过会打印 "All smoke tests passed."）
test_smoke
```

如果看到 `All smoke tests passed.`，恭喜，环境已经就绪，可以正式开始学习。

> 没有 MATLAB？本项目是纯 MATLAB 代码，需要安装 MATLAB（建议 **R2020a 或更新**，因为 demo 用了 `tiledlayout` 等较新的绘图函数；旧版本多数脚本也能跑，部分画图函数有兼容回退）。
> 需要 MATLAB 许可证（没有的话 Octave 可部分兼容，但**不保证**所有功能可用）。

---

## 六、给不同基础读者的建议

- **MATLAB 完全新手**：先老老实实做 [01-matlab-primer.md](01-matlab-primer.md) 里的所有练习，再往下走。
- **MATLAB 会用、但没做过科学计算**：跳过 01 的大部分，重点看 01 的"科学计算思维"一节和 [03-methods-and-map.md](03-methods-and-map.md)。
- **想尽快看懂论文复现**：可以先跑 [08-demos-walkthrough.md](08-demos-walkthrough.md) 里的 demo，遇到不懂的函数再回头补前面的方法章节。
- **想做自己的研究**：学完 09 后，读 [docs/roadmap.md](../docs/roadmap.md) 找方向，用 [03 章 §二](03-methods-and-map.md#二选型决策树) 检查自己的模型属于三类中的哪一类。

---

*教程维护提示：本系列是物理约定的权威来源（时间因子 $e^{ikx-i\omega t}$、归一化单位 $\varepsilon_0=\mu_0=c_0=1$ 等），编写时请保持各章一致。*
