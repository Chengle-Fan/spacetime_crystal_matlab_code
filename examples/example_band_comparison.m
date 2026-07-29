function example_band_comparison(quality)
%EXAMPLE_BAND_COMPARISON 分别绘制三种能带结构，便于逐幅比对。
%
%   example_band_comparison('quick')
%   example_band_comparison('paper')
%
% 在同一参数配置下计算并分三张图绘制：
%   图 1 — 静态光子能带（无时间调制，仅空间周期性）
%          对应原文图 (b) 中的黑色实线
%   图 2 — 空晶格近似（静态能带 + umklapp 复制）
%          对应原文图 (b) 中的黑色实线 + 灰色虚线
%   图 3 — 完整 Floquet 能带（ST-PWE 数值解）
%          对应原文图 (b) 中的彩色散点
%
% 三张图并排打开，便于对比时间调制对能带结构的影响：
%   - 图 1 → 图 2：umklapp 折叠如何引入能带交叉
%   - 图 2 → 图 3：时间调制如何打开带隙（交叉处劈裂）

% =========================================================================
% 参数设置（与 example_floquet_bands_and_fields 共用配置）
% =========================================================================

if nargin < 1
    quality = 'quick';
end
rootDir = stm_init();
p = stm_preset_modulated_slab();
% p.Lambda  = 1           — 归一化空间周期
% p.c0      = 1           — 归一化光速
% p.g       = 2π          — 倒格矢
% p.Omega   = 2π*0.20     — 调制角频率
% p.OmegaBar = 0.20       — 归一化调制频率 ω_mod/(g*c0)

switch lower(quality)
    case 'paper'
        Nspace = 20;
        Nk = 181;
    case 'quick'
        Nspace = 10;
        Nk = 101;
    otherwise
        error('quality must be ''quick'' or ''paper''.');
end
Mtime = 1;

% --- 构造时空系统矩阵 ---
epsCoeff = @(m,n) stm_fourier_modulated_slab(m,n,p);
sys = stpwe_build_system(epsCoeff, [], Nspace, Mtime, p.g, p.Omega);

% =========================================================================
% 公共数据：k 网格 + Floquet 本征解
% =========================================================================

kBar = linspace(-0.5, 0.5, Nk);    % 归一化波矢 kΛ/(2π)
kValues = p.g*kBar;                 % 物理波矢
fMax = 0.82;                        % 频率显示上限
weightThreshold = 0.035;            % m=0 参与度阈值
imagTolerance = 3e-3;               % 虚部容差

% --- 扫 k 求解 Floquet 能带 ---
allK = cell(Nk,1);
allF = cell(Nk,1);
allWeight = cell(Nk,1);

fprintf('Computing Floquet bands: %d k-points, matrix %d x %d.\n', ...
    Nk, 2*sys.S, 2*sys.S);
for ik = 1:Nk
    sol = stpwe_solve_omega(sys, kValues(ik));
    fBar = sol.omega/(p.g*p.c0);
    keep = isfinite(fBar) ...
        & real(fBar) >= 0 & real(fBar) <= fMax ...
        & abs(imag(fBar)) <= imagTolerance ...
        & sol.m0Weight >= weightThreshold;
    allK{ik} = repmat(kBar(ik), sum(keep), 1);
    allF{ik} = real(fBar(keep));
    allWeight{ik} = sol.m0Weight(keep);
end
allK = vertcat(allK{:});
allF = vertcat(allF{:});
allWeight = vertcat(allWeight{:});

% --- 静态参考能带（无时间调制） ---
epsCoeff0 = @(n) stm_fourier_modulated_slab(0,n,p);
nStaticBands = 8;
fStatic = stpwe_static_bands(kValues, epsCoeff0, Nspace, ...
    p.g, p.c0, nStaticBands);

% =========================================================================
% 图 1：静态光子能带（仅空间周期性，无时间调制）
% =========================================================================

fig1 = figure('Color','w','Position',[60 400 620 480]);
ax1 = axes(fig1);
hold(ax1, 'on');

for ib = 1:size(fStatic,1)
    plot(ax1, kBar, fStatic(ib,:), 'Color',[0.00 0.45 0.74], ...
        'LineWidth',1.2);
end

xlim(ax1, [-0.5 0.5]);
ylim(ax1, [0 fMax]);
xlabel(ax1, 'k\Lambda/(2\pi)');
ylabel(ax1, '\omega\Lambda/(2\pi c)');
title(ax1, ...
    ['(a) Static photonic bands  (no time modulation, ' ...
     num2str(nStaticBands) ' bands)']);
grid(ax1, 'on');
box(ax1, 'on');

% 标注带折叠参考线 Ω 的位置，方便与图 2 对比
yline(ax1, p.OmegaBar, '--', 'Color',[0.85 0.33 0.10], ...
    'LineWidth',1.0);
text(ax1, -0.48, p.OmegaBar - 0.02, '\Omega/(gc)', ...
    'Color',[0.85 0.33 0.10], 'FontSize',9, 'VerticalAlignment','top');

% =========================================================================
% 图 2：空晶格近似（静态能带 + umklapp 折叠复制）
% =========================================================================

fig2 = figure('Color','w','Position',[690 400 620 480]);
ax2 = axes(fig2);
hold(ax2, 'on');

% 灰色虚线：umklapp 复制（平移 ±Ω, ±2Ω）
umklappShifts = [-2 -1 1 2];
for ib = 1:size(fStatic,1)
    for shift = umklappShifts
        plot(ax2, kBar, fStatic(ib,:) + shift*p.OmegaBar, '--', ...
            'Color',[0.60 0.60 0.60], 'LineWidth',0.45);
    end
end

% 黑色实线：中心能带（shift=0）
for ib = 1:size(fStatic,1)
    plot(ax2, kBar, fStatic(ib,:), 'Color',[0.00 0.00 0.00], ...
        'LineWidth',1.0);
end

% 浅灰色水平参考线：每个 umklapp 阶的折叠中心
for shift = [-1 0 1]
    yline(ax2, shift * p.OmegaBar, ':', ...
        'Color',[0.70 0.70 0.70], 'LineWidth',0.5);
end

xlim(ax2, [-0.5 0.5]);
ylim(ax2, [0 fMax]);
xlabel(ax2, 'k\Lambda/(2\pi)');
ylabel(ax2, '\omega\Lambda/(2\pi c)');
title(ax2, ...
    '(b) Empty-lattice approx.  (static bands + umklapp replicas \pmn\Omega, \pm2\Omega)');
legend(ax2, {'umklapp \pm\Omega,\pm2\Omega', 'static (m=0)'}, ...
    'Location','northeast', 'FontSize',8);
grid(ax2, 'on');
box(ax2, 'on');

% =========================================================================
% 图 3：完整 Floquet 能带（ST-PWE 数值解，含时间调制）
% =========================================================================

fig3 = figure('Color','w','Position',[1320 400 620 480]);
ax3 = axes(fig3);
hold(ax3, 'on');

% --- 背景：静态能带（浅灰细线，便于对比 Floquet 折叠偏移） ---
for ib = 1:size(fStatic,1)
    plot(ax3, kBar, fStatic(ib,:), 'Color',[0.80 0.80 0.80], ...
        'LineWidth',0.5);
    for shift = umklappShifts
        plot(ax3, kBar, fStatic(ib,:) + shift*p.OmegaBar, ':', ...
            'Color',[0.88 0.88 0.88], 'LineWidth',0.3);
    end
end

% --- 彩色散点：Floquet 本征模，颜色 = m=0 扇区参与度 ---
scatter(ax3, allK, allF, 12, allWeight, 'filled');

% --- 注释：m=0 参与度的物理意义---
% 品红色（低参与度）→ 模式含有显著的高阶时间谐波分量 → 强调制效应
% 蓝色  （高参与度）→ 模式主要成分在 m=0 扇区 → 接近静态行为

xlim(ax3, [-0.5 0.5]);
ylim(ax3, [0 fMax]);
xlabel(ax3, 'k\Lambda/(2\pi)');
ylabel(ax3, '\omega\Lambda/(2\pi c)');
title(ax3, ...
    '(c) Full Floquet bands  (ST-PWE, colored by m=0 participation)');
grid(ax3, 'on');
box(ax3, 'on');

% 自定义色彩映射：品红(高 Floquet 混合) → 蓝(纯静态)
cmapFloquet = [linspace(1,0,256).', zeros(256,1), ones(256,1)];
colormap(ax3, cmapFloquet);
caxis(ax3, [0 1]);
cb3 = colorbar(ax3);
cb3.Label.String = 'm=0 sector participation';

% =========================================================================
% 输出
% =========================================================================

fprintf('\n=== Band comparison summary ===\n');
fprintf('Figure 1: Static bands (no modulation)  — %d bands\n', nStaticBands);
fprintf('Figure 2: Empty-lattice approximation   — static + umklapp copies\n');
fprintf('Figure 3: Full Floquet bands            — %d eigenmodes plotted\n', ...
    numel(allF));

% 保存各图
tags = {'static_bands', 'empty_lattice', 'floquet_bands'};
figs = [fig1, fig2, fig3];
for i = 1:3
    outputFile = fullfile(rootDir, 'output', ...
        ['example_band_comparison_' tags{i} '_' lower(quality) '.png']);
    save_example_figure(figs(i), outputFile);
    fprintf('Saved %s\n', outputFile);
end

end

% =========================================================================
% 辅助函数：保存图片
% =========================================================================

function save_example_figure(fig, outputFile)
try
    exportgraphics(fig, outputFile, 'Resolution',220);
catch
    print(fig, outputFile, '-dpng', '-r220');
end
end
