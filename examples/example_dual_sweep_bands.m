function example_dual_sweep_bands(quality)
%EXAMPLE_DUAL_SWEEP_BANDS 分别用"扫 k 求 ω"和"扫 ω 求 k"计算能带，并合并对比。
%
%   example_dual_sweep_bands('quick')
%   example_dual_sweep_bands('paper')
%
% 三张图：
%   图 1 — 固定 k 求 ω（ST-PWE 常规频带结构）
%           实线 = Re(ω)，虚线 = Im(ω)
%   图 2 — 固定 ω 求 k（ST-PWE 动量本征值）
%           实线 = Re(k)，虚线 = Im(k)
%   图 3 — 合并图（左：两方法的实部叠加；右：两方法的虚部叠加）
%
% 前两张图分别展示"哪里有频率禁带"和"哪里有动量禁带"；
% 第三张图将二者画在同一坐标系中，直观对比两种本征问题的互补性：
%   - 扫 k 能覆盖的区域（常规能带，无频率禁带）
%   - 扫 ω 能覆盖的区域（常规能带，无动量禁带）
%   - 二者不重合处即为禁带（一侧有本征解，另一侧无实数解）

% =========================================================================
% 参数设置
% =========================================================================

if nargin < 1
    quality = 'quick';
end
rootDir = stm_init();
p = stm_preset_modulated_slab();
% p.Lambda = 1, p.c0 = 1, p.g = 2π, p.OmegaBar = 0.20

switch lower(quality)
    case 'paper'
        Nspace = 20;
        Nk = 181;
        Nw = 161;
    case 'quick'
        Nspace = 10;
        Nk = 101;
        Nw = 81;
    otherwise
        error('quality must be ''quick'' or ''paper''.');
end
Mtime = 1;

% --- 构造系统矩阵 ---
epsCoeff = @(m,n) stm_fourier_modulated_slab(m,n,p);
sys = stpwe_build_system(epsCoeff, [], Nspace, Mtime, p.g, p.Omega);

% =========================================================================
% 公共参数
% =========================================================================

fMax = 0.82;                   % 频率显示上限（归一化 fBar）
kHalfRange = 0.5;              % kBar 显示半宽
weightThreshold = 0.035;       % m=0 参与度阈值（仅 ω 扫使用）
imagTolKscan = 5e-2;           % 扫 k 时容许的 Im(ω) 上限（需足够大以显示动量禁带）
imagTolWscan = 0.12;           % 扫 ω 时容许的 Im(k) 上限

% =========================================================================
% 扫 k 求 ω：固定 k，求解广义本征问题 A Φ = ω B Φ
% =========================================================================

kBarSweep = linspace(-kHalfRange, kHalfRange, Nk);

fprintf('=== Sweep 1: fixed-k, solve omega  (%d points) ===\n', Nk);
solKscan_k  = cell(Nk,1);
solKscan_fr = cell(Nk,1);      % Re(ω)/(g*c0)
solKscan_fi = cell(Nk,1);      % Im(ω)/(g*c0)

for ik = 1:Nk
    sol = stpwe_solve_omega(sys, p.g*kBarSweep(ik));
    fBar = sol.omega/(p.g*p.c0);
    keep = isfinite(fBar) ...
        & real(fBar) >= 0 & real(fBar) <= fMax ...
        & abs(imag(fBar)) <= imagTolKscan ...
        & sol.m0Weight >= weightThreshold;
    solKscan_k{ik}  = repmat(kBarSweep(ik), sum(keep), 1);
    solKscan_fr{ik} = real(fBar(keep));
    solKscan_fi{ik} = imag(fBar(keep));
end
kscan_k  = vertcat(solKscan_k{:});
kscan_fr = vertcat(solKscan_fr{:});
kscan_fi = vertcat(solKscan_fi{:});

% =========================================================================
% 扫 ω 求 k：固定 ω，求解本征问题 A(ω) Φ = k Φ
% =========================================================================

fBarSweep = linspace(0, fMax, Nw);

fprintf('=== Sweep 2: fixed-omega, solve k  (%d points) ===\n', Nw);
solWscan_f  = cell(Nw,1);      % 固定频率 ω/(g*c0)
solWscan_kr = cell(Nw,1);      % Re(k)/g
solWscan_ki = cell(Nw,1);      % Im(k)/g

for iw = 1:Nw
    sol = stpwe_solve_k(sys, p.g*p.c0*fBarSweep(iw));
    kNorm = sol.k/p.g;         % kBar = k/g
    keep = isfinite(kNorm) ...
        & real(kNorm) >= -kHalfRange & real(kNorm) <= kHalfRange ...
        & abs(imag(kNorm)) <= imagTolWscan;
    solWscan_f{iw}  = repmat(fBarSweep(iw), sum(keep), 1);
    solWscan_kr{iw} = real(kNorm(keep));
    solWscan_ki{iw} = imag(kNorm(keep));
end
wscan_f  = vertcat(solWscan_f{:});
wscan_kr = vertcat(solWscan_kr{:});
wscan_ki = vertcat(solWscan_ki{:});

% =========================================================================
% 图 1：固定 k → 求 ω（频带结构）
%   实心圆 = Re(ω)（传播模），空心圆 = Im(ω)（参量增长/衰减）
%   共用 y 轴，直接对比实部与虚部的量级
% =========================================================================

fig1 = figure('Color','w','Position',[60 400 680 520]);
ax1 = axes(fig1);
hold(ax1, 'on');

% 实部 Re(ω)：蓝色实心圆（传播模，构成常规能带）
scatter(ax1, kscan_k, kscan_fr, 10, 'filled', ...
    'MarkerFaceColor',[0.00 0.45 0.74], 'MarkerEdgeColor','none', ...
    'DisplayName','Re(\omega)  (propagating)');

% 虚部 Im(ω)：橙色空心圆（Im(ω) ≠ 0 ⇔ 动量禁带）
scatter(ax1, kscan_k, kscan_fi, 14, ...
    'MarkerEdgeColor',[0.85 0.33 0.10], 'LineWidth',1.0, ...
    'DisplayName','Im(\omega)  (momentum gap indicator)');

yline(ax1, 0, 'k-', 'LineWidth',0.8);
hold(ax1, 'off');

xlim(ax1, [-0.5 0.5]);
ylim(ax1, [min(kscan_fi)-0.005, fMax]);   % y 轴下延以容纳负虚部
xlabel(ax1, 'k\Lambda/(2\pi)');
ylabel(ax1, '\omega\Lambda/(2\pi c)');
title(ax1, ...
    'Figure 1: Fixed-k band structure — Re(\omega) [\bullet] and Im(\omega) [\circ]');
legend(ax1, 'Location','northeast', 'FontSize',8);
grid(ax1, 'on');
box(ax1, 'on');
fprintf('Figure 1: %d modes (fixed-k scan)\n', numel(kscan_fr));

% =========================================================================
% 图 2：固定 ω → 求 k（动量本征值）
%   实心圆 = Re(k)（传播模），空心圆 = Im(k)（凋落模）
%   共用 x 轴，直接对比实部与虚部的量级
% =========================================================================

fig2 = figure('Color','w','Position',[60 400 680 520]);
ax2 = axes(fig2);
hold(ax2, 'on');

% 实部 Re(k)：绿色实心圆（传播模，构成常规能带）
scatter(ax2, wscan_kr, wscan_f, 10, 'filled', ...
    'MarkerFaceColor',[0.00 0.62 0.45], 'MarkerEdgeColor','none', ...
    'DisplayName','Re(k)  (propagating)');

% 虚部 Im(k)：橙色空心圆（Im(k) ≠ 0 ⇔ 频率禁带）
scatter(ax2, wscan_ki, wscan_f, 14, ...
    'MarkerEdgeColor',[0.85 0.33 0.10], 'LineWidth',1.0, ...
    'DisplayName','Im(k)  (frequency gap indicator)');

xline(ax2, 0, 'k-', 'LineWidth',0.8);
hold(ax2, 'off');

xlim(ax2, [min(wscan_ki)-0.02, 0.5]);    % x 轴左延以容纳负虚部
ylim(ax2, [0 fMax]);
xlabel(ax2, 'k\Lambda/(2\pi)');
ylabel(ax2, '\omega\Lambda/(2\pi c)');
title(ax2, ...
    'Figure 2: Fixed-\omega band structure — Re(k) [\bullet] and Im(k) [\circ]');
legend(ax2, 'Location','northeast', 'FontSize',8);
grid(ax2, 'on');
box(ax2, 'on');
fprintf('Figure 2: %d modes (fixed-omega scan)\n', numel(wscan_f));

% =========================================================================
% 图 3：合并对比 —— 2×2 布局
%   (a) 左上：实部叠加（两种方法的 (k,ω) 解画在同一平面）
%   (b) 右上：Im(ω) vs k（扫 k 解，揭示动量禁带）
%   (c) 左下：Im(k) vs ω（扫 ω 解，揭示频率禁带）
%   (d) 右下：解说文字
%   注：(b)(c) 的横轴含义不同，不可混叠，故分开
% =========================================================================

fig3 = figure('Color','w','Position',[60 60 1280 680]);
tl3 = tiledlayout(fig3, 2, 2, 'TileSpacing','compact', 'Padding','compact');

% —— (a) 左上：实部叠加 ——
% 蓝色 = 扫 k → ω；绿色 = 扫 ω → k。重合 ⇔ 传播模；分歧 ⇔ 禁带
ax3a = nexttile(tl3);
hold(ax3a, 'on');
s1 = scatter(ax3a, kscan_k, kscan_fr, 10, 'filled', ...
    'MarkerFaceColor',[0.00 0.45 0.74], 'MarkerEdgeColor','none', ...
    'DisplayName','fixed-k: Re(\omega)');
s2 = scatter(ax3a, wscan_kr, wscan_f, 10, 'filled', ...
    'MarkerFaceColor',[0.00 0.62 0.45], 'MarkerEdgeColor','none', ...
    'DisplayName','fixed-\omega: Re(k)');
hold(ax3a, 'off');
xlim(ax3a, [-0.5 0.5]);
ylim(ax3a, [0 fMax]);
xlabel(ax3a, 'k\Lambda/(2\pi)');
ylabel(ax3a, '\omega\Lambda/(2\pi c)');
title(ax3a, '(a) Real parts: both methods on (k,\omega) plane');
legend(ax3a, [s1, s2], 'Location','northeast', 'FontSize',7);
grid(ax3a, 'on');
box(ax3a, 'on');

% —— (b) 右上：Im(ω) vs k（来自扫 k）——
ax3b = nexttile(tl3);
scatter(ax3b, kscan_k, kscan_fi, 10, ...
    'MarkerEdgeColor',[0.00 0.45 0.74], 'LineWidth',1.0, ...
    'DisplayName','Im(\omega)');
xlim(ax3b, [-0.5 0.5]);
xlabel(ax3b, 'k\Lambda/(2\pi)');
ylabel(ax3b, 'Im(\omega)\Lambda/(2\pi c)');
title(ax3b, '(b) Im(\omega) vs k  -- momentum gap  (\bullet = Re in Fig.1)');
grid(ax3b, 'on');
box(ax3b, 'on');
yline(ax3b, 0, 'k-', 'LineWidth',0.8);

% —— (c) 左下：Im(k) vs ω（来自扫 ω）——
ax3c = nexttile(tl3);
scatter(ax3c, wscan_ki, wscan_f, 10, ...
    'MarkerEdgeColor',[0.00 0.62 0.45], 'LineWidth',1.0, ...
    'DisplayName','Im(k)');
xlabel(ax3c, 'Im(k)\Lambda/(2\pi)');
ylabel(ax3c, '\omega\Lambda/(2\pi c)');
title(ax3c, '(c) Im(k) vs \omega  -- frequency gap  (\bullet = Re in Fig.2)');
grid(ax3c, 'on');
box(ax3c, 'on');
xline(ax3c, 0, 'k-', 'LineWidth',0.8);

% —— (d) 右下：解说 ——
ax3d = nexttile(tl3);
axis(ax3d, 'off');
hold(ax3d, 'on');
xlim(ax3d, [0 1]);
ylim(ax3d, [0 1]);

txtX = 0.05;
txtY = 0.92;
dy   = 0.045;

text(ax3d, txtX, txtY, '(d) How to read this figure:', ...
    'FontSize',10, 'FontWeight','bold', 'Units','normalized');
txtY = txtY - 1.6*dy;

text(ax3d, txtX, txtY, ...
    '\bullet Panel (a): (k, \omega) real-part overlay.', ...
    'FontSize',9, 'Units','normalized');
txtY = txtY - dy;
text(ax3d, txtX + 0.03, txtY, ...
    'Blue  \bullet = fixed-k  scan (solve \omega)', ...
    'FontSize',9, 'Color',[0.00 0.45 0.74], 'Units','normalized');
txtY = txtY - dy;
text(ax3d, txtX + 0.03, txtY, ...
    'Green \bullet = fixed-\omega scan (solve k)', ...
    'FontSize',9, 'Color',[0.00 0.62 0.45], 'Units','normalized');
txtY = txtY - 1.2*dy;

text(ax3d, txtX, txtY, ...
    '\bullet Panel (b): Im(\omega) vs k.', ...
    'FontSize',9, 'Units','normalized');
txtY = txtY - dy;
text(ax3d, txtX + 0.03, txtY, ...
    'Open circles = Im(\omega) from fixed-k scan.', ...
    'FontSize',9, 'Units','normalized');
txtY = txtY - dy;
text(ax3d, txtX + 0.03, txtY, ...
    'Im(\omega) \neq 0 \rightarrow momentum gap (growth/decay)', ...
    'FontSize',9, 'Units','normalized');
txtY = txtY - 1.2*dy;

text(ax3d, txtX, txtY, ...
    '\bullet Panel (c): Im(k) vs \omega.', ...
    'FontSize',9, 'Units','normalized');
txtY = txtY - dy;
text(ax3d, txtX + 0.03, txtY, ...
    'Open circles = Im(k) from fixed-\omega scan.', ...
    'FontSize',9, 'Units','normalized');
txtY = txtY - dy;
text(ax3d, txtX + 0.03, txtY, ...
    'Im(k) \neq 0 \rightarrow frequency gap (evanescent in space)', ...
    'FontSize',9, 'Units','normalized');
txtY = txtY - 1.2*dy;

text(ax3d, txtX, txtY, ...
    '\bullet Panel (a) mismatch:', ...
    'FontSize',9, 'Units','normalized');
txtY = txtY - dy;
text(ax3d, txtX + 0.03, txtY, ...
    'Blue-only \rightarrow \omega-gap (fixed-\omega scan misses it)', ...
    'FontSize',9, 'Units','normalized');
txtY = txtY - dy;
text(ax3d, txtX + 0.03, txtY, ...
    'Green-only \rightarrow k-gap (fixed-k scan misses it)', ...
    'FontSize',9, 'Units','normalized');

title(tl3, ...
    'Figure 3: Combined  —  complementary views of complex band structure');

% =========================================================================
% 输出
% =========================================================================

fprintf('\n=== Dual-sweep summary ===\n');
fprintf('Figure 1: Fixed-k  scan — %d points, %d modes retained\n', ...
    Nk, numel(kscan_fr));
fprintf('Figure 2: Fixed-w  scan — %d points, %d modes retained\n', ...
    Nw, numel(wscan_f));
fprintf('Figure 3: Combined overlay\n');

tags = {'fix_k_scan', 'fix_w_scan', 'combined'};
figs = [fig1, fig2, fig3];
for i = 1:3
    outputFile = fullfile(rootDir, 'output', ...
        ['example_dual_sweep_bands_' tags{i} '_' lower(quality) '.png']);
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
