%% test_temporal_interface_animation.m
% 时间界面（temporal interface / time boundary）波包演化动画。
%
% 在 t = tSwitch 时刻，介电常数发生突变 (nBefore -> nAfter)，
% 产生时间折射 (temporal refraction) 和时间反射 (temporal reflection)：
%   - 正向传播的透射波 (forward, faster/slower)
%   - 反向传播的反射波 (backward, 相位共轭)
%
% 生成 line 和 waterfall 两种模式的动画 GIF。

fprintf('=== 时间界面波包演化动画 ===\n');

%==========================================================================
% 物理参数（归一化单位: c0 = 1, lambda0 = 1）
%==========================================================================

lambda0 = 1;
k0 = 2*pi / lambda0;
c0 = 1;

% 时间界面两侧的折射率
nBefore = 1.5;
nAfter  = 2.5;
tSwitch = 15;          % 界面切换时刻
tEnd    = 30;          % 仿真总时长

% 空间网格
Lf = 32;               % 每波长采样点数
Lx = 50;               % 模拟域长度（波长数）
dx = lambda0 / Lf;
x = (0:dx:Lx*lambda0).';   % 列向量

% 时间步长 (CFL = 0.80, v_max = c0/nBefore)
dt = 0.80 * dx / c0;
nSteps = ceil(tEnd / dt);

% 高斯波包初始条件
x0 = 12 * lambda0;     % 波包初始中心位置
sigma = 3.0 * lambda0; % 波包半宽
profile = @(xq) exp(-((xq - x0)/sigma).^2) .* exp(1i*k0*(xq - x0));

E0 = profile(x).';

% H 初始条件 (t = -dt/2, H 网格) — 正向波在折射率 nBefore 介质中
xH = x(1:end-1) + dx/2;
vBefore = c0 / nBefore;
Hhalf0 = nBefore * profile(xH.' + vBefore*dt/2).';

%==========================================================================
% 时变介电常数: t < tSwitch 为 nBefore^2, t >= tSwitch 为 nAfter^2
%==========================================================================

epsFun = @(xq, tq) temporal_eps(xq, tq, tSwitch, nBefore^2, nAfter^2);

%==========================================================================
% 运行 FDTD 仿真
%==========================================================================

fprintf('\n仿真参数:\n');
fprintf('  nBefore = %.1f, nAfter = %.1f\n', nBefore, nAfter);
fprintf('  tSwitch = %.1f, tEnd = %.1f\n', tSwitch, tEnd);
fprintf('  Lx = %.0f 波长, Nx = %d 格点\n', Lx, length(x));
fprintf('  dt = %.4f, nSteps = %d\n', dt, nSteps);

cfg.x = x;
cfg.dt = dt;
cfg.nSteps = nSteps;
cfg.epsFun = epsFun;
cfg.muFun = @(xq, tq) ones(size(xq));
cfg.E0 = E0;
cfg.Hhalf0 = Hhalf0;
cfg.boundary = 'sponge';
cfg.spongeCells = 120;
cfg.spongeStrength = 0.08;
cfg.recordEvery = 2;
cfg.storeFields = true;
cfg.precision = 'double';

out = fdtd1d_db(cfg);

%==========================================================================
% 诊断: 时间界面前后的波包分解
%==========================================================================

[~, beforeId] = min(abs(out.t - (tSwitch - 2*out.dt)));
[~, afterId]  = min(abs(out.t - (tSwitch + 2*out.dt)));
Ebefore = out.E(beforeId, :);
Eafter  = out.E(afterId, :);
Hafter  = out.H(afterId, :);

% 方向分解: E_forward = (E + H/n_after)/2, E_backward = (E - H/n_after)/2
Eplus  = 0.5 * (Eafter + Hafter/nAfter);
Eminus = 0.5 * (Eafter - Hafter/nAfter);

tauNumeric = norm(Eplus) / norm(Ebefore);
rhoNumeric = norm(Eminus) / norm(Ebefore);

fprintf('\n时间界面诊断:\n');
fprintf('  |tau| (透射系数): %.5f (FDTD 数值)\n', tauNumeric);
fprintf('  |rho| (反射系数): %.5f (FDTD 数值)\n', rhoNumeric);

% 解析 Morgenthaler 系数 (用于对照)
n1 = nBefore; n2 = nAfter;
tauExact = (n1 + n2) / (2*n2);
rhoExact = (n1 - n2) / (2*n2);
fprintf('  |tau| (解析):      %.5f\n', abs(tauExact));
fprintf('  |rho| (解析):      %.5f\n', abs(rhoExact));

%==========================================================================
% 生成动画
%==========================================================================

fprintf('\n=== 生成动画 ===\n');

% --- 1. Line 模式动画: 显示 E 场实时演化 ---
fprintf('  生成 line 模式动画 (含 tSwitch 标记)...\n');
animate_fdtd1d(out, struct(...
    'saveAs',    'test_ti_line.gif', ...
    'viewMode',  'line', ...
    'fields',    'E', ...
    'fps',       15, ...
    'skipFrames', 3, ...
    'showTime',  true, ...
    'speed',     10, ...
    'title',     sprintf('Temporal Interface: n=%.1f \\rightarrow %.1f at t=%.0f', ...
                         nBefore, nAfter, tSwitch)));

% --- 2. Waterfall 模式动画: 时空演化瀑布图 ---
fprintf('  生成 waterfall 模式动画...\n');
animate_fdtd1d(out, struct(...
    'saveAs',    'test_ti_waterfall.gif', ...
    'viewMode',  'waterfall', ...
    'fields',    'E', ...
    'fps',       15, ...
    'skipFrames', 3, ...
    'speed',     10, ...
    'title',     sprintf('Temporal Interface Spacetime: n=%.1f \\rightarrow %.1f', ...
                         nBefore, nAfter)));

%==========================================================================
% 诊断图: 时间界面前后的场分布 + 方向分解
%==========================================================================

fprintf('\n=== 生成诊断快照 ===\n');

figDiag = figure('Color', 'w', 'Position', [50 50 1200 800], 'Visible', 'off');

% --- 子图 1: |E(x,t)| 时空瀑布图 (静态) ---
ax1 = subplot(2, 2, [1 2]);
imagesc(ax1, out.x, out.t, abs(out.E));
set(ax1, 'YDir', 'normal');
hold(ax1, 'on');
yline(ax1, tSwitch, 'w--', 'LineWidth', 1.5);
xlabel(ax1, 'x / \lambda_0');
ylabel(ax1, 't');
title(ax1, '|E(x,t)| — Temporal Refraction & Reflection');
colormap(ax1, parula(256));
cb1 = colorbar(ax1);
ylabel(cb1, '|E|');

% --- 子图 2: 界面前后场分布 + 方向分解 ---
ax2 = subplot(2, 2, 3);
plot(ax2, out.x, real(Ebefore), 'k-', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Before (t=%.1f)', out.t(beforeId)));
hold(ax2, 'on');
plot(ax2, out.x, real(Eplus),  'b-', 'LineWidth', 1.2, ...
    'DisplayName', 'Forward (E_+)');
plot(ax2, out.x, real(Eminus), 'r-', 'LineWidth', 1.2, ...
    'DisplayName', 'Backward (E_-)');
xlabel(ax2, 'x / \lambda_0');
ylabel(ax2, 'Re(E_y)');
title(ax2, sprintf('Directional Decomposition (t = %.1f)', out.t(afterId)));
legend(ax2, 'Location', 'best');
grid(ax2, 'on');
box(ax2, 'on');

% --- 子图 3: 能量随时间变化 ---
ax3 = subplot(2, 2, 4);
plot(ax3, out.t, out.energy / out.energy(1), 'b-', 'LineWidth', 1.4);
hold(ax3, 'on');
xline(ax3, tSwitch, 'r--', 'LineWidth', 1.3);
xlabel(ax3, 't');
ylabel(ax3, 'Energy / E(0)');
title(ax3, 'Instantaneous Energy (parametric exchange)');
grid(ax3, 'on');
box(ax3, 'on');

sgtitle(sprintf('Temporal Interface: n = %.1f \\rightarrow %.1f  (t_{switch} = %.0f)', ...
    nBefore, nAfter, tSwitch));

saveas(figDiag, 'test_ti_diagnostic.png');
close(figDiag);
fprintf('  快照已保存至: test_ti_diagnostic.png\n');

%==========================================================================
% 汇总
%==========================================================================

fprintf('\n=== 全部完成 ===\n');
fprintf('生成文件:\n');
fprintf('  test_ti_line.gif        — line 模式动画 (波包遇时间界面分裂)\n');
fprintf('  test_ti_waterfall.gif   — waterfall 模式时空演化动画\n');
fprintf('  test_ti_diagnostic.png  — 诊断快照 (分解 + 能量)\n');

%==========================================================================
% 辅助函数: 时间界面介电常数
%==========================================================================

function epsr = temporal_eps(x, t, tSwitch, epsBefore, epsAfter)
if t < tSwitch
    epsr = epsBefore * ones(size(x));
else
    epsr = epsAfter * ones(size(x));
end
end
