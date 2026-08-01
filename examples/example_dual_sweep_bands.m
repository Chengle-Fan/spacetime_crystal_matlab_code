function example_dual_sweep_bands(quality)
%EXAMPLE_DUAL_SWEEP_BANDS 分别用"扫 k 求 ω"和"扫 ω 求 k"计算能带，并合并对比。
%
%   example_dual_sweep_bands('quick')
%   example_dual_sweep_bands('paper')
%
% 三张能带图：
%   图 1 — 固定 k → ω（2D）：仅 Re(ω) 能带，颜色 = |Im(ω)|（蓝→红 = 带隙来临）
%   图 2 — 固定 ω → k（3D）：(Re(k), Im(k), ω) 三维视图
%           传播模落在 Im(k)=0 平面内，复模向 ±Im(k) 延伸 → 频率带隙指示
%   图 3 — 叠加对比（2D）：两种方法的实部解画在同一 (k,ω) 平面上
%           重合 = 传播模，分歧 = 禁带（一种方法有解、另一种无解）

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
        Nk = 191;
        Nw = 161;
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
% 图 1：固定 k → 求 ω（常规频带结构，2D）
%   仅 Re(ω) 能带，颜色编码 |Im(ω)|：蓝=纯传播，红=近带隙
% =========================================================================

fig1 = figure('Color','w','Position',[60 400 680 520]);
ax1 = axes(fig1);
hold(ax1, 'on');

scatter(ax1, kscan_k, kscan_fr, 12, abs(kscan_fi), 'filled', ...
    'MarkerEdgeColor','none');
colormap(ax1, [linspace(0,1,256).', linspace(0.447,0,256).', linspace(0.741,0,256).']);
%  深蓝(|Im|≈0, 纯传播) → 青 → 红(|Im|≈max, 近带隙)

xlim(ax1, [-0.5 0.5]);
ylim(ax1, [0 fMax]);
xlabel(ax1, 'k\Lambda/(2\pi)');
ylabel(ax1, '\omega\Lambda/(2\pi c)');
title(ax1, 'Figure 1: Fixed-k band structure  (color = |Im(\omega)|, blue→red = gap onset)');
cb1 = colorbar(ax1);
cb1.Label.String = '|Im(\omega)|\Lambda/(2\pi c)';
grid(ax1, 'on');
box(ax1, 'on');
fprintf('Figure 1: %d Re(omega) modes (fixed-k scan)\n', numel(kscan_fr));

% =========================================================================
% 图 2：固定 ω → 求 k（复动量能带结构，3D）
%   三维空间 (Re(k), Im(k), ω)：传播模落在 Im(k)=0 平面内形成能带；
%   频率带隙处 Im(k) ≠ 0，能带向 ±Im(k) 方向延伸形成复环路
% =========================================================================

fig2 = figure('Color','w','Position',[60 400 780 600]);
ax2 = axes(fig2);
hold(ax2, 'on');

% 按 |Im(k)| 分离传播模与复模（用于区分颜色/标记）
imKmag = abs(wscan_ki);
isPropK = imKmag < 1e-3;        % 传播模：Im(k) ≈ 0
isGapK   = ~isPropK;             % 复模：Im(k) ≠ 0，能带进入复平面

% 传播模 (Im(k) ≈ 0)：绿色实心，落在 Im=0 参考面上
scatter3(ax2, wscan_kr(isPropK), zeros(sum(isPropK),1), wscan_f(isPropK), ...
    8, 'filled', 'MarkerFaceColor',[0.00 0.62 0.45], 'MarkerEdgeColor','none');

% 复模 (Im(k) ≠ 0)：橙色，Im 越大颜色越深
if any(isGapK)
    scatter3(ax2, wscan_kr(isGapK), wscan_ki(isGapK), wscan_f(isGapK), ...
        14, 'filled', 'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerEdgeColor','none');
end

% 半透明参考面 Im(k) = 0：传播模的"基底"，复模的偏离方向一目了然
[XXref, ZZref] = meshgrid(linspace(-0.5, 0.5, 20), linspace(0, fMax, 20));
YYref = zeros(size(XXref));
surf(ax2, XXref, YYref, ZZref, 'FaceAlpha',0.08, 'EdgeColor','none', ...
    'FaceColor',[0.5 0.5 0.5]);

% 阴影投影：所有点投影到 Im(k)=0 面（浅灰），辅助确认能带在传播面的位置
scatter3(ax2, wscan_kr, zeros(size(wscan_ki)), wscan_f, ...
    3, 'MarkerFaceColor',[0.85 0.85 0.85], 'MarkerEdgeColor','none');

% 竖直参考线
plot3(ax2, [0 0], [0 0], [0 fMax], 'k-', 'LineWidth',1.2);

hold(ax2, 'off');
xlabel(ax2, 'Re(k)\Lambda/(2\pi)');
ylabel(ax2, 'Im(k)\Lambda/(2\pi)');
zlabel(ax2, '\omega\Lambda/(2\pi c)');
title(ax2, 'Figure 2: Fixed-\omega complex-k band structure (3D)');

xlim(ax2, [-0.5 0.5]);
imSpanK = max([abs(min(wscan_ki)), abs(max(wscan_ki)), 0.01]) * 1.3;
ylim(ax2, [-imSpanK, imSpanK]);
zlim(ax2, [0 fMax]);

view(ax2, -38, 24);      % 视角：能同时看到传播面和 Im(k) 方向
grid(ax2, 'on');
box(ax2, 'on');

fprintf('Figure 2: %d modes (3D), %d propagating + %d complex-k (fixed-omega scan)\n', ...
    numel(wscan_f), sum(isPropK), sum(isGapK));

% =========================================================================
% 图 3：叠加对比 —— 两种方法的实部解画在同一 (k, ω) 平面上
%   蓝色 = k-scan Re(ω)，绿色 = ω-scan Re(k)
%   重合 → 传播模（两种方法都有实数解）
%   分歧 → 禁带（仅一种方法有解，另一方法在该区域无传播模）
% =========================================================================

fig3 = figure('Color','w','Position',[60 400 680 520]);
ax3 = axes(fig3);
hold(ax3, 'on');

h3a = scatter(ax3, kscan_k, kscan_fr, 8, 'filled', ...
    'MarkerFaceColor',[0.00 0.45 0.74], 'MarkerEdgeColor','none');
h3b = scatter(ax3, wscan_kr, wscan_f, 8, 'filled', ...
    'MarkerFaceColor',[0.00 0.62 0.45], 'MarkerEdgeColor','none');

hold(ax3, 'off');
xlim(ax3, [-0.5 0.5]);
ylim(ax3, [0 fMax]);
xlabel(ax3, 'k\Lambda/(2\pi)');
ylabel(ax3, '\omega\Lambda/(2\pi c)');
title(ax3, 'Figure 3: Real-part overlay — fixed-k [blue] + fixed-\omega [green]');
legend(ax3, [h3a, h3b], {'fixed-k: Re(\omega)','fixed-\omega: Re(k)'}, ...
    'Location','northeast', 'FontSize',8);
grid(ax3, 'on');
box(ax3, 'on');

% =========================================================================
% 输出
% =========================================================================

fprintf('\n=== Dual-sweep summary ===\n');
fprintf('Figure 1: Fixed-k scan — %d k-points, %d modes  [2D]\n', ...
    Nk, numel(kscan_fr));
fprintf('Figure 2: Fixed-w scan — %d w-points, %d modes (%d prop + %d complex)  [3D]\n', ...
    Nw, numel(wscan_f), sum(isPropK), sum(isGapK));
fprintf('Figure 3: Overlay  [2D, (k,w) plane]\n');

tags = {'fix_k_scan', 'fix_w_scan_3d', 'overlay'};
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
