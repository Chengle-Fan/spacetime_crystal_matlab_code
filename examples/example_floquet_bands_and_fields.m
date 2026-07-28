function example_floquet_bands_and_fields(quality)
%EXAMPLE_FLOQUET_BANDS_AND_FIELDS Floquet band structure and eigenmode field patterns.
%
%   example_floquet_bands_and_fields('quick')
%   example_floquet_bands_and_fields('paper')
%
% Demonstrates the full ST-PWE workflow on a partially modulated unit cell:
%   1. Build the system matrix from analytic Fourier coefficients
%   2. Compute the Floquet band structure (fixed-k eigenvalue problem)
%   3. Select six representative eigenmodes across the Brillouin zone
%   4. Reconstruct and plot the E(x,t) field for each selected mode
%
% The color scale uses m=0 Floquet-sector participation as a transparent
% proxy for physical mode weight.

% =========================================================================
% 参数设置
% =========================================================================

if nargin < 1
    quality = 'quick';          % 默认使用 'quick'（低分辨率、快速验证）
end
rootDir = stm_init();           % 初始化工具包路径，返回输出根目录
p = stm_preset_modulated_slab(); % 获取标准参数结构体：
                                %   Lambda=1, c0=1, g=2π,
                                %   ε1=2（静态区）, εc=6（调制区背景）,
                                %   modDepth=0.6, OmegaBar=0.20,
                                %   调制区从 3Λ/4 到 Λ

% --- 精度控制 ---
switch lower(quality)
    case 'paper'
        Nspace = 20;            % 空间谐波截断：n = -20..20 → 41 个
        Nk = 181;               % k 点采样数（偶数，避开 k=0 奇异性）
    case 'quick'
        Nspace = 10;            % 空间谐波截断：n = -10..10 → 21 个
        Nk = 101;               % k 点采样数
    otherwise
        error('quality must be ''quick'' or ''paper''.');
end
Mtime = 1;                      % 时间谐波截断：m = -1, 0, 1

% --- 构造系统矩阵 ---
% epsCoeff(m,n) 返回时空介电常数的 Fourier 系数 ε_{m,n}
epsCoeff = @(m,n) stm_fourier_modulated_slab(m,n,p);
sys = stpwe_build_system(epsCoeff, [], Nspace, Mtime, p.g, p.Omega);
% sys.S = (2*Nspace+1)*(2*Mtime+1) = 总谐波数
% 矩阵尺寸 = 2*sys.S × 2*sys.S

% =========================================================================
% 扫 k 计算 Floquet 能带（固定 k，求本征频率 ω）
% =========================================================================

kBar = linspace(-0.5, 0.5, Nk);     % 归一化波矢：kΛ/(2π)
kValues = p.g*kBar;                  % 物理波矢：k = kBar * (2π/Λ)
fMax = 0.82;                         % 绘图的频率上限（归一化 fBar）
weightThreshold = 0.035;             % m=0 参与度阈值：低于此的模视为"非物理"不绘制
imagTolerance = 3e-3;                % 允许的虚部容差：衰减/增长模的筛选

% 预分配 cell 数组，每个 k 点存入选模的 k、fBar、m0Weight
allK = cell(Nk,1);
allF = cell(Nk,1);
allWeight = cell(Nk,1);

fprintf('Floquet bands: %d k points, matrix size %d x %d.\n', ...
    Nk, 2*sys.S, 2*sys.S);
for ik = 1:Nk
    % 求解固定 k 的广义本征问题 A Φ = ω B Φ
    sol = stpwe_solve_omega(sys, kValues(ik));

    % 归一化频率 fBar = ω/(g*c0) = ωΛ/(2πc)
    fBar = sol.omega/(p.g*p.c0);

    % 筛选：有限值、正频率、不超过 fMax、小虚部、有足够的 m=0 参与度
    keep = isfinite(fBar) ...
        & real(fBar) >= 0 & real(fBar) <= fMax ...
        & abs(imag(fBar)) <= imagTolerance ...
        & sol.m0Weight >= weightThreshold;

    % 收集通过的模（一个 k 点可能选到多个模）
    allK{ik} = repmat(kBar(ik), sum(keep), 1);
    allF{ik} = real(fBar(keep));
    allWeight{ik} = sol.m0Weight(keep);
end
% 将所有 k 点的结果合并为三个列向量（用于 scatter 绘图）
allK = vertcat(allK{:});
allF = vertcat(allF{:});
allWeight = vertcat(allWeight{:});

% =========================================================================
% 静态参考能带（无时间调制，m=0 仅空间调制）
% =========================================================================
% epsCoeff0(n) = ε_{0,n}：时间直流分量的 Fourier 系数
epsCoeff0 = @(n) stm_fourier_modulated_slab(0,n,p);
% 计算前 8 条静态能带（用于图 b 中的黑色实线和灰色虚线参考）
fStatic = stpwe_static_bands(kValues, epsCoeff0, Nspace, ...
    p.g, p.c0, 8);

% =========================================================================
% 选取 6 个代表性本征模，用于绘制场分布图 (c)-(h)
% =========================================================================

labels = {'c','d','e','f','g','h'};
kTargetBar = [-0.175, 0.175, 0.10, 0.29, 0.10, 0.29];
fTargetBar = [0.10, 0.10, 0.45, 0.43, 0.35, 0.32];
modes = cell(1,numel(labels));
for q = 1:numel(labels)
    % stpwe_select_mode 在 (k_target, ω_target) 附近找最近的本征模
    % 最后一个参数 [0, 0.9*p.g*p.c0] 是频率搜索范围
    modes{q} = stpwe_select_mode(sys, p.g*kTargetBar(q), ...
        p.g*p.c0*fTargetBar(q), [0, 0.9*p.g*p.c0]);
    fprintf('(%s): kBar=% .4f, fBar=% .5f%+.2ei, m0 weight=%.3f\n', ...
        labels{q}, modes{q}.k/p.g, ...
        real(modes{q}.omega/(p.g*p.c0)), ...
        imag(modes{q}.omega/(p.g*p.c0)), modes{q}.m0Weight);
end

% =========================================================================
% 绘图：3×4 拼接布局
% =========================================================================

fig = figure('Color','w','Position',[60 60 1450 850]);
tl = tiledlayout(fig, 3, 4, 'TileSpacing','compact', ...
    'Padding','compact');

% -------------------------------------------------------------------------
% 图 (a)：时空介电常数 ε(x,t) 的三维曲面
% -------------------------------------------------------------------------

axA = nexttile(tl, 1, [1 2]);           % 跨第 1、2 列
xPlot = linspace(0, 3*p.Lambda, 601);   % 3 个原胞
tauPlot = linspace(-1, 1, 161);         % 归一化时间 t/T
[XX, Tau] = meshgrid(xPlot, tauPlot);
epsPlot = stm_permittivity_modulated_slab(XX, Tau*p.T, p);
surf(axA, XX/p.Lambda, Tau, epsPlot, 'EdgeColor','none');
view(axA, -38, 28);
xlabel(axA, 'Position x/\Lambda');
ylabel(axA, 'Time t/T');
zlabel(axA, '\epsilon_r');
title(axA, '(a) Space-time permittivity');
axis(axA, 'tight');
colormap(axA, hot(256));

% -------------------------------------------------------------------------
% 图 (b)：Floquet 能带结构
% -------------------------------------------------------------------------

axB = nexttile(tl, 5, [2 2]);           % 第 5 格起，占 2×2
hold(axB, 'on');

% --- 静态参考能带：黑色实线（无时间调制的 1D 光子晶体能带）---
for ib = 1:size(fStatic,1)
    plot(axB, kBar, fStatic(ib,:), 'Color',[0.25 0.25 0.25], ...
        'LineWidth',0.8);
    % --- 灰色虚线：静态能带的 umklapp 复制（平移 ±Ω, ±2Ω）---
    for shift = [-2 -1 1 2]
        plot(axB, kBar, fStatic(ib,:) + shift*p.OmegaBar, '--', ...
            'Color',[0.65 0.65 0.65], 'LineWidth',0.45);
    end
end

% --- 彩色散点：完整 Floquet 本征模 ---
% 颜色 = m=0 扇区参与度（从品红=低到蓝=高）
scatter(axB, allK, allF, 8, allWeight, 'filled');

% --- 标记选中的 6 个模 (c)-(h) ---
for q = 1:numel(labels)
    kb = modes{q}.k/p.g;
    fb = real(modes{q}.omega/(p.g*p.c0));
    plot(axB, kb, fb, 'ko', 'MarkerFaceColor','w', 'MarkerSize',4);
    text(axB, kb + 0.012, fb + 0.012, ['(' labels{q} ')'], ...
        'FontSize',9);
end

xlim(axB, [-0.5 0.5]);
ylim(axB, [0 fMax]);
xlabel(axB, 'k\Lambda/(2\pi)');
ylabel(axB, '\omega\Lambda/(2\pi c)');
title(axB, '(b) Floquet band structure');
grid(axB, 'on');
box(axB, 'on');
% 自定义色彩映射：品红 → 蓝
cmapBand = [linspace(1,0,256).', zeros(256,1), ones(256,1)];
colormap(axB, cmapBand);
caxis(axB, [0 1]);
cb = colorbar(axB, 'northoutside');
cb.Label.String = 'm=0 sector participation';

% -------------------------------------------------------------------------
% 图 (c)-(h)：6 个本征模的时空场分布 E(x,t)
% -------------------------------------------------------------------------

xField = linspace(0, 3*p.Lambda, 451);      % 3 个原胞
tauField = linspace(0, 3, 361);              % 3 个调制周期
fieldTiles = [3 4 7 8 11 12];               % 在 3×4 布局中的 tile 位置
for q = 1:numel(labels)
    ax = nexttile(tl, fieldTiles(q));

    % 从本征向量重构时空电场 E(x,t)
    field = stpwe_reconstruct_field(modes{q}, xField, ...
        tauField*p.T, false);

    imagesc(ax, xField/p.Lambda, tauField, field);
    set(ax, 'YDir','normal');
    caxis(ax, [-1 1]);                       % 色标固定 [-1,1] 便于对比
    colormap(ax, stm_redblue(256));          % 红-白-蓝发散色图

    hold(ax, 'on');
    % 标注每个原胞中调制区起始位置（xModStart 处竖线）
    for cellNo = 0:2
        xline(ax, cellNo + p.xModStart/p.Lambda, 'k-', ...
            'LineWidth',0.45);
    end

    title(ax, ['(' labels{q} ')']);
    xlabel(ax, 'x/\Lambda');
    % 左列标 y 轴，右列不标（节省空间）
    if mod(q,2) == 1
        ylabel(ax, 't/T');
    else
        set(ax, 'YTickLabel',[]);
    end
end

% =========================================================================
% 总标题与输出
% =========================================================================

title(tl, 'Floquet bands and eigenmode fields (modulated slab)');
outputFile = fullfile(rootDir, 'output', ...
    ['example_floquet_bands_and_fields_' lower(quality) '.png']);
save_example_figure(fig, outputFile);
fprintf('Saved %s\n', outputFile);
end

% =========================================================================
% 辅助函数：保存图片（支持 exportgraphics / print 两种后端）
% =========================================================================

function save_example_figure(fig, outputFile)
try
    exportgraphics(fig, outputFile, 'Resolution',220);
catch
    print(fig, outputFile, '-dpng', '-r220');
end
end
