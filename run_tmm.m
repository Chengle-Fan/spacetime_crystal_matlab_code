%RUN_TMM  用 TMM 计算方波光学时间晶体的复准频率能带。
%
% 通常只需修改下方“用户参数”区域。脚本只考虑由两个常参数时间层组成
% 的方波调制，计算两条复准频率能带，并分别绘制其实部和虚部。采用
% 归一化单位 c0=eps0=mu0=1，因此 k 与频率具有相同量纲。
%
% TMM 在每个时间层内使用解析传播矩阵，不需要设置时间步长或截断阶数。
% 准频率约定为 exp(-i*omega*t)：正虚部对应增长，负虚部对应衰减。

%% ===================== 用户参数 =====================

% 方波的两个材料状态。每个 epsilon 和 mu 都必须是有限正实数。
% 一个周期内先取高状态 epsHigh/muHigh，再取低状态 epsLow/muLow。
epsHigh = 4;
epsLow  = 1;
muHigh  = 1;
muLow   = 1;

% 时间结构：T 是完整调制周期；dutyCycle 是高状态持续时间占比。
T = 1;
dutyCycle = 0.5;       % 高状态 epsHigh/muHigh 所占比例，范围为 (0,1)

% 波数扫描使用 k/Omega 归一化，其中 Omega=2*pi/T。
kMinNormalized = -1.65;
kMaxNormalized =  1.65;
nK = 301;

%% ===================== TMM 能带计算 =====================

if dutyCycle <= 0 || dutyCycle >= 1
    error('dutyCycle must be strictly between 0 and 1.');
end

Omega = 2*pi/T;
% 两个持续时间之和严格等于 T；数组顺序同时规定时间层的作用顺序。
epsLayers = [epsHigh epsLow];
muLayers = [muHigh muLow];
durations = T*[dutyCycle 1-dutyCycle];
kNormalized = linspace(kMinNormalized,kMaxNormalized,nK);
kScan = kNormalized*Omega;
bands = tmm_bands(kScan,epsLayers,muLayers,durations);

%% ===================== 绘制能带实部和虚部 =====================

% 每个波数对应两个准频率根，因此复制横坐标以便统一展开后绘制散点。
kPlot = repmat(kNormalized,2,1);

figure('Color','w','Position',[100 100 1050 440]);
layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

axReal = nexttile(layout);
scatter(axReal,kPlot(:),real(bands.omega(:))/Omega,12,'filled');
xlabel(axReal,'k/\Omega');
ylabel(axReal,'Re(\omega)/\Omega');
title(axReal,'准频率实部');
grid(axReal,'on');
box(axReal,'on');

axImag = nexttile(layout);
scatter(axImag,kPlot(:),imag(bands.omega(:))/Omega,12,'filled');
xlabel(axImag,'k/\Omega');
ylabel(axImag,'Im(\omega)/\Omega');
title(axImag,'准频率虚部');
grid(axImag,'on');
box(axImag,'on');

title(layout,'方波光学时间晶体');
