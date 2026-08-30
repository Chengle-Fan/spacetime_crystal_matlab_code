function result = reproduce_figure2_fdtd()
%REPRODUCE_FIGURE2_FDTD  复现论文图 2 的体带与 k-gap 两种脉冲场分布。
% 绘制论文观测量电位移 D 的自然对数幅度；不计算任何 Zak 相位。

close all;

% 本文件只改写论文样品、脉冲和绘图参数；D/B-Yee 与解析增长检查分别
% 强制调用根目录 fdtd1d 和 tmm_bands 模板，不保留求解器副本。
[templateSource,templatePathCleanup] = ...
    use_root_templates({'fdtd1d','tmm_bands'});
if isempty(templatePathCleanup)
    error('根目录模板路径清理对象创建失败。');
end

%% ===================== 论文参数与明确的数值假设 =====================

paper = struct();
paper.c0 = 0.299792458;          % 真空光速，单位 um/fs
paper.epsHigh = 3;
paper.epsLow = 1;
paper.muRelative = 1;
paper.temporalPeriod = 2;       % fs
paper.timeCrystalStart = 220;   % fs
paper.timeCrystalEnd = 340;     % fs，共 60 个时间周期
paper.finalTime = 500;          % fs；正文只给出图示范围，这里按图 2 取值
paper.pulseFwhm = 45;           % fs；采用最终发表版，而不是 arXiv v1 的 5 fs
paper.bandWavelength = 1.4;     % um，位于通带
paper.gapWavelength = 0.93;     % um，位于第一动量带隙

% 正文和所附补充材料没有报告 FDTD 网格、边界层或记录间隔。下面参数是
% 本次初步代码复现的显式假设：最短波长约 23 个网格点，每周期 40 步。
% 空间边界放在 500 fs 内任何主脉冲都到达不了的位置，并关闭海绵层。
% 这是因为普通海绵会打破空间平移对称性、散射出额外 k 分量；这些极小
% 分量在 k-gap 中也会被时间晶体指数放大，最终污染本来位于通带的算例。
% 科研定量使用前仍应继续做 dz/dt/边界距离收敛。
grid = struct();
grid.zMin = -120;               % um；使初始高斯在左端降到机器舍入以下
grid.zMax = 190;                % um；500 fs 内主脉冲不会到达右端
grid.plotZMin = 0;              % um，与论文横轴一致
grid.plotZMax = 130;            % um
grid.dz = 0.04;                 % um
grid.dt = 0.05;                 % fs
grid.recordEvery = 4;           % 每 0.2 fs 记录一次
grid.spongeWidth = 0;           % um；远边界方案不使用会散射波数的海绵层
grid.spongeRate = 0;            % 1/fs

%% ===================== 分别计算体带和 k-gap 脉冲 =====================

bandCase = simulate_paper_pulse(paper.bandWavelength,paper,grid);
gapCase = simulate_paper_pulse(paper.gapWavelength,paper,grid);

% 用同一组连续 Maxwell 方程的精确 TMM 给出中心波数理论增长率。它不是
% 用来替代 FDTD，而是用来判断长时间指数增长究竟是物理效应还是 CFL
% 数值发散。若中心波数位于通带，下面的每周期对数增长应为零。
bandTheory = floquet_growth_diagnostic(paper.bandWavelength,paper);
gapTheory = floquet_growth_diagnostic(paper.gapWavelength,paper);

%% ===================== 只绘制论文图 2 的两个场分布面板 =====================

figure('Color','w','Position',[80 80 1120 560]);
layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

bandAxes = nexttile(layout);
imagesc(bandAxes,bandCase.z,bandCase.t,bandCase.logAmplitude);
axis(bandAxes,'xy');
bandAxes.CLim = [-1.5 0.35];
colormap(bandAxes,jet(256));
bandColorbar = colorbar(bandAxes);
bandColorbar.Label.String = 'ln(|D|/D_0)';
xlabel(bandAxes,'z (\mum)');
ylabel(bandAxes,'t (fs)');
title(bandAxes,'(a) Band propagation，\lambda_0=1.4 \mum');
xlim(bandAxes,[grid.plotZMin grid.plotZMax]);
ylim(bandAxes,[0 paper.finalTime]);
draw_time_crystal_bracket(bandAxes,paper,grid.plotZMax-5);
box(bandAxes,'on');

gapAxes = nexttile(layout);
imagesc(gapAxes,gapCase.z,gapCase.t,gapCase.logAmplitude);
axis(gapAxes,'xy');

% 论文正文称 k-gap 脉冲约放大 2e4，对应 ln(2e4)≈9.9；固定到 10.5
% 可直接比较论文色条。若本实现得到更大增长，超出部分会饱和为红色，
% 实际最大值仍保存在 result.gap.maximumLogAmplitude 中。
gapAxes.CLim = [-1.5 10.5];
colormap(gapAxes,jet(256));
gapColorbar = colorbar(gapAxes);
gapColorbar.Label.String = 'ln(|D|/D_0)';
xlabel(gapAxes,'z (\mum)');
ylabel(gapAxes,'t (fs)');
title(gapAxes,'(b) Gap propagation，\lambda_0=0.93 \mum');
xlim(gapAxes,[grid.plotZMin grid.plotZMax]);
ylim(gapAxes,[0 paper.finalTime]);
draw_time_crystal_bracket(gapAxes,paper,grid.plotZMax-5);
box(gapAxes,'on');

title(layout,'论文图 2：有限时长光学时间晶体中的电位移场');

%% ===================== 返回复现数据与诊断量 =====================

result.paperParameters = paper;
result.gridParameters = grid;
result.band = bandCase;
result.gap = gapCase;
result.floquetTheory = struct('band',bandTheory,'gap',gapTheory);
result.paperReportedGapAmplification = 20000;
result.templateSource = templateSource;

fprintf('图 2 体带算例：max ln(|D|/D0) = %.6g。\n', ...
    bandCase.maximumLogAmplitude);
fprintf('图 2 k-gap 算例：max ln(|D|/D0) = %.6g，幅度比约 %.6g。\n', ...
    gapCase.maximumLogAmplitude,exp(gapCase.maximumLogAmplitude));
fprintf(['同参数精确 TMM 给出的 60 周期中心波数对数增长为 %.6g；' ...
    '论文正文 20000 倍对应 %.6g。\n'], ...
    gapTheory.logGainOverCrystal,log(result.paperReportedGapAmplification));
end

% -------------------------------------------------------------------------
function pulseCase = simulate_paper_pulse(centerWavelength,paper,grid)
% 在光行时坐标 x=z/c0 中调用当前 D/B-Yee 内核，此时真空光速归一化为 1。
z = grid.zMin:grid.dz:grid.zMax;
x = z/paper.c0;
dx = x(2)-x(1);
dt = grid.dt;
nSteps = round(paper.finalTime/dt);

if abs(nSteps*dt-paper.finalTime) > 1e-12
    error('paper.finalTime 必须是 grid.dt 的整数倍。');
end
if mod(nSteps,grid.recordEvery) ~= 0
    error('总 FDTD 步数必须能被 grid.recordEvery 整除。');
end

% 论文给出时间 FWHM。这里按强度 FWHM 解释，所以复振幅包络写成
% exp[-2*ln(2)*(z/L)^2]，其中 L=c0*pulseFwhm。
spatialFwhm = paper.c0*paper.pulseFwhm;
initialCenter = 0;
electricEnvelope = exp(-2*log(2)*((z-initialCenter)/spatialFwhm).^2);

physicalWavenumber = 2*pi/centerWavelength;
temporalWavenumber = paper.c0*physicalWavenumber;
E0 = electricEnvelope.*exp(1i*temporalWavenumber*x);

% H 位于 t=-dt/2 和空间半格。除载波的半时间步相位外，也把高斯中心
% 回退 c0*dt/2，使窄带右行初值在包络和载波两方面保持 Yee 错位一致。
xH = x(1:end-1)+dx/2;
zH = paper.c0*xH;
centerAtMinusHalfStep = initialCenter-paper.c0*dt/2;
magneticEnvelope = exp(-2*log(2)* ...
    ((zH-centerAtMinusHalfStep)/spatialFwhm).^2);
yeeArgument = (dt/dx)*sin(temporalWavenumber*dx/2);
if abs(yeeArgument) >= 1
    error('当前载波和网格不满足 Yee 离散色散稳定条件。');
end
omegaYee = (2/dt)*asin(yeeArgument);
Hhalf0 = magneticEnvelope.*exp(1i*temporalWavenumber*xH) .* ...
    exp(1i*omegaYee*dt/2);

% PTC 在 220–340 fs 内重复“高 epsilon 1 fs、低 epsilon 1 fs”。用整数
% 时间步生成 220,221,...,339 fs 的真实突变节点，避免浮点累积漂移。
halfPeriod = paper.temporalPeriod/2;
interfaceStepSpacing = round(halfPeriod/dt);
firstInterfaceStep = round(paper.timeCrystalStart/dt);
lastInterfaceStep = round((paper.timeCrystalEnd-halfPeriod)/dt);
interfaceSteps = firstInterfaceStep:interfaceStepSpacing:lastInterfaceStep;
temporalInterfaces = interfaceSteps*dt;

epsilonFunction = @(position,time) temporal_epsilon(time,paper).* ...
    ones(size(position));
muFunction = @(position,time) paper.muRelative*ones(size(position));

fdtdCfg = struct();
fdtdCfg.x = x;
fdtdCfg.dt = dt;
fdtdCfg.nSteps = nSteps;
fdtdCfg.epsFun = epsilonFunction;
fdtdCfg.muFun = muFunction;
fdtdCfg.E0 = E0;
fdtdCfg.Hhalf0 = Hhalf0;
fdtdCfg.recordEvery = grid.recordEvery;
fdtdCfg.recordField = 'D';
fdtdCfg.recordPrecision = 'single';
plotMask = z >= grid.plotZMin & z <= grid.plotZMax;
fdtdCfg.recordSpatialIndices = find(plotMask);
fdtdCfg.boundaryType = 'fixed-distant';
fdtdCfg.spongeCells = round(grid.spongeWidth/grid.dz);
fdtdCfg.spongeStrength = grid.spongeRate;
fdtdCfg.temporalInterfaces = temporalInterfaces;
fdtdCfg.certifiedMinimumEpsilon = min(paper.epsLow,paper.epsHigh);
fdtdCfg.certifiedMinimumMu = paper.muRelative;
field = fdtd1d(fdtdCfg);

% 论文画的是 D，不是 E。修复后的模板可直接记录内部连续推进变量 D，
% 并且只保存论文显示的空间 ROI，无需先保存全域 E 再做本构恢复。
amplitude = abs(field.D);
logAmplitude = single(log(max(amplitude,exp(-20))));

pulseCase = struct();
pulseCase.z = paper.c0*field.x;
pulseCase.t = field.t;
pulseCase.logAmplitude = logAmplitude;
pulseCase.maximumLogAmplitude = double(max(logAmplitude,[],'all'));
pulseCase.maximumAmplification = exp(pulseCase.maximumLogAmplitude);
pulseCase.centerWavelength = centerWavelength;
pulseCase.kOverK0 = paper.c0*paper.temporalPeriod/centerWavelength;
pulseCase.courant = field.courant;
pulseCase.cflAuditSource = field.cflAuditSource;
pulseCase.interfaceCount = numel(temporalInterfaces);
end

% -------------------------------------------------------------------------
function diagnostic = floquet_growth_diagnostic(centerWavelength,paper)
% 把真空波长转换为光行时波数，并计算一个周期 Floquet 乘子的增长指数。
kTemporal = 2*pi*paper.c0/centerWavelength;
theory = tmm_bands(kTemporal,[paper.epsHigh paper.epsLow], ...
    paper.muRelative,[paper.temporalPeriod/2 paper.temporalPeriod/2]);
omegaT = theory.omega*paper.temporalPeriod;
growthPerPeriod = max(imag(omegaT));
periodCount = round((paper.timeCrystalEnd-paper.timeCrystalStart)/ ...
    paper.temporalPeriod);

diagnostic = struct();
diagnostic.omegaT = omegaT;
diagnostic.logGainPerPeriod = growthPerPeriod;
diagnostic.periodCount = periodCount;
diagnostic.logGainOverCrystal = periodCount*growthPerPeriod;
diagnostic.amplificationOverCrystal = exp(diagnostic.logGainOverCrystal);
end

% -------------------------------------------------------------------------
function epsilon = temporal_epsilon(time,paper)
% PTC 之外保持 epsilon=1；PTC 内每周期先取高介电层，再取低介电层。
epsilon = paper.epsLow*ones(size(time));
inside = time >= paper.timeCrystalStart & time < paper.timeCrystalEnd;
phase = mod(time(inside)-paper.timeCrystalStart,paper.temporalPeriod);
epsilon(inside) = paper.epsLow+(paper.epsHigh-paper.epsLow).* ...
    double(phase < paper.temporalPeriod/2);
end

% -------------------------------------------------------------------------
function draw_time_crystal_bracket(targetAxes,paper,zPosition)
% 用白色括号标出论文中的有限时间晶体开启区间。
hold(targetAxes,'on');
capWidth = 3;
plot(targetAxes,[zPosition zPosition], ...
    [paper.timeCrystalStart paper.timeCrystalEnd],'w-','LineWidth',2);
plot(targetAxes,zPosition+[-capWidth 0], ...
    [paper.timeCrystalStart paper.timeCrystalStart],'w-','LineWidth',2);
plot(targetAxes,zPosition+[-capWidth 0], ...
    [paper.timeCrystalEnd paper.timeCrystalEnd],'w-','LineWidth',2);
text(targetAxes,zPosition-1,mean([paper.timeCrystalStart paper.timeCrystalEnd]), ...
    '\epsilon(t)','Color','w','FontWeight','bold', ...
    'HorizontalAlignment','right','VerticalAlignment','middle');
end
