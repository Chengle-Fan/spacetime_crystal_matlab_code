function result = reproduce_figure1_band()
%REPRODUCE_FIGURE1_BAND  复现论文图 1(b) 的方波光学时间晶体能带。
% 仅画能带与动量带隙；按用户要求不计算、也不标注 Zak 相位。

close all;

% 本文件只负责论文参数和作图；解析能带必须来自根目录 TMM 模板。
[templateSource,templatePathCleanup] = ...
    use_root_templates({'tmm_bands'});
if isempty(templatePathCleanup)
    error('根目录模板路径清理对象创建失败。');
end

%% ===================== 论文参数 =====================

% 论文最终发表版给出 epsilon1=3、epsilon2=1、T=2 fs，两个时间层
% 各占半周期。真空光速采用微米/飞秒单位，便于同时解释图 1 和图 2。
c0 = 0.299792458;          % 真空光速，单位 um/fs
epsilonLayers = [3 1];
muLayers = 1;
layerDurations = [1 1];    % 单位 fs，总周期 T=2 fs
T = sum(layerDurations);

% 论文横轴使用 k/k0，其中 k0=2*pi/(T*c0)。tmm_bands 在 c0=1 的
% 光行时坐标中工作，因此传入 kTemporal=c0*kPhysical。
k0 = 2*pi/(T*c0);
kNormalized = linspace(0,5,2001);
kPhysical = kNormalized*k0;
kTemporal = c0*kPhysical;

%% ===================== 精确时间传输矩阵能带 =====================

bands = tmm_bands(kTemporal,epsilonLayers,muLayers,layerDurations);
omegaT = bands.omega*T;

% 复准频率的非零虚部表示动量带隙。容差远低于绘图分辨率，只用于过滤
% 浮点舍入产生的约 1e-15 伪虚部。
gapTolerance = 1e-9;
growthT = max(abs(imag(omegaT)),[],1);
gapMask = growthT > gapTolerance;

%% ===================== 绘制论文图 1(b) 的无 Zak 版本 =====================

figure('Color','w','Position',[100 80 820 620]);
bandAxes = axes();
hold(bandAxes,'on');

% 先画灰色动量带隙，使蓝色通带始终位于前景。区间边缘取相邻 k 网格
% 中点，避免阴影系统性偏宽一个离散点。
gapIntervals = logical_intervals(gapMask);
for intervalIndex = 1:size(gapIntervals,1)
    firstIndex = gapIntervals(intervalIndex,1);
    lastIndex = gapIntervals(intervalIndex,2);
    leftEdge = interval_edge(kNormalized,firstIndex,-1);
    rightEdge = interval_edge(kNormalized,lastIndex,+1);
    patch(bandAxes,[leftEdge rightEdge rightEdge leftEdge], ...
        [-4 -4 4 4],[0.88 0.88 0.88], ...
        'EdgeColor','none','HandleVisibility','off');
end

% 论文蓝线只表示实准频率通带；带隙内部的复频率实部不画成伪通带。
bandColor = [0.16 0.34 0.68];
for branch = 1:2
    branchFrequency = real(omegaT(branch,:));
    branchFrequency(gapMask) = NaN;
    plot(bandAxes,kNormalized,branchFrequency,'Color',bandColor, ...
        'LineWidth',2.2);
end

xlabel(bandAxes,'k/k_0，k_0=2\pi/(Tc_0)');
ylabel(bandAxes,'\omega_F T');
title(bandAxes,'论文图 1(b)：方波光学时间晶体能带（未计算 Zak 相位）');
xlim(bandAxes,[0 5]);
ylim(bandAxes,[-4 4]);
xticks(bandAxes,0:1:5);
yticks(bandAxes,-4:2:4);
box(bandAxes,'on');

%% ===================== 返回最小复现数据 =====================

result.kNormalized = kNormalized;
result.packageVersion = '3.1.0';
result.omegaT = omegaT;
result.gapMask = gapMask;
result.gapIntervals = gapIntervals;
result.paperParameters = struct('c0',c0,'epsilonLayers',epsilonLayers, ...
    'muLayers',muLayers,'layerDurationsFs',layerDurations,'k0PerUm',k0);
result.templateSource = templateSource;
end

% -------------------------------------------------------------------------
function intervals = logical_intervals(mask)
% 把逻辑掩码转换为若干个闭区间下标 [first,last]。
edges = diff([false mask false]);
starts = find(edges == 1);
stops = find(edges == -1)-1;
intervals = [starts(:) stops(:)];
end

% -------------------------------------------------------------------------
function edge = interval_edge(grid,index,direction)
% 带隙端点位于相邻采样点中间；扫描边界处则直接使用首末坐标。
if direction < 0
    if index == 1
        edge = grid(1);
    else
        edge = 0.5*(grid(index-1)+grid(index));
    end
else
    if index == numel(grid)
        edge = grid(end);
    else
        edge = 0.5*(grid(index)+grid(index+1));
    end
end
end
