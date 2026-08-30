%RUN_FDTD_FFT  从有限样品复场历史重建局部 FFT 能带。
% 第二模块：读取第一模块保存的复场，只绘制局部 FFT 响应谱。

clear; clc; close all;

%% ===================== 用户参数 =====================

dataFile = 'fdtd_field_data.mat';

% 从样品两端各去掉多少个几何长度单元，避免端面反射直接进入空间 FFT。
roiMarginCells = 2;

% FFT 时间窗的起点以调制周期计。fftPeriodCount=[] 时，自动使用该起点
% 之后能够容纳的最多完整周期；也可以设置为指定的正整数周期数。
fftStartPeriod = 0;
fftPeriodCount = [];

% 时间零填充只细化显示网格，不提高真实频率分辨率；必须是 >=1 的整数。
zeroPaddingFactor = 4;
displayDynamicRangeDb = 60;

% 只显示初始中心波数附近的局部谱；设为空 [] 可显示全部有符号 k。
plotKHalfWidthNormalized = 0.25;

%% ===================== 读取并选择样品内部时空窗口 =====================

if ~isfile(dataFile)
    error('找不到 %s；请先独立运行 run_fdtd_field.m。',dataFile);
end
loadedData = load(dataFile,'fdtdData');
if ~isfield(loadedData,'fdtdData')
    error('%s 中缺少 fdtdData 结构体。',dataFile);
end
fdtdData = loadedData.fdtdData;
requiredFields = {'x','t','E','temporalPeriod','Omega','sampleEdges', ...
    'cellSize','targetRegion','kCenterNormalized'};
missingFields = requiredFields(~isfield(fdtdData,requiredFields));
if ~isempty(missingFields)
    error('fdtdData 缺少字段：%s。',strjoin(missingFields,', '));
end

if ~isscalar(roiMarginCells) || ~isreal(roiMarginCells) || ...
        ~isfinite(roiMarginCells) || roiMarginCells < 0
    error('roiMarginCells 必须是非负有限实标量。');
end
if ~isempty(plotKHalfWidthNormalized) && ...
        (~isnumeric(plotKHalfWidthNormalized) || ...
        ~isscalar(plotKHalfWidthNormalized) || ...
        ~isreal(plotKHalfWidthNormalized) || ...
        ~isfinite(plotKHalfWidthNormalized) || plotKHalfWidthNormalized <= 0)
    error('plotKHalfWidthNormalized 必须为空或有限正实标量。');
end
xMargin = roiMarginCells*fdtdData.cellSize;
xRoiMin = fdtdData.sampleEdges(1)+xMargin;
xRoiMax = fdtdData.sampleEdges(2)-xMargin;
spaceIds = find(fdtdData.x >= xRoiMin & fdtdData.x < xRoiMax);
if numel(spaceIds) < 8
    error('样品内部空间 ROI 太短；请减小 roiMarginCells 或增大样品长度。');
end

timeIds = select_time_indices(fdtdData.t,fdtdData.temporalPeriod, ...
    fftStartPeriod,fftPeriodCount);

%% ===================== FFT 能带重建 =====================

% 必须变换完整复数 E，而不是显示用的 abs(E)。所得结果包含有限样品边界、
% 初值带宽和有限观测时间的影响，是源加权响应谱，不是无限体严格本征谱。
% 一次 Gaussian 初值只能覆盖其 k 支持附近；k-gap 中也只能显示实频响应，
% 不能把谱线展宽直接解释为 Im(omega)。本入口不叠加 PWE/TMM 理论曲线。
fftCfg = struct();
fftCfg.temporalPeriod = fdtdData.temporalPeriod;
fftCfg.spatialROI = [spaceIds(1) spaceIds(end)];
fftCfg.timeROI = [timeIds(1) timeIds(end)];
fftCfg.zeroPaddingFactor = zeroPaddingFactor;
fftCfg.dynamicRangeDb = displayDynamicRangeDb;
fftBands = fdtd_fft_bands(fdtdData.E,fdtdData.x,fdtdData.t,fftCfg);

%% ===================== 只绘制 FFT 重建能带 =====================

figure('Color','w','Position',[100 80 900 650]);
imagesc(fftBands.kNormalized,fftBands.omegaNormalized,fftBands.spectrumDb);
axis xy;
colormap(parula(256));
clim([-displayDynamicRangeDb 0]);
fftColorbar = colorbar;
fftColorbar.Label.String = 'FFT 相对功率 (dB)';
xlabel('k/\Omega');
ylabel('Re(\omega)/\Omega');
ylim([-0.5 0.5]);
if ~isempty(plotKHalfWidthNormalized)
    xlim(fdtdData.kCenterNormalized+[-1 1]*plotKHalfWidthNormalized);
end
title(sprintf('有限样品 FFT 局部能带：%s 激发',fdtdData.targetRegion));
box on;

% -------------------------------------------------------------------------
function timeIds = select_time_indices(t,T,startPeriod,periodCount)
% 从任意总记录中截取无重复终点、恰含整数个调制周期的连续时间窗口。
t = t(:);
if numel(t) < 2
    error('场记录至少需要两个时间点。');
end
if ~isnumeric(startPeriod) || ~isscalar(startPeriod) || ...
        ~isreal(startPeriod) || ~isfinite(startPeriod) || startPeriod < 0
    error('fftStartPeriod 必须是非负有限实标量。');
end
dtRecord = t(2)-t(1);
samplesPerPeriodFloat = T/dtRecord;
if abs(samplesPerPeriodFloat-round(samplesPerPeriodFloat)) > 1e-9
    error('场记录采样间隔与调制周期不可公度，无法做严格 Floquet 折叠。');
end
samplesPerPeriod = round(samplesPerPeriodFloat);
startTime = startPeriod*T;
startId = round((startTime-t(1))/dtRecord)+1;
if startId < 1 || startId > numel(t) || ...
        abs(t(startId)-startTime) > 1e-9*max(1,abs(startTime))
    error('fftStartPeriod 必须落在已记录的时间网格上。');
end
availablePeriods = floor((numel(t)-startId+1)/samplesPerPeriod);
if isempty(periodCount)
    selectedPeriodCount = availablePeriods;
else
    if ~isnumeric(periodCount) || ~isscalar(periodCount) || ...
            ~isreal(periodCount) || ~isfinite(periodCount) || ...
            periodCount < 1 || periodCount ~= round(periodCount)
        error('fftPeriodCount 必须为空或正整数。');
    end
    selectedPeriodCount = periodCount;
end
if selectedPeriodCount < 1 || selectedPeriodCount > availablePeriods
    error('所选时间窗没有足够的完整调制周期。');
end
timeIds = startId:startId+selectedPeriodCount*samplesPerPeriod-1;
end
