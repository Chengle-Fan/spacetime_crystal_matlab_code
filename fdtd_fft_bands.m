function result = fdtd_fft_bands(probeSignals,time,kScan,cfg)
%FDTD_FFT_BANDS  由逐 k 有限空间 FDTD 的固定探针 E(t) 重建 Floquet 能带。
% probeSignals 的标准尺寸为 Nt-by-Nprobe-by-Nk。MATLAB 会省略末尾的
% 单例维：nK=1 时二维 Nt-by-Nprobe 表示多探针单列；nK>1 时二维
% Nt-by-Nk 才表示单探针多列。
% 函数只做时间窗、FFT、探针功率合并及按 Omega 的 Floquet 折叠，不求解
% Maxwell 方程，也不把高斯激发中心 k 误认为每个场样本的精确本征波数。

if nargin ~= 4 || ~isstruct(cfg) || ~isscalar(cfg)
    error('必须以 probeSignals、time、kScan 和标量结构体 cfg 调用。');
end
required = {'temporalPeriod','analysisTimeRange'};
missing = required(~isfield(cfg,required));
if ~isempty(missing)
    error('cfg 缺少字段：%s。',strjoin(missing,', '));
end

%% ===================== 数据尺寸与采样网格检查 =====================

if ~isnumeric(time) || ~isvector(time) || numel(time) < 4 || ...
        ~isreal(time) || any(~isfinite(time))
    error('time 必须是至少含 4 点的有限实向量。');
end
time = time(:);
dtVector = diff(time);
if any(dtVector <= 0)
    error('time 必须严格递增。');
end
dtRecord = mean(dtVector);
timeRoundoff = 128*eps(max(abs(time)));
uniformTimeTolerance = max(1e-9*dtRecord,timeRoundoff);
if max(abs(dtVector-dtRecord)) > uniformTimeTolerance
    error('time 必须是均匀采样网格。');
end

if ~isnumeric(kScan) || ~isvector(kScan) || isempty(kScan) || ...
        ~isreal(kScan) || any(~isfinite(kScan))
    error('kScan 必须是非空有限实向量。');
end
kScan = kScan(:).';
nK = numel(kScan);

if ~isnumeric(probeSignals) || isempty(probeSignals) || ...
        any(~isfinite(probeSignals(:)))
    error('probeSignals 必须是非空有限数值数组。');
end
if ismatrix(probeSignals)
    if size(probeSignals,1) ~= numel(time)
        error('probeSignals 的第一维必须等于 numel(time)。');
    end
    % MATLAB 会省略尺寸为 1 的末维，所以 Nt-by-Nprobe-by-1
    % 会表现为普通二维矩阵。对单 k 扫描必须把第二维解释为
    % 探针，否则多探针单列数据会被错认成多个 k。
    if nK == 1
        probeSignals = reshape(probeSignals,numel(time),size(probeSignals,2),1);
    elseif size(probeSignals,2) == nK
        probeSignals = reshape(probeSignals,numel(time),1,nK);
    else
        error('二维 probeSignals 在多 k 时必须具有 Nt-by-Nk 尺寸。');
    end
elseif ndims(probeSignals) == 3
    if size(probeSignals,1) ~= numel(time) || size(probeSignals,3) ~= nK
        error('三维 probeSignals 必须具有 Nt-by-Nprobe-by-Nk 尺寸。');
    end
else
    error('probeSignals 必须是二维 Nt-by-Nk 或三维 Nt-by-Nprobe-by-Nk 数组。');
end
nProbe = size(probeSignals,2);

T = read_positive_scalar(cfg.temporalPeriod,'cfg.temporalPeriod');
Omega = 2*pi/T;
samplesPerPeriodFloat = T/dtRecord;
samplesPerPeriod = round(samplesPerPeriodFloat);
if samplesPerPeriod < 2 || ...
        abs(samplesPerPeriodFloat-samplesPerPeriod) > ...
        1e-10*max(1,samplesPerPeriod)
    error('探针记录间隔必须与 temporalPeriod 可公度，且每周期至少记录两点。');
end

analysisTimeRange = cfg.analysisTimeRange;
if ~isnumeric(analysisTimeRange) || ~isvector(analysisTimeRange) || ...
        numel(analysisTimeRange) ~= 2 || ~isreal(analysisTimeRange) || ...
        any(~isfinite(analysisTimeRange)) || ...
        analysisTimeRange(1) < time(1)-uniformTimeTolerance || ...
        analysisTimeRange(2) <= analysisTimeRange(1) || ...
        analysisTimeRange(2) > time(end)+uniformTimeTolerance
    error('cfg.analysisTimeRange 必须是记录范围内严格递增的两个有限时刻。');
end
analysisNodeFloat = (analysisTimeRange-time(1))/dtRecord;
analysisNodes = round(analysisNodeFloat);
nodeAlignmentTolerance = max(1e-8, ...
    128*eps(max(abs(analysisNodeFloat))));
if any(abs(analysisNodeFloat-analysisNodes) > nodeAlignmentTolerance)
    error('cfg.analysisTimeRange 的两个端点必须对齐到探针时间网格。');
end
analysisDuration = analysisTimeRange(2)-analysisTimeRange(1);
analysisPeriodCountFloat = analysisDuration/T;
analysisPeriodCount = round(analysisPeriodCountFloat);
if analysisPeriodCount < 2 || ...
        abs(analysisPeriodCountFloat-analysisPeriodCount) > ...
        1e-10*max(1,analysisPeriodCount)
    error('FFT 分析窗口必须覆盖至少两个完整整数调制周期。');
end

% 使用 [start,end) 的无重复端点序列。入口可保留 end 时刻做场图，但 FFT
% 不同时使用相邻周期的同一物理端点。
analysisRows = (analysisNodes(1)+1):analysisNodes(2);
nTime = numel(analysisRows);
if nTime ~= analysisPeriodCount*samplesPerPeriod
    error('内部采样计数与整数周期分析窗口不一致。');
end
analysisSignals = double(probeSignals(analysisRows,:,:));

%% ===================== FFT 参数与有限时间窗 =====================

if isfield(cfg,'zeroPaddingFactor') && ~isempty(cfg.zeroPaddingFactor)
    zeroPaddingFactor = read_integer(cfg.zeroPaddingFactor, ...
        'cfg.zeroPaddingFactor',1);
else
    zeroPaddingFactor = 1;
end
if isfield(cfg,'dynamicRangeDb') && ~isempty(cfg.dynamicRangeDb)
    dynamicRangeDb = read_positive_scalar(cfg.dynamicRangeDb, ...
        'cfg.dynamicRangeDb');
else
    dynamicRangeDb = 60;
end
if isfield(cfg,'activeColumnRelativeThreshold') && ...
        ~isempty(cfg.activeColumnRelativeThreshold)
    activeColumnRelativeThreshold = cfg.activeColumnRelativeThreshold;
    if ~isnumeric(activeColumnRelativeThreshold) || ...
            ~isscalar(activeColumnRelativeThreshold) || ...
            ~isreal(activeColumnRelativeThreshold) || ...
            ~isfinite(activeColumnRelativeThreshold) || ...
            activeColumnRelativeThreshold < 0 || ...
            activeColumnRelativeThreshold >= 1
        error('cfg.activeColumnRelativeThreshold 必须位于 [0,1)。');
    end
else
    activeColumnRelativeThreshold = 1e-12;
end
if isfield(cfg,'returnProbeSignals') && ~isempty(cfg.returnProbeSignals)
    returnProbeSignals = cfg.returnProbeSignals;
    if ~islogical(returnProbeSignals) || ~isscalar(returnProbeSignals)
        error('cfg.returnProbeSignals 必须是逻辑标量。');
    end
else
    returnProbeSignals = false;
end

sampleIndex = (0:nTime-1).';
window = 0.5-0.5*cos(2*pi*sampleIndex/nTime);
nFft = zeroPaddingFactor*nTime;
windowedSignals = analysisSignals.*reshape(window,[],1,1);

% 物理约定 exp(i*k*x-i*omega*t) 下，ifft 的正指数核把正 omega 放在
% 正频率坐标。多个固定探针不相干相加功率，避免空间节点造成相消。
% MATLAB 的 ifft 含 1/nFft。乘回 nFft 并除以窗函数总和，使一个
% 落在频率格点上的复指数峰值等于其相干振幅，且不随零填充改变。
% raw/folded power 因而是连续谱的离散采样；跨频率积分时仍应乘 bin 宽。
complexAmplitudeScale = nFft/sum(window);
complexSpectrum = ifft(windowedSignals,nFft,1)*complexAmplitudeScale;
rawPowerUnshifted = abs(complexSpectrum).^2;
rawPowerUnshifted = reshape(sum(rawPowerUnshifted,2),nFft,nK);

%% ===================== 把物理频率副本折叠到第一时间区 =====================

% nFold 是一个 Omega 宽度内的显示 bin 数。对未移位 FFT bin 按模 nFold
% 显式累加功率，适用于奇数/偶数周期数；不能只重排一个物理频率副本。
nFold = zeroPaddingFactor*analysisPeriodCount;
if mod(nFft,nFold) ~= 0 || nFft/nFold ~= samplesPerPeriod
    error('FFT 点数与每周期探针采样数不一致，无法严格执行 Floquet 折叠。');
end
foldedPowerUnshifted = zeros(nFold,nK);
for rawIndex = 1:nFft
    foldedIndex = mod(rawIndex-1,nFold)+1;
    foldedPowerUnshifted(foldedIndex,:) = ...
        foldedPowerUnshifted(foldedIndex,:)+rawPowerUnshifted(rawIndex,:);
end

frequencyOrders = (-floor(nFold/2):ceil(nFold/2)-1).';
centeredFoldedIndices = mod(frequencyOrders,nFold)+1;
foldedPower = foldedPowerUnshifted(centeredFoldedIndices,:);
omegaNormalized = frequencyOrders/nFold;

rawFrequencyOrders = (-floor(nFft/2):ceil(nFft/2)-1).';
centeredRawIndices = mod(rawFrequencyOrders,nFft)+1;
rawPower = rawPowerUnshifted(centeredRawIndices,:);
rawOmegaNormalized = rawFrequencyOrders* ...
    (T/(nFft*dtRecord));

%% ===================== 可靠列掩码与显示归一化 =====================

rawColumnPower = max(foldedPower,[],1);
globalMaximum = max(rawColumnPower);
if ~isfinite(globalMaximum) || globalMaximum <= 0
    error('探针 FFT 功率为空或发生溢出；请检查 FDTD 场与分析窗口。');
end
activeKMask = isfinite(rawColumnPower) & rawColumnPower > 0 & ...
    rawColumnPower >= activeColumnRelativeThreshold*globalMaximum;

relativePower = zeros(size(foldedPower));
relativePower(:,activeKMask) = foldedPower(:,activeKMask)./ ...
    rawColumnPower(activeKMask);
globalNormalizedPower = foldedPower/globalMaximum;
displayFloor = 10^(-dynamicRangeDb/10);
if displayFloor == 0
    error('cfg.dynamicRangeDb 过大，已超出双精度 dB 显示范围。');
end
spectrumDb = -dynamicRangeDb*ones(size(relativePower));
spectrumDb(:,activeKMask) = 10*log10(max( ...
    relativePower(:,activeKMask),displayFloor));

%% ===================== 可审计输出 =====================

result.kNormalized = kScan/Omega;
result.omegaNormalized = omegaNormalized;
result.spectrumDb = spectrumDb;
result.foldedPower = foldedPower;
result.rawColumnPower = rawColumnPower;
result.activeKMask = activeKMask;
result.columnNormalizedPower = relativePower;
result.globalNormalizedPower = globalNormalizedPower;
result.rawOmegaNormalized = rawOmegaNormalized;
result.rawPower = rawPower;
result.probeCount = nProbe;
result.recordTimeStep = dtRecord;
result.samplesPerPeriod = samplesPerPeriod;
result.analysisTimeRange = analysisTimeRange(:).';
result.analysisRows = analysisRows;
result.analysisPeriodCount = analysisPeriodCount;
result.temporalPeriod = T;
result.nFft = nFft;
result.nFoldedFrequency = nFold;
result.nativeOmegaResolutionNormalized = 1/analysisPeriodCount;
result.displayOmegaSpacingNormalized = 1/nFold;
result.windowName = 'periodic-hann';
result.windowCoherentGain = mean(window);
result.windowEnergyGain = mean(window.^2);
result.windowEnbwBins = mean(window.^2)/mean(window)^2;
result.window = struct('name','periodic-hann', ...
    'coherentGain',mean(window), ...
    'energyGain',mean(window.^2), ...
    'enbwBins',mean(window.^2)/mean(window)^2);
result.spectralNormalization = struct( ...
    'name','window-coherent-amplitude-squared', ...
    'complexAmplitudeScale',complexAmplitudeScale, ...
    'powerMeaning','squared coherent amplitude sampled on the FFT grid', ...
    'frequencyIntegralRequiresBinWidth',true);
result.activeColumnRelativeThreshold = activeColumnRelativeThreshold;
if returnProbeSignals
    result.probeSignals = probeSignals;
    result.time = time;
end
end

% -------------------------------------------------------------------------
function value = read_positive_scalar(value,label)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value <= 0
    error('%s 必须是有限正实标量。',label);
end
end

% -------------------------------------------------------------------------
function value = read_integer(value,label,minimumValue)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value ~= round(value) || value < minimumValue
    error('%s 必须是不小于 %d 的有限整数。',label,minimumValue);
end
end
