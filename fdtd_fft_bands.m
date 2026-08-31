function result = fdtd_fft_bands(probeSignals,time,kScan,cfg)
%FDTD_FFT_BANDS  由逐 k 有限空间 FDTD 的固定探针 E(t) 重建 Floquet 能带。
% probeSignals 的尺寸为 Nt-by-Nprobe-by-Nk；二维 Nt-by-Nk 输入表示单探针。
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
if max(abs(dtVector-dtRecord)) > 1e-10*max(1,abs(dtRecord))
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
    if size(probeSignals,1) ~= numel(time) || size(probeSignals,2) ~= nK
        error('二维 probeSignals 必须具有 Nt-by-Nk 尺寸。');
    end
    probeSignals = reshape(probeSignals,numel(time),1,nK);
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
        analysisTimeRange(1) < time(1) || ...
        analysisTimeRange(2) <= analysisTimeRange(1) || ...
        analysisTimeRange(2) > time(end)+ ...
        1e-10*max(1,abs(time(end)))
    error('cfg.analysisTimeRange 必须是记录范围内严格递增的两个有限时刻。');
end
analysisNodeFloat = (analysisTimeRange-time(1))/dtRecord;
analysisNodes = round(analysisNodeFloat);
if any(abs(analysisNodeFloat-analysisNodes) > ...
        1e-10*max(1,numel(time)))
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
complexSpectrum = ifft(windowedSignals,nFft,1);
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
result.rawOmegaNormalized = rawOmegaNormalized;
result.rawPower = rawPower;
result.probeCount = nProbe;
result.recordTimeStep = dtRecord;
result.samplesPerPeriod = samplesPerPeriod;
result.analysisTimeRange = analysisTimeRange(:).';
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
result.activeColumnRelativeThreshold = activeColumnRelativeThreshold;
if returnProbeSignals
    result.probeSignals = probeSignals;
    result.time = time;
    result.analysisRows = analysisRows;
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
