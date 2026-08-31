%RUN_FDTD_FFT  扫描高斯波包中心 k，并由有限空间 FDTD 探针数据重建能带。
% 每个 k 都构造真实的有限宽度波包、调用 fdtd1d 推进，并在固定空间探针
% 记录 E(t)。fdtd_fft_bands 只处理这些可测时间序列，不再推进单一 Fourier 模式。

clear; clc; close all;

%% ===================== 用户参数：时间晶体材料 =====================

modulationType = 'square';       % 'square' 或 'sinusoidal'
epsHigh = 4;
epsLow = 1;
muRelative = 1;
backgroundEpsilon = 1;

temporalPeriod = 1;
dutyCycle = 0.5;                 % 仅用于方波
squarePhaseFraction = 0.25;      % t=0 在方波周期中的相位，范围 [0,1)
sinePhase = 0;                   % 仅用于正弦调制，单位 rad

%% ===================== 用户参数：扫描与高斯初始激发 =====================

% kNormalized = k/Omega（归一化单位 c0=1）。每一列能带数据来自一次独立
% 的有限空间 FDTD；高斯包的有限谱宽意味着结果是源加权响应，不是本征值。
kMinNormalized = -1.1;
kMaxNormalized = 1.1;
nK = 81;

fieldAmplitude = 1;
pulseIntensityFwhm = 10;         % 强度 |E|^2 的空间 FWHM

%% ===================== 用户参数：有限样品、时间与探针 =====================

% 纯时间晶体没有空间晶格。cellSize 只作为方便定义有限几何的长度单位，
% sampleCellCount*cellSize 是实际有限样品长度。
cellSize = 1;
sampleCellCount = 72;
pointsPerCell = 12;
backgroundCellCount = 6;
boundaryLayerCellCount = 4;

stepsPerPeriod = 100;
simulationPeriodCount = 16;
recordEvery = 1;

boundaryType = 'fixed-distant';  % 'fixed-distant' 或 'sponge'
spongeStrength = 8;              % 仅用于 sponge；该边界不是 PML
initialTailTolerance = 1e-4;     % 边界安全检查采用的相对振幅阈值

% “特殊点”相对波包/样品中心的位置。各探针先分别 FFT，再非相干相加功率，
% 以免单个点恰好落在场节点而漏掉能带分支。
probeOffsets = [-4 -2 0 2 4]*cellSize;

% 保存最接近该 k 的一次完整样品内部 E(x,t)，用于核查真实时空场。
representativeKNormalized = 0.70;
showRepresentativeFieldMap = true;
recordPrecision = 'single';      % 只影响历史数组；内部推进保持 double
progressPrintEvery = 10;

%% ===================== 用户参数：FFT 重建 =====================

analysisStartPeriod = 0;
analysisPeriodCount = simulationPeriodCount;
zeroPaddingFactor = 4;           % 只细化显示网格，不提高真实分辨率
displayDynamicRangeDb = 60;
activeColumnRelativeThreshold = 1e-12;
returnProbeSignalsInResult = false;

%% ===================== 参数检查与材料函数 =====================

modulationType = validate_choice(modulationType, ...
    {'square','sinusoidal'},'modulationType');
boundaryType = validate_choice(boundaryType, ...
    {'fixed-distant','sponge'},'boundaryType');
recordPrecision = validate_choice(recordPrecision, ...
    {'double','single'},'recordPrecision');

positiveParameters = {epsHigh,epsLow,muRelative,backgroundEpsilon, ...
    temporalPeriod,cellSize,pulseIntensityFwhm,fieldAmplitude, ...
    initialTailTolerance,displayDynamicRangeDb};
validPositive = cellfun(@(value) isnumeric(value) && isscalar(value) && ...
    isreal(value) && isfinite(value) && value > 0,positiveParameters);
if ~all(validPositive) || epsHigh < epsLow || initialTailTolerance >= 1
    error(['介电常数、mu、周期、长度、场幅和动态范围必须为有限正数；' ...
        '同时要求 epsHigh>=epsLow 且 initialTailTolerance<1。']);
end

integerParameters = {nK,sampleCellCount,pointsPerCell, ...
    backgroundCellCount,boundaryLayerCellCount,stepsPerPeriod, ...
    simulationPeriodCount,recordEvery,analysisStartPeriod, ...
    analysisPeriodCount,zeroPaddingFactor,progressPrintEvery};
validIntegers = cellfun(@(value) isnumeric(value) && isscalar(value) && ...
    isreal(value) && isfinite(value) && value == round(value),integerParameters);
if ~all(validIntegers) || nK < 2 || sampleCellCount < 1 || ...
        pointsPerCell < 2 || backgroundCellCount < 0 || ...
        boundaryLayerCellCount < 0 || stepsPerPeriod < 4 || ...
        simulationPeriodCount < 1 || recordEvery < 1 || ...
        analysisStartPeriod < 0 || analysisPeriodCount < 2 || ...
        zeroPaddingFactor < 1 || progressPrintEvery < 1
    error('空间、时间、扫描、记录与 FFT 计数参数不是允许范围内的整数。');
end
if mod(stepsPerPeriod,recordEvery) ~= 0
    error('recordEvery 必须整除 stepsPerPeriod，保证探针采样与调制周期可公度。');
end
if analysisStartPeriod+analysisPeriodCount > simulationPeriodCount
    error('FFT 分析窗口不能超出实际 FDTD 仿真时间。');
end
if ~isnumeric(kMinNormalized) || ~isscalar(kMinNormalized) || ...
        ~isreal(kMinNormalized) || ~isfinite(kMinNormalized) || ...
        ~isnumeric(kMaxNormalized) || ~isscalar(kMaxNormalized) || ...
        ~isreal(kMaxNormalized) || ~isfinite(kMaxNormalized) || ...
        kMinNormalized >= kMaxNormalized
    error('kMinNormalized 和 kMaxNormalized 必须是递增的有限实标量。');
end
if ~isnumeric(probeOffsets) || ~isvector(probeOffsets) || ...
        isempty(probeOffsets) || ~isreal(probeOffsets) || ...
        any(~isfinite(probeOffsets)) || any(diff(probeOffsets(:)) <= 0)
    error('probeOffsets 必须是严格递增的非空有限实向量。');
end
if ~isnumeric(representativeKNormalized) || ...
        ~isscalar(representativeKNormalized) || ...
        ~isreal(representativeKNormalized) || ...
        ~isfinite(representativeKNormalized)
    error('representativeKNormalized 必须是有限实标量。');
end
if ~islogical(showRepresentativeFieldMap) || ...
        ~isscalar(showRepresentativeFieldMap) || ...
        ~islogical(returnProbeSignalsInResult) || ...
        ~isscalar(returnProbeSignalsInResult)
    error('两个显示/返回开关必须是逻辑标量。');
end
if ~isnumeric(activeColumnRelativeThreshold) || ...
        ~isscalar(activeColumnRelativeThreshold) || ...
        ~isreal(activeColumnRelativeThreshold) || ...
        ~isfinite(activeColumnRelativeThreshold) || ...
        activeColumnRelativeThreshold < 0 || activeColumnRelativeThreshold >= 1
    error('activeColumnRelativeThreshold 必须位于 [0,1)。');
end

T = temporalPeriod;
Omega = 2*pi/T;
dt = T/stepsPerPeriod;
nSteps = simulationPeriodCount*stepsPerPeriod;

switch modulationType
    case 'square'
        highStepsFloat = dutyCycle*stepsPerPeriod;
        phaseStepsFloat = squarePhaseFraction*stepsPerPeriod;
        if ~isnumeric(dutyCycle) || ~isscalar(dutyCycle) || ...
                ~isreal(dutyCycle) || ~isfinite(dutyCycle) || ...
                dutyCycle <= 0 || dutyCycle >= 1 || ...
                abs(highStepsFloat-round(highStepsFloat)) > 1e-10
            error('dutyCycle 位于 (0,1)，且 dutyCycle*stepsPerPeriod 必须为整数。');
        end
        if ~isnumeric(squarePhaseFraction) || ...
                ~isscalar(squarePhaseFraction) || ...
                ~isreal(squarePhaseFraction) || ...
                ~isfinite(squarePhaseFraction) || ...
                squarePhaseFraction < 0 || squarePhaseFraction >= 1 || ...
                abs(phaseStepsFloat-round(phaseStepsFloat)) > 1e-10
            error(['squarePhaseFraction 必须位于 [0,1)，并对应整数个' ...
                ' FDTD 时间步。']);
        end
        if round(highStepsFloat) < 2 || ...
                stepsPerPeriod-round(highStepsFloat) < 2
            error('方波高、低介电时间层均至少需要两个 FDTD 步。');
        end
    case 'sinusoidal'
        if ~isnumeric(sinePhase) || ~isscalar(sinePhase) || ...
                ~isreal(sinePhase) || ~isfinite(sinePhase)
            error('sinePhase 必须是有限实标量。');
        end
end

epsilonTime = @(time) periodic_temporal_epsilon(time,modulationType, ...
    epsHigh,epsLow,T,dutyCycle,squarePhaseFraction,sinePhase);

% 自动识别所有方波时间界面，供 fdtd1d 的 D/B 中心修正使用。
candidateNodes = 0:nSteps-1;
candidateTimes = candidateNodes*dt;
epsilonBeforeNode = epsilonTime(candidateTimes-dt/4);
epsilonAfterNode = epsilonTime(candidateTimes+dt/4);
interfaceTolerance = 128*eps(max([1 epsHigh epsLow]));
temporalInterfaceNodes = candidateNodes( ...
    abs(epsilonAfterNode-epsilonBeforeNode) > interfaceTolerance);
temporalInterfaces = temporalInterfaceNodes*dt;

%% ===================== 有限空间样品与测量点 =====================

dx = cellSize/pointsPerCell;
nSample = sampleCellCount*pointsPerCell;
nBackground = backgroundCellCount*pointsPerCell;
nBoundary = boundaryLayerCellCount*pointsPerCell;
Nx = nSample+2*nBackground+2*nBoundary;
if Nx < 5
    error('有限空间网格至少需要 5 个 E 节点。');
end
x = (0:Nx-1)*dx;
iSampleFirst = nBoundary+nBackground+1;
iSampleLast = iSampleFirst+nSample-1;
sampleEdges = [x(iSampleFirst)-dx/2,x(iSampleLast)+dx/2];
xCenter = mean(sampleEdges);
inSample = @(position) position >= sampleEdges(1) & ...
    position < sampleEdges(2);

epsilonFunction = @(position,time) backgroundEpsilon+ ...
    (epsilonTime(time)-backgroundEpsilon).*double(inSample(position));
muFunction = @(position,time) muRelative*ones(size(position));

probeTargets = xCenter+probeOffsets(:).';
if any(probeTargets <= sampleEdges(1)) || any(probeTargets >= sampleEdges(2))
    error('所有 probeOffsets 对应的测点都必须严格位于有限样品内部。');
end
probeIndices = zeros(size(probeTargets));
for probeIndex = 1:numel(probeTargets)
    [~,probeIndices(probeIndex)] = min(abs(x-probeTargets(probeIndex)));
end
if numel(unique(probeIndices)) ~= numel(probeIndices)
    error('多个 probeOffsets 落在同一 Yee 节点；请增大探针间距或细化网格。');
end
if any(probeIndices < iSampleFirst) || any(probeIndices > iSampleLast)
    error('探针对齐到 Yee 网格后越出样品，请调整 probeOffsets。');
end

minimumEpsilon = min([epsLow epsHigh backgroundEpsilon]);
maximumWaveSpeed = 1/sqrt(minimumEpsilon*muRelative);
amplitudeGaussianWidth = pulseIntensityFwhm/sqrt(2*log(2));
sampleEdgeAmplitude = exp(-((sampleCellCount*cellSize/2) / ...
    amplitudeGaussianWidth)^2);
if sampleEdgeAmplitude > initialTailTolerance
    error(['初始高斯波包在样品端面仍大于 initialTailTolerance；' ...
        '请增大样品或减小 pulseIntensityFwhm。']);
end

if strcmp(boundaryType,'fixed-distant')
    distanceToBoundary = min(xCenter-x(1),x(end)-xCenter);
    resolvableTailRadius = amplitudeGaussianWidth* ...
        sqrt(log(1/initialTailTolerance));
    if maximumWaveSpeed*nSteps*dt+resolvableTailRadius >= distanceToBoundary
        error(['fixed-distant 边界距离不足：在指定尾场阈值下，场可能在 FFT ' ...
            '窗口结束前到达端点。请增加背景/边界长度、缩短仿真时间，' ...
            '或显式选择 sponge 并做边界收敛。']);
    end
    spongeCellsForSolver = 0;
    spongeStrengthForSolver = 0;
else
    if nBoundary < 1 || ~isnumeric(spongeStrength) || ...
            ~isscalar(spongeStrength) || ~isreal(spongeStrength) || ...
            ~isfinite(spongeStrength) || spongeStrength <= 0
        error('sponge 边界要求正的边界层宽度和 spongeStrength。');
    end
    warning(['当前 sponge 不是 PML；请改变边界距离、宽度和强度检查' ...
        '能带脊线的边界收敛。']);
    spongeCellsForSolver = nBoundary;
    spongeStrengthForSolver = spongeStrength;
end

%% ===================== 逐 k 构造波包并运行完整 FDTD =====================

kNormalized = linspace(kMinNormalized,kMaxNormalized,nK);
kScan = kNormalized*Omega;

scanCfg = struct();
scanCfg.x = x;
scanCfg.dt = dt;
scanCfg.nSteps = nSteps;
scanCfg.epsFun = epsilonFunction;
scanCfg.muFun = muFunction;
scanCfg.kScan = kScan;
scanCfg.xCenter = xCenter;
scanCfg.pulseIntensityFwhm = pulseIntensityFwhm;
scanCfg.fieldAmplitude = fieldAmplitude;
scanCfg.initialEpsilon = epsilonTime(-dt/2);
scanCfg.muRelative = muRelative;
scanCfg.probeIndices = probeIndices;
scanCfg.representativeSpatialIndices = iSampleFirst:iSampleLast;
scanCfg.representativeK = representativeKNormalized*Omega;
scanCfg.recordEvery = recordEvery;
scanCfg.recordPrecision = recordPrecision;
scanCfg.boundaryType = boundaryType;
scanCfg.spongeCells = spongeCellsForSolver;
scanCfg.spongeStrength = spongeStrengthForSolver;
scanCfg.temporalInterfaces = temporalInterfaces;
scanCfg.certifiedMinimumEpsilon = minimumEpsilon;
scanCfg.certifiedMinimumMu = muRelative;
scanCfg.progressEvery = progressPrintEvery;
scanCfg.progressLabel = 'FDTD k 扫描';
scanCfg.progressKScale = Omega;
scanCfg.progressKLabel = 'k/Omega';

scanResult = fdtd_gaussian_k_scan(scanCfg);
probeSignals = scanResult.probeSignals;
probeTime = scanResult.time;
probePositions = scanResult.probePositions;
representativeField = scanResult.representativeField;
representativeField.kNormalized = representativeField.k/Omega;
representativeField.sampleEdges = sampleEdges;

%% ===================== 探针 FFT 与 Floquet 折叠 =====================

fftCfg = struct();
fftCfg.temporalPeriod = T;
fftCfg.analysisTimeRange = ...
    [analysisStartPeriod,analysisStartPeriod+analysisPeriodCount]*T;
fftCfg.zeroPaddingFactor = zeroPaddingFactor;
fftCfg.dynamicRangeDb = displayDynamicRangeDb;
fftCfg.activeColumnRelativeThreshold = activeColumnRelativeThreshold;
fftCfg.returnProbeSignals = returnProbeSignalsInResult;
fftBands = fdtd_fft_bands(probeSignals,probeTime,kScan,fftCfg);
fftBands.probePositions = probePositions;
fftBands.probeOffsets = probePositions-xCenter;
fftBands.sampleEdges = sampleEdges;
fftBands.pulseIntensityFwhm = pulseIntensityFwhm;
fftBands.representativeKNormalized = representativeField.kNormalized;

%% ===================== 绘制实测式场图与重建能带 =====================

if showRepresentativeFieldMap
    figure('Color','w','Position',[80 80 980 600]);
    imagesc(representativeField.x,representativeField.t/T, ...
        abs(representativeField.E));
    axis xy;
    colormap(parula(256));
    fieldColorbar = colorbar;
    fieldColorbar.Label.String = '|E(x,t)|';
    hold on;
    for probePosition = probePositions
        xline(probePosition,'w:','LineWidth',0.9);
    end
    xlabel('x');
    ylabel('归一化时间 t/T');
    title(sprintf('代表性高斯波包 FDTD 场：k_c/\\Omega=%.4g', ...
        representativeField.kNormalized));
    box on;
end

figure('Color','w','Position',[110 90 900 650]);
imagesc(fftBands.kNormalized,fftBands.omegaNormalized,fftBands.spectrumDb);
axis xy;
colormap(parula(256));
bandAxes = gca;
bandAxes.CLim = [-displayDynamicRangeDb 0];
bandColorbar = colorbar;
bandColorbar.Label.String = '每个激发 k 独立归一化的探针 FFT 功率 (dB)';
xlabel('高斯激发中心 k_c/\Omega');
ylabel('Re(\omega_F)/\Omega');
ylim([-0.5 0.5]);
title(sprintf(['有限空间高斯激发/探针采样重建能带' ...
    '（N_T=%d，N_p=%d）'],analysisPeriodCount,numel(probeIndices)));
box on;

fprintf(['FDTD–FFT 完成：%d 次有限空间仿真，%d 个 E 探针，' ...
    'Courant=%.6g，原生 Delta(omega_F/Omega)=%.6g；' ...
    '图为有限波包和有限时间窗的源加权响应。\n'], ...
    nK,numel(probeIndices),representativeField.courant, ...
    fftBands.nativeOmegaResolutionNormalized);

% -------------------------------------------------------------------------
function value = validate_choice(value,allowedValues,label)
% 把字符向量或标量 string 统一为允许列表中的标准字符文本。
if isstring(value) && isscalar(value)
    value = char(value);
elseif ~ischar(value) || size(value,1) ~= 1
    error('%s 必须是文本标量。',label);
end
matched = strcmpi(value,allowedValues);
if ~any(matched)
    error('%s 必须是：%s。',label,strjoin(allowedValues,', '));
end
value = allowedValues{find(matched,1)};
end

% -------------------------------------------------------------------------
function epsilon = periodic_temporal_epsilon(time,modulationType, ...
        epsHigh,epsLow,T,dutyCycle,squarePhaseFraction,sinePhase)
% 周期材料在仿真开始前已经存在，避免把泵浦开启瞬态混入能带 FFT。
switch modulationType
    case 'square'
        phase = mod(time/T+squarePhaseFraction,1);
        epsilon = epsLow+(epsHigh-epsLow).*double(phase < dutyCycle);
    case 'sinusoidal'
        epsMean = (epsHigh+epsLow)/2;
        epsAmplitude = (epsHigh-epsLow)/2;
        epsilon = epsMean+epsAmplitude*cos(2*pi*time/T+sinePhase);
end
end
