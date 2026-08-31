function result = reproduce_fdtd_band(overrides)
%REPRODUCE_FDTD_BAND  用有限空间高斯波包、固定 E 探针和 FFT 重建论文能带。
% 默认覆盖论文图 1 的 0<=k/k0<=5；可传入标量 overrides 结构体做网格、
% k 点数、时间周期数和波包宽度收敛。普通调用不需要任何输入。

if nargin < 1
    overrides = struct();
elseif ~isstruct(overrides) || ~isscalar(overrides)
    error('可选 overrides 必须是标量结构体。');
end
close all;

% 波包扫描、D/B-Yee 推进、探针 FFT 和 TMM 对照均强制解析到根目录。
[templateSource,templatePathCleanup] = use_root_templates( ...
    {'fdtd_gaussian_k_scan','fdtd1d','fdtd_fft_bands','tmm_bands'});
if isempty(templatePathCleanup)
    error('根目录模板路径清理对象创建失败。');
end

%% ===================== 论文参数与独立数值假设 =====================

paper = struct();
paper.c0 = 0.299792458;             % um/fs
paper.epsHigh = 3;
paper.epsLow = 1;
paper.muRelative = 1;
paper.temporalPeriod = 2;           % fs
paper.k0PerUm = 2*pi/(paper.temporalPeriod*paper.c0);
paper.figure2PulseFwhmFs = 45;      % 图 2 参数；论文没有指定能带扫描波包

T = paper.temporalPeriod;
Omega = 2*pi/T;

% 旧两状态回归可廉价使用 501 k × 64 周期；有限空间 FDTD 不能沿用同等
% 计算量。下面是新实验式数据流的默认验证配置，所有值都写入返回结构体。
numerics = struct();
numerics.kMinOverK0 = 0;
numerics.kMaxOverK0 = 5;
numerics.nK = read_override_integer(overrides,'nK',101,3);
numerics.stepsPerPeriod = read_override_integer( ...
    overrides,'stepsPerPeriod',100,20);
numerics.analysisPeriodCount = read_override_integer( ...
    overrides,'analysisPeriodCount',16,4);
numerics.zeroPaddingFactor = read_override_integer( ...
    overrides,'zeroPaddingFactor',4,1);
numerics.dxLightTimeFs = read_override_positive( ...
    overrides,'dxLightTimeFs',0.04);
numerics.sampleLengthLightTimeFs = read_override_positive( ...
    overrides,'sampleLengthLightTimeFs',192);
numerics.backgroundLengthLightTimeFs = read_override_nonnegative( ...
    overrides,'backgroundLengthLightTimeFs',4);
numerics.boundaryLengthLightTimeFs = read_override_nonnegative( ...
    overrides,'boundaryLengthLightTimeFs',4);
numerics.pulseIntensityFwhmFs = read_override_positive( ...
    overrides,'pulseIntensityFwhmFs',24);
numerics.probeOffsetsLightTimeFs = read_override_vector( ...
    overrides,'probeOffsetsLightTimeFs',[-8 -4 0 4 8]);
numerics.initialTailTolerance = read_override_fraction( ...
    overrides,'initialTailTolerance',1e-4);
numerics.displayDynamicRangeDb = read_override_positive( ...
    overrides,'displayDynamicRangeDb',60);
numerics.activeColumnRelativeThreshold = read_override_fraction_or_zero( ...
    overrides,'activeColumnRelativeThreshold',1e-12);
numerics.representativeKOverK0 = read_override_finite( ...
    overrides,'representativeKOverK0',0.65);
numerics.recordPrecision = 'single';
numerics.boundaryType = 'fixed-distant';

dx = numerics.dxLightTimeFs;
dt = T/numerics.stepsPerPeriod;
nSteps = numerics.analysisPeriodCount*numerics.stepsPerPeriod;
if dt/(dx*sqrt(min([paper.epsLow 1])*paper.muRelative)) >= 1
    error('当前 dx/dt 使最坏 Courant 数不小于 1。');
end

%% ===================== 方波时间材料与有限样品 =====================

highSteps = numerics.stepsPerPeriod/2;
if highSteps ~= round(highSteps)
    error('stepsPerPeriod 必须为偶数，使两个 1 fs 时间层严格对齐网格。');
end
phaseShift = T/4;
epsilonTime = @(time) paper.epsLow+(paper.epsHigh-paper.epsLow).* ...
    double(mod(time+phaseShift,T) < T/2);

candidateNodes = 0:nSteps-1;
candidateTimes = candidateNodes*dt;
epsilonBefore = epsilonTime(candidateTimes-dt/4);
epsilonAfter = epsilonTime(candidateTimes+dt/4);
interfaceTolerance = 128*eps(max([1 paper.epsHigh paper.epsLow]));
temporalInterfaces = candidateTimes( ...
    abs(epsilonAfter-epsilonBefore) > interfaceTolerance);

nSample = integer_grid_count( ...
    numerics.sampleLengthLightTimeFs,dx,'sampleLengthLightTimeFs');
nBackground = integer_grid_count( ...
    numerics.backgroundLengthLightTimeFs,dx, ...
    'backgroundLengthLightTimeFs');
nBoundary = integer_grid_count( ...
    numerics.boundaryLengthLightTimeFs,dx, ...
    'boundaryLengthLightTimeFs');
Nx = nSample+2*nBackground+2*nBoundary;
x = (0:Nx-1)*dx;                  % 光行时坐标 x=z/c0，单位 fs
iSampleFirst = nBoundary+nBackground+1;
iSampleLast = iSampleFirst+nSample-1;
sampleEdgesX = [x(iSampleFirst)-dx/2,x(iSampleLast)+dx/2];
xCenter = mean(sampleEdgesX);
inSample = @(position) position >= sampleEdgesX(1) & ...
    position < sampleEdgesX(2);
epsilonFunction = @(position,time) 1+ ...
    (epsilonTime(time)-1).*double(inSample(position));
muFunction = @(position,time) paper.muRelative*ones(size(position));

probeTargets = xCenter+numerics.probeOffsetsLightTimeFs;
if any(probeTargets <= sampleEdgesX(1)) || ...
        any(probeTargets >= sampleEdgesX(2))
    error('所有能带探针都必须严格位于有限时间晶体样品内部。');
end
probeIndices = zeros(size(probeTargets));
for probeIndex = 1:numel(probeTargets)
    [~,probeIndices(probeIndex)] = min(abs(x-probeTargets(probeIndex)));
end
if any(diff(probeIndices) <= 0)
    error('探针对齐网格后必须保持严格递增且互不重复。');
end

amplitudeGaussianWidth = ...
    numerics.pulseIntensityFwhmFs/sqrt(2*log(2));
sampleEdgeAmplitude = exp(-(numerics.sampleLengthLightTimeFs/2 / ...
    amplitudeGaussianWidth)^2);
if sampleEdgeAmplitude > numerics.initialTailTolerance
    error('高斯初值在有限样品端面仍高于 initialTailTolerance。');
end
maximumSpeed = 1/sqrt(min([paper.epsLow 1])*paper.muRelative);
resolvableTailRadius = amplitudeGaussianWidth* ...
    sqrt(log(1/numerics.initialTailTolerance));
distanceToBoundary = min(xCenter-x(1),x(end)-xCenter);
if maximumSpeed*nSteps*dt+resolvableTailRadius >= distanceToBoundary
    error(['固定远端距离不足：主场或指定阈值以上的初始高斯尾场可能在' ...
        '分析结束前到达端点。请增大空间长度或缩短记录周期数。']);
end

%% ===================== 扫描 k 并生成真实探针 E(t) =====================

kOverK0 = linspace(numerics.kMinOverK0,numerics.kMaxOverK0,numerics.nK);
kPhysical = kOverK0*paper.k0PerUm;
kSolver = paper.c0*kPhysical;     % = (k/k0)*Omega，单位 1/fs

scanCfg = struct();
scanCfg.x = x;
scanCfg.dt = dt;
scanCfg.nSteps = nSteps;
scanCfg.epsFun = epsilonFunction;
scanCfg.muFun = muFunction;
scanCfg.kScan = kSolver;
scanCfg.xCenter = xCenter;
scanCfg.pulseIntensityFwhm = numerics.pulseIntensityFwhmFs;
scanCfg.fieldAmplitude = 1;
scanCfg.initialEpsilon = epsilonTime(-dt/2);
scanCfg.muRelative = paper.muRelative;
scanCfg.probeIndices = probeIndices;
scanCfg.representativeSpatialIndices = iSampleFirst:iSampleLast;
scanCfg.representativeK = numerics.representativeKOverK0*Omega;
scanCfg.recordEvery = 1;
scanCfg.recordPrecision = numerics.recordPrecision;
scanCfg.boundaryType = numerics.boundaryType;
scanCfg.spongeCells = 0;
scanCfg.spongeStrength = 0;
scanCfg.temporalInterfaces = temporalInterfaces;
scanCfg.certifiedMinimumEpsilon = min([paper.epsLow 1]);
scanCfg.certifiedMinimumMu = paper.muRelative;
scanCfg.progressEvery = max(1,ceil(numerics.nK/10));
scanCfg.progressLabel = '论文参数有限空间 FDTD';
scanCfg.progressKScale = Omega;
scanCfg.progressKLabel = 'k/k0';
scanResult = fdtd_gaussian_k_scan(scanCfg);

%% ===================== 新探针接口 FFT 与 Floquet 折叠 =====================

fftCfg = struct();
fftCfg.temporalPeriod = T;
fftCfg.analysisTimeRange = [0 numerics.analysisPeriodCount*T];
fftCfg.zeroPaddingFactor = numerics.zeroPaddingFactor;
fftCfg.dynamicRangeDb = numerics.displayDynamicRangeDb;
fftCfg.activeColumnRelativeThreshold = ...
    numerics.activeColumnRelativeThreshold;
fftCfg.returnProbeSignals = false;
fftBands = fdtd_fft_bands( ...
    scanResult.probeSignals,scanResult.time,kSolver,fftCfg);

%% ===================== 第一 Floquet 区响应谱 =====================

figure('Color','w','Position',[100 80 900 650]);
imagesc(kOverK0,2*pi*fftBands.omegaNormalized,fftBands.spectrumDb);
axis xy;
bandAxes = gca;
bandAxes.CLim = [-numerics.displayDynamicRangeDb 0];
colormap(parula(256));
bandColorbar = colorbar;
bandColorbar.Label.String = '每个高斯激发 k_c 独立归一化的 E 探针 FFT 功率 (dB)';
xlabel('高斯激发中心 k_c/k_0');
ylabel('\omega_F T');
xlim([numerics.kMinOverK0 numerics.kMaxOverK0]);
ylim([-pi pi]);
title(sprintf(['论文参数的有限空间波包/探针 FDTD–FFT' ...
    '（N_k=%d，N_T=%d）'],numerics.nK,numerics.analysisPeriodCount));
box on;

%% ===================== 对 Yee-TMM 的有限窗峰位核对 =====================

% 用中心 k 的 Yee 离散波数隔离空间离散误差；理论曲线不叠加到主图。
kYee = (2/dx)*sin(kSolver*dx/2);
tmmYee = tmm_bands(kYee,[paper.epsHigh paper.epsLow], ...
    paper.muRelative,[T/2 T/2]);
tmmNormalized = tmmYee.omega/Omega;
stableMask = max(abs(imag(tmmNormalized)),[],1) < 1e-9;
searchHalfWidth = max(0.12, ...
    1.5*fftBands.nativeOmegaResolutionNormalized);
[visibleMaximumError,visibleMedianError,visibleSampleCount, ...
    minimumVisiblePeakDb,visibleMatches] = ...
    compare_stable_ridges(fftBands.omegaNormalized,fftBands.spectrumDb, ...
    tmmNormalized,stableMask,searchHalfWidth,-35);
[strongMaximumError,strongMedianError,strongSampleCount, ...
    minimumStrongPeakDb,strongMatches] = ...
    compare_stable_ridges(fftBands.omegaNormalized,fftBands.spectrumDb, ...
    tmmNormalized,stableMask,searchHalfWidth,-20);
strongP90Error = empirical_percentile( ...
    strongMatches.errorOverOmega,0.90);
visibleP90Error = empirical_percentile( ...
    visibleMatches.errorOverOmega,0.90);

fprintf(['有限空间高斯/探针 FDTD–FFT：%d 个 k，%d 个探针，' ...
    '%d 个周期，Courant=%.6g。\n'], ...
    numerics.nK,numel(probeIndices),numerics.analysisPeriodCount, ...
    scanResult.courant);
fprintf(['对 Yee-TMM 稳定通带的强峰（>=-20 dB）：样本数 %d，' ...
    'median/P90/max |Delta omega|/Omega = %.6g / %.6g / %.6g。\n'], ...
    strongSampleCount,strongMedianError,strongP90Error,strongMaximumError);
fprintf(['全部可见峰（>=-35 dB）：样本数 %d，median/max ' ...
    '|Delta omega|/Omega = %.6g / %.6g，最弱匹配峰 %.3f dB。\n'], ...
    visibleSampleCount,visibleMedianError,visibleMaximumError, ...
    minimumVisiblePeakDb);

%% ===================== 返回复现与实验式观测数据 =====================

representativeField = scanResult.representativeField;
representativeField.zUm = paper.c0*representativeField.x;
representativeField.kOverK0 = representativeField.k/Omega;

result.kNormalized = kOverK0;
result.kPhysicalPerUm = kPhysical;
result.omegaT = 2*pi*fftBands.omegaNormalized;
result.spectrumDb = fftBands.spectrumDb;
result.foldedPower = fftBands.foldedPower;
result.rawOmegaNormalized = fftBands.rawOmegaNormalized;
result.rawPower = fftBands.rawPower;
result.rawColumnPower = fftBands.rawColumnPower;
result.activeKMask = fftBands.activeKMask;
result.probeSignals = scanResult.probeSignals;
result.probeTimeFs = scanResult.time;
result.probePositionsUm = paper.c0*scanResult.probePositions;
result.representativeField = representativeField;
result.validation = struct( ...
    'stableStrongSampleCount',strongSampleCount, ...
    'strongMedianErrorOverOmega',strongMedianError, ...
    'strongP90ErrorOverOmega',strongP90Error, ...
    'strongMaximumErrorOverOmega',strongMaximumError, ...
    'minimumStrongPeakDb',minimumStrongPeakDb, ...
    'stableVisibleSampleCount',visibleSampleCount, ...
    'visibleMedianErrorOverOmega',visibleMedianError, ...
    'visibleP90ErrorOverOmega',visibleP90Error, ...
    'visibleMaximumErrorOverOmega',visibleMaximumError, ...
    'minimumVisiblePeakDb',minimumVisiblePeakDb, ...
    'searchHalfWidthOverOmega',searchHalfWidth, ...
    'strongMatches',strongMatches,'visibleMatches',visibleMatches);
result.paperParameters = paper;
result.numericalParameters = numerics;
result.numericalParameters.dtFs = dt;
result.numericalParameters.dxPhysicalUm = paper.c0*dx;
result.numericalParameters.sampleEdgesUm = paper.c0*sampleEdgesX;
result.numericalParameters.probePositionsUm = ...
    paper.c0*scanResult.probePositions;
result.samplingMetadata = struct('nFft',fftBands.nFft, ...
    'nFoldedFrequency',fftBands.nFoldedFrequency, ...
    'nativeOmegaResolutionNormalized', ...
    fftBands.nativeOmegaResolutionNormalized, ...
    'displayOmegaSpacingNormalized', ...
    fftBands.displayOmegaSpacingNormalized, ...
    'windowName',fftBands.windowName, ...
    'windowEnbwBins',fftBands.windowEnbwBins, ...
    'courantNumber',scanResult.courant);
result.templateSource = templateSource;
end

% -------------------------------------------------------------------------
function [maximumError,medianError,sampleCount,minimumPeakDb,matches] = ...
        compare_stable_ridges(omegaAxis,spectrumDb,theoryOmega,stableMask, ...
        searchHalfWidth,visibilityFloorDb)
% 仅统计理论邻域中确实高于可见阈值的有限波包响应峰。
errors = [];
peakLevels = [];
kIndices = [];
branches = [];
targetFrequencies = [];
observedFrequencies = [];
for kIndex = find(stableMask)
    spectrumColumn = spectrumDb(:,kIndex);
    previousValue = spectrumColumn([end;(1:numel(spectrumColumn)-1).']);
    nextValue = spectrumColumn([(2:numel(spectrumColumn)).';1]);
    localMaximumMask = spectrumColumn >= previousValue & ...
        spectrumColumn >= nextValue;
    for branch = 1:size(theoryOmega,1)
        target = real(theoryOmega(branch,kIndex));
        circularOffset = mod(omegaAxis-target+0.5,1)-0.5;
        candidateRows = find(abs(circularOffset) <= searchHalfWidth & ...
            localMaximumMask);
        if isempty(candidateRows)
            continue;
        end
        [peakLevel,localMaximum] = max(spectrumDb(candidateRows,kIndex));
        if peakLevel < visibilityFloorDb
            continue;
        end
        peakRow = candidateRows(localMaximum);
        errors(end+1) = abs(circularOffset(peakRow)); %#ok<AGROW>
        peakLevels(end+1) = peakLevel; %#ok<AGROW>
        kIndices(end+1) = kIndex; %#ok<AGROW>
        branches(end+1) = branch; %#ok<AGROW>
        targetFrequencies(end+1) = target; %#ok<AGROW>
        observedFrequencies(end+1) = omegaAxis(peakRow); %#ok<AGROW>
    end
end
if isempty(errors)
    error('没有找到高于 %.3g dB 的稳定通带探针峰。',visibilityFloorDb);
end
maximumError = max(errors);
medianError = median(errors);
sampleCount = numel(errors);
minimumPeakDb = min(peakLevels);
matches = struct('errorOverOmega',errors,'peakDb',peakLevels, ...
    'kIndex',kIndices,'branch',branches, ...
    'targetOmegaNormalized',targetFrequencies, ...
    'observedOmegaNormalized',observedFrequencies);
end

% -------------------------------------------------------------------------
function value = empirical_percentile(samples,fraction)
% Base MATLAB 的确定性最近秩分位数，不依赖 Statistics Toolbox。
samples = sort(samples(:));
index = max(1,min(numel(samples),ceil(fraction*numel(samples))));
value = samples(index);
end

% -------------------------------------------------------------------------
function count = integer_grid_count(lengthValue,dx,label)
countFloat = lengthValue/dx;
count = round(countFloat);
if count < 0 || abs(countFloat-count) > 1e-10*max(1,count)
    error('%s 必须是 dx 的非负整数倍。',label);
end
end

% -------------------------------------------------------------------------
function value = read_override_integer(overrides,name,defaultValue,minimumValue)
if isfield(overrides,name) && ~isempty(overrides.(name))
    value = overrides.(name);
else
    value = defaultValue;
end
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value ~= round(value) || value < minimumValue
    error('overrides.%s 必须是不小于 %d 的有限整数。',name,minimumValue);
end
end

% -------------------------------------------------------------------------
function value = read_override_positive(overrides,name,defaultValue)
value = read_override_finite(overrides,name,defaultValue);
if value <= 0
    error('overrides.%s 必须是有限正实数。',name);
end
end

% -------------------------------------------------------------------------
function value = read_override_nonnegative(overrides,name,defaultValue)
value = read_override_finite(overrides,name,defaultValue);
if value < 0
    error('overrides.%s 必须是有限非负实数。',name);
end
end

% -------------------------------------------------------------------------
function value = read_override_fraction(overrides,name,defaultValue)
value = read_override_finite(overrides,name,defaultValue);
if value <= 0 || value >= 1
    error('overrides.%s 必须位于 (0,1)。',name);
end
end

% -------------------------------------------------------------------------
function value = read_override_fraction_or_zero(overrides,name,defaultValue)
value = read_override_finite(overrides,name,defaultValue);
if value < 0 || value >= 1
    error('overrides.%s 必须位于 [0,1)。',name);
end
end

% -------------------------------------------------------------------------
function value = read_override_finite(overrides,name,defaultValue)
if isfield(overrides,name) && ~isempty(overrides.(name))
    value = overrides.(name);
else
    value = defaultValue;
end
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value)
    error('overrides.%s 必须是有限实标量。',name);
end
end

% -------------------------------------------------------------------------
function value = read_override_vector(overrides,name,defaultValue)
if isfield(overrides,name) && ~isempty(overrides.(name))
    value = overrides.(name);
else
    value = defaultValue;
end
if ~isnumeric(value) || ~isvector(value) || isempty(value) || ...
        ~isreal(value) || any(~isfinite(value)) || any(diff(value(:)) <= 0)
    error('overrides.%s 必须是严格递增的非空有限实向量。',name);
end
value = value(:).';
end
