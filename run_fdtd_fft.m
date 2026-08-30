%RUN_FDTD_FFT  逐 k 进行有限时间原胞采样，并用 FFT 重建 Floquet 能带。
% 本脚本独立运行，只绘制有限采样得到的第一时间 Floquet 区能带映射。

clear; clc; close all;

%% ===================== 用户参数：时间晶体材料 =====================

% 与场分布模块相同，介电调制可选方波或正弦；样品在空间上保持均匀。
modulationType = 'square';
epsHigh = 4;
epsLow = 1;
muRelative = 1;

temporalPeriod = 1;
dutyCycle = 0.5;       % 仅用于方波
sinePhase = 0;         % 仅用于正弦调制，单位为弧度

%% ===================== 用户参数：k 扫描 =====================

% 纯时间晶体保持 k 守恒，因此每个 k 可以独立推进。k 不属于空间
% Brillouin 区，这里只是人为选择需要显示的 k/Omega 扫描范围。
kMinNormalized = -1.2;
kMaxNormalized = 1.2;
nK = 301;

%% ===================== 用户参数：FDTD 与时间原胞采样 =====================

% dx 只定义 Yee 空间差分的离散波数，不代表实际空间原胞。应通过减小 dx
% 检查数值色散收敛；dt 由每个时间周期的步数确定。
dx = 1/40;
stepsPerPeriod = 200;

% 这里只使用有限个时间原胞采样点，并不直接求无限时间本征值。在每个
% 原胞的同一调制相位记录一次 E；真实归一化频率分辨率约为
% 1/temporalCellCount，有限 Hann 窗还会使谱峰展宽。
temporalCellCount = 64;

% 零填充只细化纵轴显示网格，不能提升由 temporalCellCount 决定的分辨率。
zeroPaddingFactor = 4;
displayDynamicRangeDb = 60;

% 默认不返回两组初态的完整同相位复场样本；核心仍会返回未归一化列功率、
% 每周期主导增长率和离散元数据。科研排查时可临时改为 true。
returnTemporalSamples = false;

%% ===================== 参数检查与一个周期的材料定义 =====================

materialParameters = {epsHigh,epsLow,muRelative,temporalPeriod,dx};
validMaterialParameters = cellfun(@(value) isnumeric(value) && ...
    isscalar(value) && isreal(value) && isfinite(value) && value > 0, ...
    materialParameters);
if ~all(validMaterialParameters) || epsHigh < epsLow
    error('必须满足 epsHigh >= epsLow > 0，且 mu、T、dx 均为有限正实数。');
end
integerParameters = {nK,stepsPerPeriod,temporalCellCount,zeroPaddingFactor};
validIntegerParameters = cellfun(@(value) isnumeric(value) && ...
    isscalar(value) && isreal(value) && isfinite(value) && ...
    value == round(value),integerParameters);
if ~all(validIntegerParameters) || nK < 2 || ...
        stepsPerPeriod < 4 || temporalCellCount < 4 || zeroPaddingFactor < 1
    error('nK、stepsPerPeriod、temporalCellCount 和零填充因子必须是允许范围内的整数。');
end
if ~isnumeric(kMinNormalized) || ~isscalar(kMinNormalized) || ...
        ~isreal(kMinNormalized) || ~isfinite(kMinNormalized) || ...
        ~isnumeric(kMaxNormalized) || ~isscalar(kMaxNormalized) || ...
        ~isreal(kMaxNormalized) || ~isfinite(kMaxNormalized) || ...
        kMinNormalized >= kMaxNormalized
    error('kMinNormalized 和 kMaxNormalized 必须是有限实标量，且前者严格小于后者。');
end
if ~isnumeric(displayDynamicRangeDb) || ~isscalar(displayDynamicRangeDb) || ...
        ~isreal(displayDynamicRangeDb) || ~isfinite(displayDynamicRangeDb) || ...
        displayDynamicRangeDb <= 0
    error('displayDynamicRangeDb 必须是有限正实标量。');
end

T = temporalPeriod;
Omega = 2*pi/T;
dt = T/stepsPerPeriod;

modulationType = validate_modulation_type(modulationType);
switch lower(modulationType)
    case 'square'
        if ~isnumeric(dutyCycle) || ~isscalar(dutyCycle) || ...
                ~isreal(dutyCycle) || ~isfinite(dutyCycle)
            error('dutyCycle 必须是有限实标量。');
        end
        highStepsFloat = dutyCycle*stepsPerPeriod;
        if dutyCycle <= 0 || dutyCycle >= 1 || ...
                abs(highStepsFloat-round(highStepsFloat)) > 1e-10
            error('方波 dutyCycle*stepsPerPeriod 必须为整数，且 dutyCycle 位于 (0,1)。');
        end
        highSteps = round(highStepsFloat);
        lowSteps = stepsPerPeriod-highSteps;
        if highSteps < 2 || lowSteps < 2
            error('方波高、低介电状态均至少需要两个 FDTD 时间步。');
        end

        % 让 t=0 位于高介电层内部。两个 interfaceSteps 是一个周期内发生
        % 高->低和低->高突变的整数时间节点，供 D/B 中心修正使用。
        phaseSteps = floor(highSteps/2);
        phaseShift = phaseSteps*dt;
        highDuration = highSteps*dt;
        epsilonTime = @(time) epsLow+(epsHigh-epsLow).* ...
            double(mod(time+phaseShift,T) < highDuration);
        interfaceSteps = [highSteps-phaseSteps,stepsPerPeriod-phaseSteps];

    case 'sinusoidal'
        if ~isnumeric(sinePhase) || ~isscalar(sinePhase) || ...
                ~isreal(sinePhase) || ~isfinite(sinePhase)
            error('sinePhase 必须是有限实标量。');
        end
        epsMean = (epsHigh+epsLow)/2;
        epsAmplitude = (epsHigh-epsLow)/2;
        epsilonTime = @(time) epsMean+epsAmplitude*cos(Omega*time+sinePhase);
        interfaceSteps = [];

    otherwise
        error('modulationType 必须是 ''square'' 或 ''sinusoidal''。');
end

%% ===================== 逐 k 推进并重建 Floquet 准频率 =====================

kNormalized = linspace(kMinNormalized,kMaxNormalized,nK);
kScan = kNormalized*Omega;

bandCfg = struct();
bandCfg.temporalPeriod = T;
bandCfg.stepsPerPeriod = stepsPerPeriod;
bandCfg.temporalCellCount = temporalCellCount;
bandCfg.spatialStep = dx;
bandCfg.temporalInterfaceSteps = interfaceSteps;
bandCfg.zeroPaddingFactor = zeroPaddingFactor;
bandCfg.dynamicRangeDb = displayDynamicRangeDb;
bandCfg.returnTemporalSamples = returnTemporalSamples;

% fdtd_fft_bands 把每个 k 当作独立的 Yee Fourier 模式，同时用两组线性
% 独立 D/B 初态降低漏支风险。每个 k 都只保留有限个同相位时间采样点。
fftBands = fdtd_fft_bands(kScan,epsilonTime,muRelative,bandCfg);

%% ===================== 只绘制第一时间 Floquet 区能带 =====================

figure('Color','w','Position',[100 80 900 650]);
imagesc(fftBands.kNormalized,fftBands.omegaNormalized,fftBands.spectrumDb);
axis xy;
colormap(parula(256));
bandAxes = gca;
bandAxes.CLim = [-displayDynamicRangeDb 0];
bandColorbar = colorbar;
bandColorbar.Label.String = '每个 k 独立归一化的 FFT 功率 (dB)';
xlabel('k/\Omega');
ylabel('Re(\omega)/\Omega');
ylim([-0.5 0.5]);
title(sprintf('有限时间原胞采样重建的 Floquet 能带（N_T=%d）', ...
    temporalCellCount));
box on;

fprintf(['FDTD–FFT：Courant=%.6g，原生 Delta(omega/Omega)=%.6g，' ...
    '零填充显示点=%d；增长率仅作为返回诊断，不由谱宽读取。\n'], ...
    fftBands.courant,fftBands.nativeOmegaResolutionNormalized,fftBands.nFft);

% -------------------------------------------------------------------------
function modulationType = validate_modulation_type(modulationType)
% 接受字符行向量或标量 string，再统一转换成 switch 使用的字符文本。
if isstring(modulationType) && isscalar(modulationType)
    modulationType = char(modulationType);
elseif ~ischar(modulationType) || size(modulationType,1) ~= 1
    error('modulationType 必须是 ''square'' 或 ''sinusoidal'' 文本。');
end
end
