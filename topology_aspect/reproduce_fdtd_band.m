function result = reproduce_fdtd_band()
%REPRODUCE_FDTD_BAND  用逐 k、有限时间原胞同相位采样重建论文参数能带。
% 只绘制 FDTD–FFT 能带映射；TMM 仅在后台量化初步峰位误差，不叠加曲线。

close all;

% 本文件只设置论文扫描参数和验证指标；逐 k 推进及 TMM 对照均强制来自
% 根目录模板，避免 topology_aspect 中的内核副本与主代码发生漂移。
[templateSource,templatePathCleanup] = ...
    use_root_templates({'fdtd_fft_bands','tmm_bands'});
if isempty(templatePathCleanup)
    error('根目录模板路径清理对象创建失败。');
end

%% ===================== 论文参数与有限采样参数 =====================

c0 = 0.299792458;          % um/fs
epsHigh = 3;
epsLow = 1;
muRelative = 1;
T = 2;                     % fs
OmegaModulation = 2*pi/T;
k0 = 2*pi/(T*c0);          % 论文图 1 横轴的参考波数，单位 1/um

% 与论文图 1 相同，扫描 0<=k/k0<=5。fdtd_fft_bands 使用光行时空间
% 坐标 x=z/c0，因此内部波数 kTemporal=c0*kPhysical=(k/k0)*2*pi/T。
kNormalized = linspace(0,5,501);
kPhysical = kNormalized*k0;
kTemporal = c0*kPhysical;

% 每周期 200 个时间步、空间步长 0.02 fs，真空 Courant 数为 0.5。
% 64 个同相位样本只构成有限时间窗；原生 Delta(omega/Omega)=1/64。
dxLightTime = 0.02;        % fs，对应物理空间步 c0*dx≈0.006 um
stepsPerPeriod = 200;
temporalCellCount = 64;
zeroPaddingFactor = 4;
displayDynamicRangeDb = 60;
dt = T/stepsPerPeriod;

%% ===================== 一个方波时间周期 =====================

highSteps = stepsPerPeriod/2;
phaseSteps = highSteps/2;
phaseShift = phaseSteps*dt;
highDuration = T/2;

% 令 t=0 位于 epsilon=3 时间层中心，与论文推导采用的时间反演对称原点
% 一致。两个界面严格位于 0.5 fs 和 1.5 fs 的整数 FDTD 节点。
epsilonTime = @(time) epsLow+(epsHigh-epsLow).* ...
    double(mod(time+phaseShift,T) < highDuration);
interfaceSteps = [highSteps-phaseSteps,stepsPerPeriod-phaseSteps];

bandCfg = struct();
bandCfg.temporalPeriod = T;
bandCfg.stepsPerPeriod = stepsPerPeriod;
bandCfg.temporalCellCount = temporalCellCount;
bandCfg.spatialStep = dxLightTime;
bandCfg.temporalInterfaceSteps = interfaceSteps;
bandCfg.zeroPaddingFactor = zeroPaddingFactor;
bandCfg.dynamicRangeDb = displayDynamicRangeDb;

fftBands = fdtd_fft_bands(kTemporal,epsilonTime,muRelative,bandCfg);

%% ===================== 只绘制有限采样 FDTD–FFT 能带 =====================

figure('Color','w','Position',[100 80 880 640]);
imagesc(fftBands.kNormalized,2*pi*fftBands.omegaNormalized, ...
    fftBands.spectrumDb);
axis xy;
bandAxes = gca;
bandAxes.CLim = [-displayDynamicRangeDb 0];
colormap(parula(256));
bandColorbar = colorbar;
bandColorbar.Label.String = '每个 k 独立归一化的有限窗 FFT 功率 (dB)';
xlabel('k/k_0，k_0=2\pi/(Tc_0)');
ylabel('\omega_F T');
xlim([0 5]);
ylim([-pi pi]);
title(sprintf('论文参数的有限时间原胞 FDTD–FFT 能带（N_T=%d）', ...
    temporalCellCount));
box on;

%% ===================== 后台初步峰位验证 =====================

% 精确 TMM 只用于报告误差，不进入图中。比较采用模 1 的循环准频率距离，
% 因为 +Omega/2 与 -Omega/2 是同一个 Floquet 区边界。
tmm = tmm_bands(kTemporal,[epsHigh epsLow],muRelative,[T/2 T/2]);
tmmNormalized = tmm.omega/OmegaModulation;
stableMask = max(abs(imag(tmmNormalized)),[],1) < 1e-9;
[maximumError,medianError,sampleCount] = compare_stable_ridges( ...
    fftBands.omegaNormalized,fftBands.spectrumDb,tmmNormalized,stableMask);

% 两组完备初态的后半窗对数振幅斜率给出正增长支估计。这里用与 FDTD
% 完全相同的 Yee 离散波数做 TMM 对照，把“有限窗拟合误差”和连续空间
% 数值色散分开；只比较增长大于 0.05/周期的带隙内部点。
kYee = (2/dxLightTime)*sin(kTemporal*dxLightTime/2);
tmmYee = tmm_bands(kYee,[epsHigh epsLow],muRelative,[T/2 T/2]);
theoryGrowthPerPeriod = max(imag(tmmYee.omega*T),[],1);
growthComparisonMask = theoryGrowthPerPeriod > 0.05;
growthErrors = abs(fftBands.dominantGrowthPerPeriod(growthComparisonMask)- ...
    theoryGrowthPerPeriod(growthComparisonMask));
medianGrowthError = median(growthErrors);
maximumGrowthError = max(growthErrors);

fprintf(['FDTD–FFT 对 TMM 通带峰位：样本数 %d，median |Delta omega|/Omega' ...
    ' = %.6g，max = %.6g。\n'],sampleCount,medianError,maximumError);
fprintf(['FDTD 后半窗增长拟合对 Yee-TMM k-gap：样本数 %d，median/max ' ...
    '|Delta Im(omega*T)| = %.6g / %.6g。\n'], ...
    nnz(growthComparisonMask),medianGrowthError,maximumGrowthError);

%% ===================== 返回复现数据 =====================

result.kNormalized = kNormalized;
result.kPhysicalPerUm = kPhysical;
result.omegaT = 2*pi*fftBands.omegaNormalized;
result.spectrumDb = fftBands.spectrumDb;
result.rawColumnPower = fftBands.rawColumnPower;
result.dominantGrowthPerPeriod = fftBands.dominantGrowthPerPeriod;
result.validation = struct('stableSampleCount',sampleCount, ...
    'medianErrorOverOmega',medianError,'maximumErrorOverOmega',maximumError, ...
    'growthSampleCount',nnz(growthComparisonMask), ...
    'medianGrowthErrorPerPeriod',medianGrowthError, ...
    'maximumGrowthErrorPerPeriod',maximumGrowthError);
result.paperParameters = struct('c0',c0,'epsHigh',epsHigh,'epsLow',epsLow, ...
    'muRelative',muRelative,'temporalPeriodFs',T,'k0PerUm',k0);
result.templateSource = templateSource;
result.samplingParameters = struct('dxLightTimeFs',dxLightTime, ...
    'dtFs',dt,'temporalCellCount',temporalCellCount, ...
    'zeroPaddingFactor',zeroPaddingFactor, ...
    'nFft',fftBands.nFft, ...
    'courantNumber',fftBands.courant, ...
    'nativeOmegaResolutionNormalized', ...
    fftBands.nativeOmegaResolutionNormalized, ...
    'windowName',fftBands.windowName);
end

% -------------------------------------------------------------------------
function [maximumError,medianError,sampleCount] = compare_stable_ridges( ...
        omegaAxis,spectrumDb,theoryOmega,stableMask)
% 在每条 TMM 通带附近独立寻找最强有限窗峰，并计算循环频率距离。
errors = [];
searchHalfWidth = 0.04;

for kIndex = find(stableMask)
    for branch = 1:size(theoryOmega,1)
        target = real(theoryOmega(branch,kIndex));
        circularOffset = mod(omegaAxis-target+0.5,1)-0.5;
        candidateRows = find(abs(circularOffset) <= searchHalfWidth);
        if isempty(candidateRows)
            continue;
        end
        [~,localMaximum] = max(spectrumDb(candidateRows,kIndex));
        peakRow = candidateRows(localMaximum);
        errors(end+1) = abs(circularOffset(peakRow)); %#ok<AGROW>
    end
end

if isempty(errors)
    error('没有找到可用于 FDTD–FFT/TMM 初步验证的通带峰。');
end
maximumError = max(errors);
medianError = median(errors);
sampleCount = numel(errors);
end
