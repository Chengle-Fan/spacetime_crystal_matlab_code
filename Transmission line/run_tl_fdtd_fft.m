%RUN_TL_FDTD_FFT Finite-chain Gaussian-packet FDTD/FFT band response.
% Every horizontal column is an independent finite-chain simulation whose
% complex voltage packet has centre k_c and a stated intensity FWHM.  Fixed
% node-voltage probes are transformed independently and their powers are
% added before explicit Floquet folding.  The result is a source-weighted
% finite-sample response, not an infinite-bulk eigenspectrum.

clear; clc; close all;

%% Physical finite transmission-line sample

physicalCfg = struct();
physicalCfg.topology = 'sspp';
physicalCfg.allowAssumptions = true;
physicalCfg.cellCount = 192;
model = tl_build_model(physicalCfg);

T = model.modulation.period;
stepsPerPeriod = max(128,ceil(T/(0.45*model.derived.maximumLeapfrogDt)));
dt = T/stepsPerPeriod;
simulationPeriodCount = 16;
recordEvery = 1;

%% Gaussian source centres, finite width, and fixed probes

kScan = linspace(0,pi/model.cell.a,81);
pulseIntensityFwhmCells = 16;
voltageAmplitude = 1;
centerNode = round((model.finite.cellCount+1)/2);
probeOffsets = [-6 -3 0 3 6];
probeNodeIndices = centerNode+probeOffsets;
representativeK = 0.65*pi/model.cell.a;

scanCfg = struct();
scanCfg.kScan = kScan;
scanCfg.dt = dt;
scanCfg.nSteps = simulationPeriodCount*stepsPerPeriod;
scanCfg.recordEvery = recordEvery;
scanCfg.pulseIntensityFwhmCells = pulseIntensityFwhmCells;
scanCfg.voltageAmplitude = voltageAmplitude;
scanCfg.centerNode = centerNode;
scanCfg.probeNodeIndices = probeNodeIndices;
scanCfg.targetInitialFrequencyHz = model.modulation.fmHz/2;
scanCfg.representativeK = representativeK;
scanCfg.precision = 'single';
scanCfg.boundaryType = 'open';
scanCfg.modulationEnabled = true;
scanCfg.modulationStart = 0;
scanCfg.modulationEnd = Inf;
scanCfg.initialTailTolerance = 1e-4;
scanCfg.requireNoBoundaryArrival = true;
scanCfg.zeroKPropagationDirection = 1;
scanCfg.progressEvery = 10;
scanResult = tl_fdtd_gaussian_k_scan(model,scanCfg);

%% Endpoint-free integer-period FFT and explicit Floquet folding

fftCfg = struct();
fftCfg.temporalPeriod = T;
fftCfg.analysisTimeRange = [0 simulationPeriodCount*T];
fftCfg.zeroPaddingFactor = 4;
fftCfg.dynamicRangeDb = 60;
fftCfg.activeColumnRelativeThreshold = 1e-12;
fftCfg.returnProbePower = false;
fftCfg.returnProbeSignals = false;
fftBands = tl_fdtd_fft_bands( ...
    scanResult.probeSignals,scanResult.time,kScan,fftCfg);
fftBands.kaOverPi = kScan*model.cell.a/pi;

%% Representative real-space field and reconstructed response

representative = scanResult.representativeField;
figure('Color','w','Position',[100 100 960 620]);
imagesc(representative.node.x/model.cell.a, ...
    representative.node.t/T,real(representative.node.V));
axis xy; colorbar;
xlabel('node position x/a');
ylabel('t/T');
title(sprintf(['Representative finite-chain V(x,t), ' ...
    'k_c a/\\pi=%.3f'],representative.k*model.cell.a/pi));

figure('Color','w','Position',[120 100 900 650]);
imagesc(fftBands.kaOverPi,fftBands.omegaOverOmega,fftBands.spectrumDb);
axis xy;
ylim([-0.5 0.5]);
colormap(parula(256));
colorbar;
caxis([-fftCfg.dynamicRangeDb 0]); %#ok<CAXIS> R2020a compatibility
xlabel('Gaussian source centre k_c a/\pi');
ylabel('Re(\omega)/\Omega');
title('Finite-chain multi-probe V(t) FFT (per-k display normalization)');

fprintf(['FDTD-FFT complete: %d k centres, %d probes, intensity FWHM ' ...
    '%.3g cells, %d periods. Native Delta(omega/Omega)=%.4g; zero ' ...
    'padding changes only the plotted spacing.\n'],numel(kScan), ...
    fftBands.probeCount,pulseIntensityFwhmCells, ...
    fftBands.analysisPeriodCount, ...
    fftBands.nativeOmegaResolutionNormalized);
fprintf(['Boundary audit: tail amplitude %.3g, distance %.3g cells, ' ...
    'travel scale %.3g cells. Spectrum linewidth is not Im(omega).\n'], ...
    scanResult.initialEdgeAmplitudeBound, ...
    scanResult.distanceToBoundary/model.cell.a, ...
    scanResult.travelDistanceScale/model.cell.a);
