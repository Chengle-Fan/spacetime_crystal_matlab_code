%RUN_TL_FDTD_FFT Independent finite-chain x-t FFT band reconstruction.
% This script launches both a modulated run and a same-source, same-port
% unmodulated reference.  It does not read run_tl_fdtd_field workspace data.

clear; clc; close all;

physicalCfg = struct();
physicalCfg.topology = 'sspp';
physicalCfg.allowAssumptions = true;
physicalCfg.cellCount = 96;
model = tl_build_model(physicalCfg);

period = model.modulation.period;
stepsPerPeriod = max(128,ceil(period/(0.45*model.derived.maximumLeapfrogDt)));
dt = period/stepsPerPeriod;
recordPeriods = 16;
pulseCenter = 0.65*period;
pulseWidth = 0.12*period;
sourceVoltage = @(time) ((time-pulseCenter)/pulseWidth).* ...
    exp(-0.5*((time-pulseCenter)/pulseWidth).^2);

bandCfg = struct();
bandCfg.dt = dt;
bandCfg.nSteps = recordPeriods*stepsPerPeriod;
bandCfg.recordEvery = 1;
bandCfg.boundaryType = 'matched';
bandCfg.modulationEnabled = true;
bandCfg.precision = 'single';
bandCfg.source = struct('type','thevenin','node',1, ...
    'impedance',model.ports.Zsource,'waveformFcn',sourceVoltage);
bandField = tl_fdtd1d(model,bandCfg);

referenceCfg = bandCfg;
referenceCfg.modulationEnabled = false;
referenceField = tl_fdtd1d(model,referenceCfg);

roi = 9:(bandField.grid.branchCount-8);
timeRows = 1:(numel(bandField.node.t)-1); % remove repeated endpoint phase
cellIndex = bandField.branch.startNode(roi)-1;
channelOffset = 0.5*model.cell.a*ones(size(roi));
observation = struct( ...
    'data',bandField.branch.IAtNodeTime(timeRows,roi), ...
    'x',bandField.branch.x(roi), ...
    't',bandField.node.t(timeRows), ...
    'cellIndex',cellIndex, ...
    'channelId',ones(size(roi)), ...
    'channelOffset',channelOffset, ...
    'observableName','centered signed series current');
referenceObservation = observation;
referenceObservation.data = ...
    referenceField.branch.IAtNodeTime(timeRows,roi);

fftCfg = struct();
fftCfg.cellPeriod = model.cell.a;
fftCfg.temporalPeriod = period;
fftCfg.zeroPaddingTime = 2;
fftCfg.zeroPaddingSpace = 2;
fftCfg.dynamicRangeDb = 60;
fftCfg.referenceThresholdDb = 45;
fftCfg.ridgeThresholdDb = 25;
fftCfg.ridgeCount = 2;
fftBands = tl_xt_fft_bands(observation,referenceObservation,fftCfg);

figure('Color','w','Position',[100 100 900 650]);
imagesc(fftBands.kaOverPi,fftBands.omegaOverOmega,fftBands.spectrumDb);
axis xy;
ylim([-0.5 0.5]);
colormap(parula(256));
colorbar;
clim([-fftCfg.dynamicRangeDb 0]);
xlabel('ka/\pi');
ylabel('Re(\omega)/\Omega');
title('Finite-chain signed x-t FFT (source-supported column normalization)');

fprintf(['x-t FFT complete: %d periods, native Delta(omega/Omega)=%.4g, ' ...
    'native Delta(ka/pi)=%.4g.  Ridge linewidth is not Im(omega).\n'], ...
    fftBands.periodCount,fftBands.nativeOmegaResolution/ ...
    model.modulation.OmegaRadPerSec, ...
    fftBands.nativeKResolution*model.cell.a/pi);
