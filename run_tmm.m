%RUN_TMM  TMM bands of a square-wave photonic time crystal.
%
% Edit only the parameter block below. The script computes and plots the real
% and imaginary parts of the two quasifrequency bands. Normalized units use
% c0=eps0=mu0=1.

%% ===================== USER PARAMETERS =====================

% Square-wave material values.
epsHigh = 4;
epsLow  = 1;
muHigh  = 1;
muLow   = 1;

% Temporal structure.
T = 1;
dutyCycle = 0.5;       % fraction of one period at epsHigh/muHigh

% Wavenumber scan, normalized as k/Omega.
kMinNormalized = -1.65;
kMaxNormalized =  1.65;
nK = 301;

%% ===================== TMM CALCULATION =====================

if dutyCycle <= 0 || dutyCycle >= 1
    error('dutyCycle must be strictly between 0 and 1.');
end

Omega = 2*pi/T;
epsLayers = [epsHigh epsLow];
muLayers = [muHigh muLow];
durations = T*[dutyCycle 1-dutyCycle];
kNormalized = linspace(kMinNormalized,kMaxNormalized,nK);
kScan = kNormalized*Omega;
bands = tmm_bands(kScan,epsLayers,muLayers,durations);

%% ===================== REAL AND IMAGINARY BANDS =====================

kPlot = repmat(kNormalized,2,1);

figure('Color','w','Position',[100 100 1050 440]);
layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

axReal = nexttile(layout);
scatter(axReal,kPlot(:),real(bands.omega(:))/Omega,12,'filled');
xlabel(axReal,'k/\Omega');
ylabel(axReal,'Re(\omega)/\Omega');
title(axReal,'Real bands');
grid(axReal,'on');
box(axReal,'on');

axImag = nexttile(layout);
scatter(axImag,kPlot(:),imag(bands.omega(:))/Omega,12,'filled');
xlabel(axImag,'k/\Omega');
ylabel(axImag,'Im(\omega)/\Omega');
title(axImag,'Imaginary bands');
grid(axImag,'on');
box(axImag,'on');

title(layout,'Square-wave time crystal');
