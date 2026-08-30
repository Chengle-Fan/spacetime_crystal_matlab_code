%RUN_PWE  PWE bands of a one-dimensional photonic time crystal.
%
% Edit only the parameter block below. The script computes and plots the real
% and imaginary parts of the two quasifrequency bands. Normalized units use
% c0=eps0=mu0=1 and mu_r=1.

%% ===================== USER PARAMETERS =====================

% Material modulation: choose 'square' or 'sinusoidal'.
modulationType = 'square';

% Permittivity range shared by both modulation types.
epsHigh = 5;
epsLow  = 1;

% Square-wave parameters.
dutyCycle = 0.5;       % fraction of one period at epsHigh
timeShift = 0;         % temporal shift in the same units as T

% Sinusoidal parameter.
sinePhase = 0;         % radians; eps=mean+amplitude*cos(Omega*t+sinePhase)

% Period, Fourier truncation and material sampling.
T = 1;
Mtime = 19;            % retained orders; increase to check band convergence
Nt = 4096;             % increase together with Mtime for a square wave

% Wavenumber scan, normalized as k/Omega.
kMinNormalized = -2;
kMaxNormalized =  2;
nK = 301;

%% ===================== MATERIAL DEFINITION =====================

Omega = 2*pi/T;
switch lower(modulationType)
    case 'square'
        if dutyCycle <= 0 || dutyCycle >= 1
            error('dutyCycle must be strictly between 0 and 1.');
        end
        epsFun = @(t) epsLow + (epsHigh-epsLow) .* ...
            (mod(t-timeShift,T) < dutyCycle*T);
        profileName = 'Square-wave time crystal';

    case 'sinusoidal'
        epsMean = (epsHigh+epsLow)/2;
        epsAmplitude = (epsHigh-epsLow)/2;
        epsFun = @(t) epsMean + epsAmplitude*cos(Omega*t+sinePhase);
        profileName = 'Sinusoidal time crystal';

    otherwise
        error('modulationType must be ''square'' or ''sinusoidal''.');
end

if epsLow <= 0 || epsHigh < epsLow
    error('Require epsHigh >= epsLow > 0.');
end

%% ===================== PWE CALCULATION =====================

pweCfg = struct('T',T,'Mtime',Mtime,'Nt',Nt);
fourier = pwe_fourier(epsFun,pweCfg);
kNormalized = linspace(kMinNormalized,kMaxNormalized,nK);
kScan = kNormalized*Omega;
bands = pwe_bands(fourier,kScan);

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

title(layout,profileName);
