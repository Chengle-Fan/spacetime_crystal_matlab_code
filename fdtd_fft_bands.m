function result = fdtd_fft_bands(D, x, t, fftCfg)
%FDTD_FFT_BANDS  2-D FFT band structure of a D(x,t) field history.
%
% Converts a broadband FDTD field history into the Floquet band structure:
% spatial FFT along x, per-k Hann-windowed temporal FFT along t, power
% normalization, and Floquet folding of the spectrum into the first Brillouin
% zone. The FFT gives the REAL spectrum ridges only: it cannot recover the
% imaginary part of the quasifrequency inside a momentum/energy gap. For that
% one needs TMM (tmm_bands) or PWE (pwe_bands) — this function never calls
% either; theory curves may only be overlaid by the entry script.
%
% Folding:
%   * frequency is always folded modulo Omega = 2*pi/T into [-Omega/2, Omega/2];
%   * wavenumber is folded modulo g = 2*pi/Lambda into [-g/2, g/2] only when
%     fftCfg.spatialPeriod is supplied (spatially uniform media may leave it
%     empty to disable k folding).
%
% The FFT must be computed on grids spanning an INTEGER number of periods and
% WITHOUT a repeated endpoint. This function accepts either an endpoint-free
% periodic grid or a grid carrying one repeated endpoint (the natural FDTD
% record t = 0 .. nSteps*dt); in the latter case the repeated last sample is
% dropped automatically. Zero-padding refines the plot grid only — it never
% improves the true frequency resolution.
%
% INPUTS
%   D      : nTime x Nx field history (rows = time records, cols = space).
%   x      : length-Nx uniform periodic spatial grid (endpoint-free or with
%            one repeated endpoint).
%   t      : length-nTime uniform time grid.
%   fftCfg : struct with required field
%              temporalPeriod T
%            and optional fields
%              spatialPeriod     Lambda (enable k folding; default: none)
%              zeroPaddingFactor (default 4; plotting grid refinement only).
%
% OUTPUT (struct)
%   result.spectralDb            : nFold x nK folded spectral power in dB,
%                                  normalized to 0 dB at each k column.
%   result.omegaT                : nFold x 1 frequency grid in units of omega*T.
%   result.kGrid                 : nK x 1 folded (or raw) wavenumber grid.
%   result.kNormalized           : nK x 1 wavenumber in units of k*T/(2*pi).
%   result.power                 : nFold x nK folded power (linear, normalized).
%   result.ridgePositive         : nK x 1 positive-omega ridge (omega*T units).
%   result.ridgeNegative         : nK x 1 negative-omega ridge (omega*T units).
%   result.nativeFreqResolutionOmegaT     : 2*pi/nativePeriodCount.
%   result.zeroPaddedGridSpacingOmegaT    : 2*pi/nFold (display only).
%   result.spatialBinSpacingKOverK0       : T/(Nx*dx).
%   result.zeroPaddingFactor, result.temporalPeriod, result.spatialPeriod.

%==========================================================================
% Validate inputs
%==========================================================================
if nargin < 4 || isempty(fftCfg)
    error('fftCfg with field temporalPeriod is required.');
end
if ~isfield(fftCfg, 'temporalPeriod') || isempty(fftCfg.temporalPeriod)
    error('fftCfg.temporalPeriod (T) is required.');
end
T = fftCfg.temporalPeriod;
if ~isscalar(T) || ~isfinite(T) || T <= 0
    error('fftCfg.temporalPeriod must be a finite positive scalar.');
end
if isfield(fftCfg, 'spatialPeriod') && ~isempty(fftCfg.spatialPeriod)
    Lambda = fftCfg.spatialPeriod;
    if ~isscalar(Lambda) || ~isfinite(Lambda) || Lambda <= 0
        error('fftCfg.spatialPeriod must be a finite positive scalar.');
    end
    doFoldK = true;
else
    Lambda = [];
    doFoldK = false;
end
if isfield(fftCfg, 'zeroPaddingFactor') && ~isempty(fftCfg.zeroPaddingFactor)
    zeroPaddingFactor = fftCfg.zeroPaddingFactor;
else
    zeroPaddingFactor = 4;
end
if ~isscalar(zeroPaddingFactor) || zeroPaddingFactor < 1 || ...
        zeroPaddingFactor ~= round(zeroPaddingFactor)
    error('fftCfg.zeroPaddingFactor must be a positive integer.');
end

[nTime, Nx] = size(D);
if nTime < 2 || Nx < 2
    error('D must have at least 2 time rows and 2 spatial columns.');
end
if numel(t) ~= nTime || numel(x) ~= Nx
    error('D, x, and t must have consistent sizes.');
end
dtRecord = t(2) - t(1);
dx = x(2) - x(1);
if dtRecord <= 0 || any(abs(diff(t) - dtRecord) > 1e-10*max(1, abs(dtRecord)))
    error('t must be uniformly spaced with positive dt.');
end
if dx <= 0 || any(abs(diff(x) - dx) > 1e-10*max(1, abs(dx)))
    error('x must be uniformly spaced with positive dx.');
end

%==========================================================================
% Endpoint-free integer-period temporal grid (drop one repeated endpoint)
%==========================================================================
periodCountFull = nTime*dtRecord/T;
periodCountDrop = (nTime - 1)*dtRecord/T;
if abs(periodCountFull - round(periodCountFull)) < 1e-9
    nTimeEff = nTime;                       % already endpoint-free
    dropLast = false;
elseif abs(periodCountDrop - round(periodCountDrop)) < 1e-9
    nTimeEff = nTime - 1;                   % drop the repeated endpoint
    dropLast = true;
else
    error(['The FFT grid must span an integer number of modulation periods ', ...
        '(endpoint-free or with one repeated endpoint).']);
end
nativePeriodCount = round((nTimeEff)*dtRecord/T);
if dropLast
    D = D(1:nTimeEff, :);
end

%==========================================================================
% Spatial FFT, per-k normalization, temporal Hann window + FFT
%==========================================================================
Dk = fft(D, [], 2);
kFull = (0:Nx-1)'*(2*pi/(Nx*dx));

% Per-k column normalization: prevents unstable gap columns from drowning the
% pass bands in the display.
columnNorm = sqrt(sum(abs(Dk).^2, 1));
columnNorm(columnNorm == 0) = 1;
Dk = Dk./columnNorm;

window = 0.5 - 0.5*cos(2*pi*(0:nTimeEff-1)'/nTimeEff);
Dk = Dk.*window;

nFft = zeroPaddingFactor*nTimeEff;
spectrum = fftshift(ifft(Dk, nFft, 1), 1);
power = abs(spectrum).^2;                   % nFft x Nx
clear spectrum Dk;

%==========================================================================
% Floquet folding in frequency (always)
%==========================================================================
nFoldFloat = nFft*dtRecord/T;
nFold = round(nFoldFloat);
if abs(nFold - nFoldFloat) > 1e-10 || mod(nFold, 2) ~= 0 || ...
        mod(nFft, nFold) ~= 0
    error('FFT sampling is not commensurate with the modulation period.');
end
nHarmonicBlocks = nFft/nFold;
powerFold = sum(reshape(power, nFold, nHarmonicBlocks, Nx), 2);
powerFold = reshape(powerFold, nFold, Nx);
powerFold = fftshift(powerFold, 1);
omegaT = (-nFold/2:nFold/2-1)'*(2*pi/nFold);

%==========================================================================
% Floquet folding in wavenumber (optional)
%==========================================================================
if doFoldK
    nCells = Nx*dx/Lambda;
    nCellsFull = Nx*dx/Lambda;
    nCellsDrop = (Nx-1)*dx/Lambda;
    if abs(nCellsFull - round(nCellsFull)) < 1e-9
        nCellsEff = round(nCellsFull);
    elseif abs(nCellsDrop - round(nCellsDrop)) < 1e-9
        nCellsEff = round(nCellsDrop);
        powerFold = powerFold(:, 1:Nx-1);
        Nx = Nx - 1;
    else
        error(['The spatial grid must span an integer number of spatial ', ...
            'periods to fold k (endpoint-free or with one repeated endpoint).']);
    end
    nFoldK = nCellsEff;
    nK = nFoldK;
    kReplicas = Nx/nFoldK;
    if kReplicas ~= round(kReplicas)
        error('Spatial period count does not divide the spatial grid.');
    end
    powerFold = sum(reshape(powerFold, nFold, nFoldK, kReplicas), 3);
    powerFold = fftshift(powerFold, 2);
    g = 2*pi/Lambda;
    kGrid = (-nFoldK/2:nFoldK/2-1)'*(g/nFoldK);
else
    nK = Nx;
    kGrid = kFull;
end
kNormalized = kGrid*T/(2*pi);

%==========================================================================
% Per-k normalization, dB, and ridge extraction
%==========================================================================
columnMaximum = max(powerFold, [], 1);
columnMaximum(columnMaximum == 0) = 1;
powerFold = powerFold./columnMaximum;
spectralDb = 10*log10(max(powerFold, 1e-8));

positiveRows = find(omegaT >= 0);
negativeRows = find(omegaT <= 0);
[~, posPeak] = max(powerFold(positiveRows, :), [], 1);
[~, negPeak] = max(powerFold(negativeRows, :), [], 1);
ridgePositive = omegaT(positiveRows(posPeak)).';
ridgeNegative = omegaT(negativeRows(negPeak)).';

%==========================================================================
% Output
%==========================================================================
result.spectralDb    = spectralDb;
result.omegaT        = omegaT;
result.kGrid         = kGrid;
result.kNormalized   = kNormalized;
result.power         = powerFold;
result.ridgePositive = ridgePositive;
result.ridgeNegative = ridgeNegative;
result.nativeFreqResolutionOmegaT  = 2*pi/nativePeriodCount;
result.zeroPaddedGridSpacingOmegaT = 2*pi/nFold;
result.spatialBinSpacingKOverK0    = T/(Nx*dx);
result.zeroPaddingFactor = zeroPaddingFactor;
result.temporalPeriod = T;
result.spatialPeriod  = Lambda;
end
