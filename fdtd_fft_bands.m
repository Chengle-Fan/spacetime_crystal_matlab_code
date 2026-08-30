function result = fdtd_fft_bands(fieldData, x, t, fftCfg)
%FDTD_FFT_BANDS  2-D FFT band structure of an E(x,t) (or D) field history.
%
% Converts a broadband FDTD field history into the Floquet band structure:
% spatial FFT along x, per-k windowed temporal FFT along t, power
% normalization, and (optional) Floquet folding into the first Brillouin
% zone. The FFT gives the REAL spectrum ridges only: it cannot recover the
% imaginary part of the quasifrequency inside a momentum/energy gap. For that
% one needs TMM (tmm_bands) or PWE (pwe_bands) — this function never calls
% either; theory curves may only be overlaid by the entry script.
%
% Experimental observable: the acceptance observable is the ELECTRIC FIELD
% E(x,t) sampled at the integer Yee nodes (fieldComponent = 'E', the
% default). D may be passed only as an optional diagnostic. The recorded E
% value at an ideal temporal switch is the after-switch constitutive value
% D/eps(t+), because D is continuous and E is recovered at the integer node.
%
% Grids must span an INTEGER number of modulation periods WITHOUT a repeated
% endpoint. The time grid may carry one repeated endpoint (the natural FDTD
% record t = 0 .. nSteps*dt); it is dropped first. A repeated SPATIAL
% endpoint, when present, is also detected and dropped BEFORE the spatial FFT
% (never after). In the no-k-fold mode the endpoint semantics cannot be
% inferred reliably: pass fftCfg.hasRepeatedSpatialEndpoint = true to drop
% one repeated spatial endpoint, otherwise the grid is assumed endpoint-free.
% Zero-padding refines the plot grid only — it never improves the true
% resolution; the effective resolution is set by the window ENBW, not by the
% DFT bin spacing.
%
% Measurement ROI contract (P1-10): a finite-sample band analysis must NOT
% FFT the whole recorded domain (sample edges, background zones and the
% source would all leak into every k column). Select the sample-interior
% measurement window with fftCfg.spatialROI = [iStart iEnd] (grid indices on
% the cleaned spatial grid) and the steady-state time window with
% fftCfg.timeROI = [iStart iEnd] (record indices on the cleaned time grid,
% i.e. after any repeated temporal endpoint is dropped). Windows are applied
% inside the ROI, before the FFTs. For a folded transform the ROI must be
% endpoint-free and span an integer number of periods: the spatial ROI an
% integer number of spatial periods (when k folding is enabled) and the time
% ROI an integer number of modulation periods (when temporal folding is
% enabled). The effective resolution is the window ENBW times the DFT bin
% spacing (2*pi/L_ROI, 2*pi/nPeriodsROI); zero-padding only refines the
% plot grid.
%
% INPUTS
%   fieldData : nTime x Nx field history (rows = time records, cols = space).
%   x         : length-Nx uniform grid (endpoint-free, or with one repeated
%               endpoint when that is explicitly marked/detected).
%   t         : length-nTime uniform time grid.
%   fftCfg    : struct with required field
%                 temporalPeriod T
%               and optional fields
%                 spatialPeriod          Lambda (enable k folding; default none)
%                 zeroPaddingFactor      (default 4; FFT length factor.
%                                        Fractional values e.g. 1/4 are
%                                        allowed; the FFT length nFft =
%                                        zeroPaddingFactor*nTimeEff must stay
%                                        a positive integer commensurate with
%                                        the modulation period)
%                 fieldComponent         ('E' default, 'D'/'H' diagnostic)
%                 foldTemporal           (logical, default true)
%                 temporalWindow         'hann' (default) | 'rect' | 'tukey'
%                 temporalWindowTukeyR   (default 0.5)
%                 spatialWindow          'rect' (default) | 'hann' | 'tukey'
%                 spatialWindowTukeyR    (default 0.5)
%                 noiseFloor             relative per-k active threshold
%                                        (default 1e-9, in [0,1))
%                 hasRepeatedSpatialEndpoint (logical, default false; only
%                                        meaningful when k folding is off)
%                 spatialROI             [iStart iEnd] spatial grid indices of
%                                        the measurement ROI (default: full
%                                        cleaned grid)
%                 timeROI                [iStart iEnd] time record indices of
%                                        the measurement window (default: full
%                                        cleaned record)
%                 measurementROI         optional 2-vector [x1 x2] physical
%                                        ROI, recorded in the output
%                 timeWindowStart        optional scalar t of window start
%                 modulationPhase        optional scalar phase, recorded
%
% OUTPUT (struct)
%   result.spectralDb            : nFold x nK folded spectral power in dB,
%                                  normalized to 0 dB per ACTIVE k column.
%   result.power                 : nFold x nK folded power (linear, per-k
%                                  normalized; inactive columns stay 0).
%   result.powerRaw              : nFold x nK folded raw physical power
%                                  (sum of |FFT|^2 per Floquet class).
%   result.rawSpectralPower      : nTimeEff x nK raw signed (k,omega)
%                                  response spectrum (native resolution, no
%                                  temporal folding, no zero-padding).
%   result.omegaT                : nFold x 1 frequency grid (omega*T units).
%   result.rawOmegaGridT         : nTimeEff x 1 signed native omega*T grid.
%   result.kGrid                 : nK x 1 signed folded (or raw) wavenumber.
%   result.kNormalized           : nK x 1 wavenumber in units of k*T/(2*pi).
%   result.ridgePositive         : nK x 1 positive-omega ridge (omega*T units),
%                                  NaN on inactive k columns.
%   result.ridgeNegative         : nK x 1 negative-omega ridge, NaN inactive.
%   result.activeKMask           : nK x 1 logical; true on reliable columns.
%   result.rawColumnPower        : nK x 1 per-k spatial-spectrum energy.
%   result.dftBinSpacingOmegaT   : 2*pi/nativePeriodCount (bin spacing only).
%   result.effectiveFreqResolutionOmegaT : ENBW * bin spacing (real limit).
%   result.zeroPaddedGridSpacingOmegaT   : 2*pi/nFold (display only).
%   result.spatialDftBinSpacingK       : 2*pi/L_ROI (bin spacing only).
%   result.spatialDftBinSpacingKOverK0 : T/(Nx*dx) (bin spacing only).
%   result.spatialEffectiveResolutionKOverK0 : ENBW_spatial * bin spacing.
%   result.spatialROI, result.timeROI  : index ranges actually used.
%   result.roiXStart/.roiXEnd/.roiXLength, result.roiTStart/.roiTEnd/
%       .roiTLength : physical ROI extents; result.roiPhaseAtStart =
%       mod(t(roiTStart), T)/T (modulation phase at the time-window start).
%   result.temporalWindow/.temporalWindowGains(.coherentGain,.energyGain,
%       .enbwBins), .spatialWindow/.spatialWindowGains.
%   result.timeEndpointDropped, result.spatialEndpointDropped (logicals).
%   result.fieldComponent, result.zeroPaddingFactor, result.temporalPeriod,
%   result.spatialPeriod, result.measurementROI, result.timeWindowStart,
%   result.modulationPhase, result.interfaceSideConvention.

%==========================================================================
% Validate configuration
%==========================================================================
if nargin < 4 || isempty(fftCfg)
    error('fftCfg with field temporalPeriod is required.');
end
if ~isfield(fftCfg, 'temporalPeriod') || isempty(fftCfg.temporalPeriod)
    error('fftCfg.temporalPeriod (T) is required.');
end
T = fftCfg.temporalPeriod;
if ~isscalar(T) || ~isreal(T) || ~isfinite(T) || T <= 0
    error('fftCfg.temporalPeriod must be a finite positive scalar.');
end

if isfield(fftCfg, 'spatialPeriod') && ~isempty(fftCfg.spatialPeriod)
    Lambda = fftCfg.spatialPeriod;
    if ~isscalar(Lambda) || ~isreal(Lambda) || ~isfinite(Lambda) || Lambda <= 0
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
if ~isscalar(zeroPaddingFactor) || ~isreal(zeroPaddingFactor) || ...
        ~isfinite(zeroPaddingFactor) || zeroPaddingFactor <= 0
    error(['fftCfg.zeroPaddingFactor must be a finite positive scalar. A ', ...
        'fractional factor (e.g. 1/4) is allowed; the FFT length ', ...
        'nFft = zeroPaddingFactor*nTimeEff must still be a positive ', ...
        'integer commensurate with the modulation period, which is ', ...
        'checked in the fold step.']);
end

if isfield(fftCfg, 'foldTemporal') && ~isempty(fftCfg.foldTemporal)
    foldTemporal = logical_scalar(fftCfg.foldTemporal, 'fftCfg.foldTemporal');
else
    foldTemporal = true;
end
if isfield(fftCfg, 'hasRepeatedSpatialEndpoint') && ...
        ~isempty(fftCfg.hasRepeatedSpatialEndpoint)
    hasRepeatedSpatialEndpoint = logical_scalar( ...
        fftCfg.hasRepeatedSpatialEndpoint, 'fftCfg.hasRepeatedSpatialEndpoint');
else
    hasRepeatedSpatialEndpoint = false;
end

if isfield(fftCfg, 'fieldComponent') && ~isempty(fftCfg.fieldComponent)
    fieldComponent = char(fftCfg.fieldComponent);
else
    fieldComponent = 'E';
end
if ~ismember(fieldComponent, {'E','D','H'})
    error('fftCfg.fieldComponent must be ''E'' (default), ''D'', or ''H''.');
end

if isfield(fftCfg, 'noiseFloor') && ~isempty(fftCfg.noiseFloor)
    noiseFloor = fftCfg.noiseFloor;
    if ~isscalar(noiseFloor) || ~isreal(noiseFloor) || ~isfinite(noiseFloor) || ...
            noiseFloor < 0 || noiseFloor >= 1
        error('fftCfg.noiseFloor must be a finite real scalar in [0, 1).');
    end
else
    noiseFloor = 1e-9;
end

[temporalWindow, temporalWindowTukeyR] = parse_window(fftCfg, 'temporalWindow', ...
    'temporalWindowTukeyR', 'hann');
[spatialWindow, spatialWindowTukeyR] = parse_window(fftCfg, 'spatialWindow', ...
    'spatialWindowTukeyR', 'rect');

if isfield(fftCfg, 'measurementROI') && ~isempty(fftCfg.measurementROI)
    measurementROI = fftCfg.measurementROI(:).';
    if numel(measurementROI) ~= 2 || ~all(isfinite(measurementROI)) || ...
            measurementROI(1) >= measurementROI(2)
        error('fftCfg.measurementROI must be a finite 2-vector [x1 x2] with x1 < x2.');
    end
else
    measurementROI = [];
end
if isfield(fftCfg, 'timeWindowStart') && ~isempty(fftCfg.timeWindowStart)
    timeWindowStart = fftCfg.timeWindowStart;
    if ~isscalar(timeWindowStart) || ~isfinite(timeWindowStart)
        error('fftCfg.timeWindowStart must be a finite scalar.');
    end
else
    timeWindowStart = [];
end
if isfield(fftCfg, 'modulationPhase') && ~isempty(fftCfg.modulationPhase)
    modulationPhase = fftCfg.modulationPhase;
    if ~isscalar(modulationPhase) || ~isfinite(modulationPhase)
        error('fftCfg.modulationPhase must be a finite scalar.');
    end
else
    modulationPhase = [];
end

% Measurement ROI as index ranges on the cleaned grids (default: full grid /
% full record). spatialROI = [iStart iEnd] spatial columns, timeROI = [iStart
% iEnd] time rows. Indices are interpreted AFTER any repeated endpoint is
% dropped, so the caller must account for the endpoint semantics.
if isfield(fftCfg, 'spatialROI') && ~isempty(fftCfg.spatialROI)
    spatialROI = fftCfg.spatialROI(:).';
    if numel(spatialROI) ~= 2 || any(~isfinite(spatialROI)) || ...
            any(spatialROI ~= round(spatialROI)) || spatialROI(1) < 1 || ...
            spatialROI(2) <= spatialROI(1)
        error(['fftCfg.spatialROI must be a 2-vector [iStart iEnd] of ', ...
            'integer grid indices with 1 <= iStart < iEnd (at least two ', ...
            'points are required for an FFT).']);
    end
else
    spatialROI = [];
end
if isfield(fftCfg, 'timeROI') && ~isempty(fftCfg.timeROI)
    timeROI = fftCfg.timeROI(:).';
    if numel(timeROI) ~= 2 || any(~isfinite(timeROI)) || ...
            any(timeROI ~= round(timeROI)) || timeROI(1) < 1 || ...
            timeROI(2) <= timeROI(1)
        error(['fftCfg.timeROI must be a 2-vector [iStart iEnd] of integer ', ...
            'record indices with 1 <= iStart < iEnd (at least two points ', ...
            'are required for an FFT).']);
    end
else
    timeROI = [];
end

%==========================================================================
% Validate field data and grids
%==========================================================================
if ~isnumeric(fieldData) || ~all(isfinite(fieldData), 'all')
    error('fieldData must be a finite numeric array.');
end
[nTime, Nx] = size(fieldData);
if nTime < 2 || Nx < 2
    error('fieldData must have at least 2 time rows and 2 spatial columns.');
end
if numel(t) ~= nTime || numel(x) ~= Nx
    error('fieldData, x, and t must have consistent sizes.');
end
if ~isreal(t) || ~all(isfinite(t)) || ~isreal(x) || ~all(isfinite(x))
    error('x and t must be finite real vectors.');
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
% Temporal grid: integer number of periods, endpoint-free (drop one repeated
% endpoint before any transform)
%==========================================================================
periodCountFull = nTime*dtRecord/T;
periodCountDrop = (nTime - 1)*dtRecord/T;
if abs(periodCountFull - round(periodCountFull)) < 1e-9
    nTimeEff = nTime;
    timeEndpointDropped = false;
elseif abs(periodCountDrop - round(periodCountDrop)) < 1e-9
    nTimeEff = nTime - 1;
    timeEndpointDropped = true;
    fieldData = fieldData(1:nTimeEff, :);
else
    error(['The FFT grid must span an integer number of modulation periods ', ...
        '(endpoint-free or with one repeated endpoint).']);
end

%==========================================================================
% Spatial grid: detect/drop ONE repeated spatial endpoint BEFORE the FFT
% (P1-06) so the DFT and the k grid are built on the true physical samples.
% The k-fold integer-cell gate for doFoldK is applied AFTER the spatial ROI
% crop, on the ROI columns themselves (an interior ROI may span an integer
% number of unit cells even when the full grid does not).
%==========================================================================
spatialEndpointDropped = false;
repeatedMismatch = [];
repeatedScale = [];
if ~doFoldK && hasRepeatedSpatialEndpoint
    spatialEndpointDropped = true;
    repeatedCol = fieldData(:, end);   % same physical point as column 1
    fieldData = fieldData(:, 1:Nx-1);
    x = x(1:Nx-1);
    Nx = Nx - 1;
    repeatedMismatch = max(abs(fieldData(:,1) - repeatedCol), [], 'all');
    repeatedScale = max(1, max(abs(fieldData), [], 'all'));
end

%==========================================================================
% Measurement ROI crops (P1-10). timeROI rows are on the cleaned time grid
% (after the repeated temporal endpoint is dropped), spatialROI columns on
% the cleaned spatial grid (after any repeated spatial endpoint is dropped).
%==========================================================================
if ~isempty(timeROI)
    if timeROI(2) > nTimeEff
        error(['fftCfg.timeROI indices must lie on the cleaned record ', ...
            '(1..%d); account for a dropped repeated temporal endpoint.'], ...
            nTimeEff);
    end
    fieldData = fieldData(timeROI(1):timeROI(2), :);
    tROI = t(timeROI(1):timeROI(2));
    nTimeEff = numel(tROI);
    dtRecord = tROI(2) - tROI(1);   % unchanged by an integer-index crop
else
    tROI = t(1:nTimeEff);
    timeROI = [1 nTimeEff];
end
if ~isempty(spatialROI)
    if spatialROI(2) > Nx
        error(['fftCfg.spatialROI indices must lie on the cleaned spatial ', ...
            'grid (1..%d); account for a dropped repeated spatial endpoint.'], ...
            Nx);
    end
    fieldData = fieldData(:, spatialROI(1):spatialROI(2));
    xROI = x(spatialROI(1):spatialROI(2));
    Nx = numel(xROI);
    dx = xROI(2) - xROI(1);         % unchanged by an integer-index crop
else
    xROI = x;
    spatialROI = [1 Nx];
end
% k-fold integer-cell gate on the FFT domain: the spatial ROI (or the full
% grid when no spatialROI is given) must span an integer number of unit
% cells to fold k (endpoint-free or with one repeated endpoint inside).
if doFoldK
    nCellsFull = Nx*dx/Lambda;
    nCellsDrop = (Nx-1)*dx/Lambda;
    if abs(nCellsFull - round(nCellsFull)) < 1e-9
        % endpoint-free ROI
    elseif abs(nCellsDrop - round(nCellsDrop)) < 1e-9
        spatialEndpointDropped = true;
        repeatedCol = fieldData(:, end);   % same physical point as ROI col 1
        fieldData = fieldData(:, 1:Nx-1);
        xROI = xROI(1:Nx-1);
        Nx = Nx - 1;
        spatialROI(2) = spatialROI(2) - 1; % FFT domain now excludes the copy
        repeatedMismatch = max(abs(fieldData(:,1) - repeatedCol), [], 'all');
        repeatedScale = max(1, max(abs(fieldData), [], 'all'));
    else
        error(['The spatial FFT domain (the spatial ROI, or the full grid ', ...
            'when no spatialROI is given) must span an integer number of ', ...
            'spatial periods to fold k; got %.6f periods over %d cells ', ...
            '(Lambda %g, dx %g).'], nCellsFull, Nx, Lambda, dx);
    end
end
if spatialEndpointDropped
    % Sanity: the repeated endpoint and column 1 are the SAME physical point
    % up to solver error (the dropped column is the +1-period copy).
    if repeatedMismatch > 1e-2*repeatedScale
        warning('fdtd_fft_bands:endpointMismatch', ...
            ['Dropping a repeated spatial endpoint whose first/last field ', ...
            'samples differ by %.3e (field scale %.3e); the dropped point ', ...
            'is not consistent with periodicity.'], repeatedMismatch, ...
            repeatedScale);
    end
end
% Number of modulation periods in the (possibly cropped) time ROI. Used for
% the native frequency grid and the DFT-bin-spacing reporting. Kept EXACT
% (may be non-integer): the native omega*T bin spacing is 2*pi/periodCount.
nativePeriodCount = nTimeEff*dtRecord/T;

%==========================================================================
% Spatial window (applied inside the ROI, before the spatial FFT) and gains
%==========================================================================
wSpatial = make_window(spatialWindow, spatialWindowTukeyR, Nx, 'row');
[spatialCG, spatialEG, spatialENBW] = window_gains(wSpatial);
fieldData = fieldData .* wSpatial;

%==========================================================================
% Spatial FFT (unshifted), per-k column energy, temporal window + FFT
%==========================================================================
Dk = fft(fieldData, [], 2);                       % unshifted spatial bins
rawColumnPower = sum(abs(Dk).^2, 1);              % 1 x Nx per-k energy
rawColumnPower = rawColumnPower(:);               % nK0 x 1 (pre-fold)

wTemporal = make_window(temporalWindow, temporalWindowTukeyR, nTimeEff, 'col');
[temporalCG, temporalEG, temporalENBW] = window_gains(wTemporal);
Dk = Dk .* wTemporal;

% Native (no zero-padding, no temporal folding) signed spectrum.
specNative = abs(ifft(Dk, nTimeEff, 1)).^2;       % nTimeEff x Nx
freqPermNative = centered_index(nTimeEff);
specNative = specNative(freqPermNative, :);
rawOmegaGridT = (-floor(nTimeEff/2):ceil(nTimeEff/2)-1)' ...
    *(2*pi/nativePeriodCount);

%==========================================================================
% Optional temporal Floquet folding (P1-07): accumulate harmonic copies on
% the UNSHIFTED bins modulo nFold, then apply a single centered reorder.
% This is parity-correct for odd/even period counts and samples per period.
%==========================================================================
if foldTemporal
    roiPeriodCount = nTimeEff*dtRecord/T;
    if abs(roiPeriodCount - round(roiPeriodCount)) > 1e-9
        error(['The time ROI must span an integer number of modulation ', ...
            'periods to fold omega (got %.6f). Select a timeROI that starts ', ...
            'and ends on period boundaries.'], roiPeriodCount);
    end
    nFft = zeroPaddingFactor*nTimeEff;
    if nFft < 1 || abs(nFft - round(nFft)) > 1e-10
        error(['fftCfg.zeroPaddingFactor*nTimeEff = %.6g must be a positive ', ...
            'integer (zeroPaddingFactor %.4g, %d time samples).'], ...
            nFft, zeroPaddingFactor, nTimeEff);
    end
    nFft = round(nFft);
    nFoldFloat = nFft*dtRecord/T;
    nFold = round(nFoldFloat);
    if nFold < 1 || abs(nFold - nFoldFloat) > 1e-10 || mod(nFft, nFold) ~= 0
        error(['FFT sampling is not commensurate with the modulation period ', ...
            '(nFold = %.6g is not a positive integer, or nFft = %d is not a ', ...
            'multiple of it).'], nFoldFloat, nFft);
    end
    specPad = abs(ifft(Dk, nFft, 1)).^2; % nFft x Nx
    nHarmonicBlocks = nFft/nFold;
    % Sum the Floquet harmonic copies directly on the UNSHIFTED bins: bin b
    % and b+nFold differ by exactly Omega, so mod(b, nFold) is the class.
    % A single centered reorder at the end fixes odd/even parity (P1-07).
    powerFold = reshape(specPad, nFold, nHarmonicBlocks, Nx);
    powerFold = sum(powerFold, 2);
    powerFold = reshape(powerFold, nFold, Nx);
    powerFold = powerFold(centered_index(nFold), :);
    omegaT = (-floor(nFold/2):ceil(nFold/2)-1)'*(2*pi/nFold);
    nFreqOut = nFold;
else
    powerFold = specNative;
    omegaT = rawOmegaGridT;
    nFreqOut = nTimeEff;
end

%==========================================================================
% Spatial Floquet folding (optional, P1-08): sum spatial replicas, then a
% centered k-grid — parity-correct for odd/even period counts.
%==========================================================================
if doFoldK
    nCellsActual = Nx*dx/Lambda;
    if abs(nCellsActual - round(nCellsActual)) > 1e-9
        error(['The spatial ROI must span an integer number of spatial ', ...
            'periods to fold k (got %.6f cells). Select a spatialROI with ', ...
            'integer unit-cell extent.'], nCellsActual);
    end
    nFoldK = round(nCellsActual);
    kReplicas = Nx/nFoldK;
    if abs(kReplicas - round(kReplicas)) > 1e-9
        error('Spatial period count does not divide the spatial ROI.');
    end
    specNative = sum(reshape(specNative, nTimeEff, nFoldK, kReplicas), 3);
    powerFold  = sum(reshape(powerFold, nFreqOut, nFoldK, kReplicas), 3);
    % Per-k column energy in spatial-fold mode must be the FOLDED
    % per-Floquet-class energy (sum over the kReplicas replicas), consistent
    % with powerFold, not the energy of a single replica (P1-09/P1-10).
    % rawColumnPower is in unshifted-bin order here; the class sum keeps the
    % same class indexing as the reshape/sum above.
    rawColumnPower = sum(reshape(rawColumnPower, nFoldK, kReplicas), 2);
    kPerm = centered_index(nFoldK);
    specNative = specNative(:, kPerm);
    powerFold  = powerFold(:, kPerm);
    nK = nFoldK;
    g = 2*pi/Lambda;
    kGrid = (-floor(nFoldK/2):ceil(nFoldK/2)-1)'*(g/nFoldK);
else
    % Signed Nyquist k grid (P3-05): reorder columns to ascending signed k.
    kPerm = centered_index(Nx);
    specNative = specNative(:, kPerm);
    powerFold  = powerFold(:, kPerm);
    nK = Nx;
    kGrid = (-floor(Nx/2):ceil(Nx/2)-1)'*(2*pi/(Nx*dx));
end
kNormalized = kGrid*T/(2*pi);
rawColumnPower = rawColumnPower(kPerm);

%==========================================================================
% Active-k mask (P1-09): only normalize columns with real excitation. The
% threshold combines a ROBUST reference column energy, machine precision and
% the configurable relative noiseFloor. Inactive columns keep zero power and
% NaN ridges instead of amplifying float noise into fake bands.
%
% The reference is the STRONGEST column energy, not the median. A column is
% active iff its per-k energy is within noiseFloor of the strongest column.
% For a sparse pure +/-k mode (P1-09 / repair.md 7.317) the median lands
% INSIDE the float-noise floor, so a median reference would mark every
% column active and amplify float noise into 0 dB pseudo-ridges; the max
% reference keeps the signal columns and rejects the noise columns. A
% parametric-growth (momentum-gap) column can only raise this relative
% floor: columns more than 1/noiseFloor below the dominant signal are not
% resolvable in the finite-time-window spectrum anyway, so the mask stays
% faithful in that regime as well.
%==========================================================================
colEnergy = rawColumnPower;
posEn = colEnergy(colEnergy > 0);
if isempty(posEn)
    activeKMask = false(nK, 1);
else
    refEnergy = max(posEn);
    activeThreshold = max(eps(refEnergy)*1e3, noiseFloor*refEnergy);
    activeKMask = colEnergy >= activeThreshold;
end

% Raw physical power (never per-k normalized) and per-k-normalized display.
powerRaw = powerFold;
columnMax = max(powerFold, [], 1);
columnMax(~activeKMask) = 1;                    % inactive stays 0
powerNorm = powerFold ./ columnMax;
powerNorm(:, ~activeKMask) = 0;

spectralDb = 10*log10(max(powerNorm, 1e-8));
spectralDb(:, ~activeKMask) = -inf;

%==========================================================================
% Ridge extraction on active columns only (P3-03: nK x 1 columns)
%==========================================================================
ridgePositive = nan(nK, 1);
ridgeNegative = nan(nK, 1);
posRows = find(omegaT >= 0);
negRows = find(omegaT <= 0);
for j = find(activeKMask).'
    [~, ip] = max(powerNorm(posRows, j), [], 1);
    [~, in] = max(powerNorm(negRows, j), [], 1);
    ridgePositive(j) = omegaT(posRows(ip));
    ridgeNegative(j) = omegaT(negRows(in));
end

%==========================================================================
% Output
%==========================================================================
result.spectralDb    = spectralDb;
result.omegaT        = omegaT;
result.kGrid         = kGrid;
result.kNormalized   = kNormalized;
result.power         = powerNorm;
result.powerRaw      = powerRaw;
result.rawSpectralPower = specNative;
result.rawOmegaGridT    = rawOmegaGridT;
result.ridgePositive = ridgePositive;
result.ridgeNegative = ridgeNegative;
result.activeKMask   = activeKMask;
result.rawColumnPower = rawColumnPower;
result.dftBinSpacingOmegaT = 2*pi/nativePeriodCount;
result.effectiveFreqResolutionOmegaT = temporalENBW*2*pi/nativePeriodCount;
result.zeroPaddedGridSpacingOmegaT = 2*pi/nFreqOut;
result.spatialDftBinSpacingK = 2*pi/(Nx*dx);
result.spatialDftBinSpacingKOverK0 = T/(Nx*dx);
result.spatialEffectiveResolutionKOverK0 = spatialENBW*T/(Nx*dx);
result.spatialROI = spatialROI;
result.timeROI = timeROI;
result.roiXStart = xROI(1);
result.roiXEnd = xROI(end);
result.roiXLength = xROI(end) - xROI(1);
result.roiTStart = tROI(1);
result.roiTEnd = tROI(end);
result.roiTLength = tROI(end) - tROI(1);
result.roiPhaseAtStart = mod(tROI(1), T)/T;
result.temporalWindow = temporalWindow;
result.temporalWindowGains = struct('coherentGain', temporalCG, ...
    'energyGain', temporalEG, 'enbwBins', temporalENBW);
result.spatialWindow = spatialWindow;
result.spatialWindowGains = struct('coherentGain', spatialCG, ...
    'energyGain', spatialEG, 'enbwBins', spatialENBW);
result.timeEndpointDropped = timeEndpointDropped;
result.spatialEndpointDropped = spatialEndpointDropped;
result.fieldComponent = fieldComponent;
result.interfaceSideConvention = ...
    'E sampled on integer Yee nodes with the after-switch constitutive value D/eps(t+) at ideal switches';
result.zeroPaddingFactor = zeroPaddingFactor;
result.temporalPeriod = T;
result.spatialPeriod  = Lambda;
result.measurementROI = measurementROI;
result.timeWindowStart = timeWindowStart;
result.modulationPhase = modulationPhase;
result.nativePeriodCount = nativePeriodCount;
result.nFold = nFreqOut;
end

% -------------------------------------------------------------------------
function perm = centered_index(N)
%CENTERED_INDEX 1-based reorder such that new(j) <-> old(perm(j)) gives the
%signed centered grid (-floor(N/2):ceil(N/2)-1) from unshifted bins 0:N-1.
%Replaces fftshift and is parity-correct for odd and even N.
perm = mod((-floor(N/2):ceil(N/2)-1), N) + 1;
end

% -------------------------------------------------------------------------
function w = make_window(kind, tukeyR, M, orientation)
%MAKE_WINDOW Build a Hann/rect/Tukey window of length M.
switch lower(kind)
    case 'hann'
        % Periodic (DFT-form) Hann: w(n) = 0.5 - 0.5*cos(2*pi*n/M) has zero
        % leakage on bin edges for the length-M DFT.
        w = 0.5 - 0.5*cos(2*pi*(0:M-1)'/M);
    case 'rect'
        w = ones(M, 1);
    case 'tukey'
        if ~isscalar(tukeyR) || ~isreal(tukeyR) || ~isfinite(tukeyR) || ...
                tukeyR < 0 || tukeyR > 1
            error('Tukey taper parameter must be a finite real scalar in [0,1].');
        end
        w = tukey_window(M, tukeyR);   % base MATLAB, no Toolbox dependency
    otherwise
        error('Unsupported window type ''%s''.', kind);
end
if strcmp(orientation, 'row')
    w = w.';
end
end

% -------------------------------------------------------------------------
function w = tukey_window(M, r)
%TUKEY_WINDOW Base-MATLAB Tukey window, matching the Signal Processing
%Toolbox tukeywin(M, r): r = 0 -> rectangular, r = 1 -> Hann, 0 < r < 1 ->
%cosine-tapered edges around a flat middle. Implemented locally so the
%package needs no Toolbox (README: "Base MATLAB, no Toolbox dependency").
%The right taper must be MIRRORED (flipud): assigning w(1:nr+1) directly to
%w(M-nr:M) puts a second RISING taper at the end and never tapers to zero,
%and any w(w==0)=1 fill would clobber the genuine zero taper endpoints.
if M <= 1
    w = ones(M, 1);
    return;
end
if r <= 0
    w = ones(M, 1);
elseif r >= 1
    w = 0.5 - 0.5*cos(2*pi*(0:M-1)'/(M-1));   % symmetric Hann
else
    N = M - 1;
    nr = floor(r*N/2);                   % taper samples per side (0-based)
    m = (0:nr).';
    w = zeros(M, 1);
    w(1:nr+1) = 0.5*(1 + cos(pi*(-1 + 2*m/(r*N))));  % left taper 0 -> 1
    w(M-nr:M) = flipud(w(1:nr+1));       % right taper mirrors down to 0
    w(nr+2:M-nr-1) = 1;                  % flat middle (empty if tapers meet)
end
end

% -------------------------------------------------------------------------
function [kind, tukeyR] = parse_window(fftCfg, kindField, rField, defaultKind)
%PARSE_WINDOW Read a window type + optional Tukey taper from fftCfg.
if isfield(fftCfg, kindField) && ~isempty(fftCfg.(kindField))
    kind = lower(char(fftCfg.(kindField)));
else
    kind = defaultKind;
end
if ~ismember(kind, {'hann','rect','tukey'})
    error('fftCfg.%s must be ''hann'', ''rect'', or ''tukey''.', kindField);
end
if isfield(fftCfg, rField) && ~isempty(fftCfg.(rField))
    tukeyR = fftCfg.(rField);
else
    tukeyR = 0.5;
end
end

% -------------------------------------------------------------------------
function [cg, eg, enbw] = window_gains(w)
%WINDOW_GAINS Coherent gain, energy gain, and ENBW (in bins) of a window.
M = numel(w);
cg = sum(w)/M;
eg = sum(w.^2)/M;
if sum(w) ~= 0
    enbw = (M*sum(w.^2))/sum(w)^2;
else
    enbw = nan;
end
end

% -------------------------------------------------------------------------
function v = logical_scalar(x, label)
%LOGICAL_SCALAR Strict scalar-logical validation (0/1 or logical).
if ~islogical(x) && ~(isnumeric(x) && isscalar(x) && (x == 0 || x == 1))
    error('%s must be a scalar logical (true/false or 0/1).', label);
end
v = logical(x);
end
