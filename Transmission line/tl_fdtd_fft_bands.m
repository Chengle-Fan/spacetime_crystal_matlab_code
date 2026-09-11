function result = tl_fdtd_fft_bands(probeSignals,time,kScan,fftCfg)
%TL_FDTD_FFT_BANDS Reconstruct finite-chain bands from fixed V(t) probes.
%
%   result = tl_fdtd_fft_bands(probeSignals,time,kScan,fftCfg)
%
% probeSignals is normally Nt-by-Nprobe-by-Nk.  Because MATLAB drops a
% trailing singleton dimension, a two-dimensional array is Nt-by-Nprobe
% when nK=1 and Nt-by-Nk (one probe) when nK>1.  The routine applies an
% endpoint-free integer-period
% time gate, a periodic Hann window, a full physical-frequency FFT, and an
% explicit incoherent power sum of replicas separated by Omega into
% [-Omega/2,Omega/2).  It never advances a Bloch mode or calls PWE/TMM.
%
% Required fftCfg fields: temporalPeriod, analysisTimeRange=[start end].
% The horizontal coordinate is the Gaussian source centre k_c supplied by
% the caller.  It is not an exact wave number for every packet component.

if nargin ~= 4 || ~isstruct(fftCfg) || ~isscalar(fftCfg)
    error('Use tl_fdtd_fft_bands(probeSignals,time,kScan,fftCfg).');
end
required = {'temporalPeriod','analysisTimeRange'};
missing = required(~isfield(fftCfg,required));
if ~isempty(missing)
    error('fftCfg is missing field(s): %s.',strjoin(missing,', '));
end

time = validate_time(time);
dtRecord = time(2)-time(1);
kScan = validate_k_scan(kScan);
nK = numel(kScan);
probeSignals = validate_signals(probeSignals,numel(time),nK);
nProbe = size(probeSignals,2);

T = tl_option('positive',fftCfg,'temporalPeriod',[]);
Omega = 2*pi/T;
samplesPerPeriodFloat = T/dtRecord;
samplesPerPeriod = round(samplesPerPeriodFloat);
if samplesPerPeriod < 2 || abs(samplesPerPeriodFloat-samplesPerPeriod) > ...
        1e-9*max(1,samplesPerPeriod)
    error(['The probe sampling interval must be commensurate with the ' ...
        'modulation period and provide at least two samples per period.']);
end

[analysisRows,analysisRange,analysisPeriodCount] = ...
    analysis_gate(time,dtRecord,T,fftCfg.analysisTimeRange);
nTime = numel(analysisRows);
if nTime ~= analysisPeriodCount*samplesPerPeriod
    error('The integer-period analysis gate has an inconsistent sample count.');
end
signals = double(probeSignals(analysisRows,:,:));

zeroPaddingFactor = tl_option('integer',fftCfg,'zeroPaddingFactor',4,1);
dynamicRangeDb = tl_option('positive',fftCfg,'dynamicRangeDb',60);
activeThreshold = tl_option('fraction-or-zero', ...
    fftCfg,'activeColumnRelativeThreshold',1e-12);
returnProbePower = tl_option('logical',fftCfg,'returnProbePower',false);
returnProbeSignals = tl_option('logical',fftCfg,'returnProbeSignals',false);

sampleIndex = (0:nTime-1).';
window = 0.5-0.5*cos(2*pi*sampleIndex/nTime);
nFft = zeroPaddingFactor*nTime;
windowed = signals.*reshape(window,[],1,1);

% Under exp(i*k*x-i*omega*t), ifft places positive omega on positive bins.
% Probe powers are added, not complex amplitudes, so a field node at one
% probe cannot coherently cancel a branch visible at another probe.
% Undo ifft's 1/nFft and correct the Hann coherent gain.  A grid-aligned
% complex tone then has zero-padding-invariant peak amplitude.
complexAmplitudeScale = nFft/sum(window);
complexSpectrum = ifft(windowed,nFft,1)*complexAmplitudeScale;
rawProbePowerUnshifted = abs(complexSpectrum).^2;
rawPowerUnshifted = reshape(sum(rawProbePowerUnshifted,2),nFft,nK);

% Accumulate every sampled physical-frequency replica modulo one Omega.
% Folding unshifted integer DFT bins avoids odd/even half-zone errors.
nFold = zeroPaddingFactor*analysisPeriodCount;
if mod(nFft,nFold) ~= 0 || nFft/nFold ~= samplesPerPeriod
    error('The FFT and Floquet folding grids are not commensurate.');
end
[foldedPower,foldedOrders] = tl_fold_power( ...
    rawPowerUnshifted,0:nFft-1,nFold);
foldedOmega = foldedOrders*Omega/nFold;

rawOrders = (-floor(nFft/2):ceil(nFft/2)-1).';
rawIndices = mod(rawOrders,nFft)+1;
rawPower = rawPowerUnshifted(rawIndices,:);
rawOmega = rawOrders*2*pi/(nFft*dtRecord);

rawColumnPower = max(foldedPower,[],1);
globalMaximum = max(rawColumnPower);
if ~isfinite(globalMaximum) || globalMaximum <= 0
    error('The probe FFT returned empty or non-finite spectral power.');
end
activeKMask = isfinite(rawColumnPower) & rawColumnPower > 0 & ...
    rawColumnPower >= activeThreshold*globalMaximum;
columnNormalizedPower = zeros(size(foldedPower));
columnNormalizedPower(:,activeKMask) = ...
    foldedPower(:,activeKMask)./rawColumnPower(activeKMask);
globalNormalizedPower = foldedPower/globalMaximum;
displayFloor = 10^(-dynamicRangeDb/10);
if displayFloor == 0
    error('fftCfg.dynamicRangeDb exceeds double-precision display range.');
end
spectrumDb = -dynamicRangeDb*ones(size(foldedPower));
spectrumDb(:,activeKMask) = 10*log10(max( ...
    columnNormalizedPower(:,activeKMask),displayFloor));

result.k = kScan;
result.foldedOmega = foldedOmega;
result.omegaOverOmega = foldedOmega/Omega;
result.spectrumDb = spectrumDb;
result.foldedPower = foldedPower;
result.rawColumnPower = rawColumnPower;
result.activeKMask = activeKMask;
result.columnNormalizedPower = columnNormalizedPower;
result.globalNormalizedPower = globalNormalizedPower;
result.rawOmega = rawOmega;
result.rawOmegaOverOmega = rawOmega/Omega;
result.rawPower = rawPower;
result.probeCount = nProbe;
result.recordTimeStep = dtRecord;
result.samplesPerPeriod = samplesPerPeriod;
result.analysisTimeRange = analysisRange;
result.analysisRows = analysisRows;
result.analysisPeriodCount = analysisPeriodCount;
result.temporalPeriod = T;
result.nFft = nFft;
result.nFoldedFrequency = nFold;
result.nativeOmegaResolution = 2*pi/(nTime*dtRecord);
result.nativeOmegaResolutionNormalized = 1/analysisPeriodCount;
result.zeroPaddedOmegaSpacing = 2*pi/(nFft*dtRecord);
result.displayOmegaSpacingNormalized = 1/nFold;
result.window = struct('name','periodic-hann', ...
    'coherentGain',mean(window), ...
    'energyGain',mean(window.^2), ...
    'enbwBins',mean(window.^2)/mean(window)^2);
result.spectralNormalization = struct( ...
    'name','window-coherent-amplitude-squared', ...
    'complexAmplitudeScale',complexAmplitudeScale, ...
    'powerMeaning','squared coherent amplitude sampled on the FFT grid', ...
    'frequencyIntegralRequiresBinWidth',true);
result.activeColumnRelativeThreshold = activeThreshold;
result.fftConfig = fftCfg;
result.interpretation = ['Finite-chain, finite-source, finite-window real-' ...
    'frequency response. Linewidth is not Im(omega).'];
if returnProbePower
    result.rawProbePower = rawProbePowerUnshifted(rawIndices,:,:);
    result.foldedProbePower = tl_fold_power( ...
        rawProbePowerUnshifted,0:nFft-1,nFold);
end
if returnProbeSignals
    result.probeSignals = probeSignals;
    result.time = time;
end
end

% -------------------------------------------------------------------------
function time = validate_time(value)
if ~isnumeric(value) || ~isvector(value) || numel(value) < 4 || ...
        ~isreal(value) || any(~isfinite(value))
    error('time must contain at least four finite real samples.');
end
time = value(:);
steps = diff(time);
if any(steps <= 0)
    error('time must be strictly increasing.');
end
dt = steps(1);
timeRoundoff = 128*eps(max(abs(time)));
if max(abs(steps-dt)) > max(1e-9*dt,timeRoundoff)
    error('time must be uniformly sampled.');
end
end

function k = validate_k_scan(value)
if ~isnumeric(value) || isempty(value) || ~isvector(value) || ...
        ~isreal(value) || any(~isfinite(value))
    error('kScan must be a nonempty finite real vector.');
end
k = value(:).';
end

function signals = validate_signals(value,nTime,nK)
if ~isnumeric(value) || isempty(value) || any(~isfinite(value(:)))
    error('probeSignals must be a finite nonempty numeric array.');
end
if ismatrix(value)
    if size(value,1) ~= nTime
        error('probeSignals must have one row per time sample.');
    end
    % MATLAB drops a trailing singleton dimension.  For nK=1, an
    % Nt-by-Nprobe-by-1 array therefore arrives as Nt-by-Nprobe and must not
    % be mistaken for several k columns.
    if nK == 1
        signals = reshape(value,nTime,size(value,2),1);
    elseif size(value,2) == nK
        signals = reshape(value,nTime,1,nK);
    else
        error('For multiple k values, 2-D probeSignals must be Nt-by-Nk.');
    end
elseif ndims(value) == 3
    if size(value,1) ~= nTime || size(value,3) ~= nK
        error('A three-dimensional probeSignals input must be Nt-by-Nprobe-by-Nk.');
    end
    signals = value;
else
    error('probeSignals must be Nt-by-Nk or Nt-by-Nprobe-by-Nk.');
end
end

function [rows,timeRange,periodCount] = analysis_gate(time,dt,T,value)
if ~isnumeric(value) || ~isvector(value) || numel(value) ~= 2 || ...
        ~isreal(value) || any(~isfinite(value)) || value(2) <= value(1)
    error('fftCfg.analysisTimeRange must contain two increasing finite times.');
end
timeRange = value(:).';
tolerance = max(1e-7*dt,128*eps(max(abs(time))));
if timeRange(1) < time(1)-tolerance || timeRange(2) > time(end)+tolerance
    error('fftCfg.analysisTimeRange lies outside the recorded time range.');
end
nodeFloat = (timeRange-time(1))/dt;
nodes = round(nodeFloat);
nodeTolerance = max(1e-8,128*eps(max(abs(nodeFloat))));
if any(abs(nodeFloat-nodes) > nodeTolerance)
    error('Both analysisTimeRange endpoints must lie on the probe time grid.');
end
periodCountFloat = (timeRange(2)-timeRange(1))/T;
periodCount = round(periodCountFloat);
if periodCount < 2 || abs(periodCountFloat-periodCount) > ...
        1e-9*max(1,periodCount)
    error('The FFT analysis gate must span at least two complete periods.');
end
% MATLAB row nodes are one-based.  [start,end) excludes the repeated end.
rows = (nodes(1)+1):nodes(2);
end
