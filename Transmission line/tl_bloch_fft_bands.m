function result = tl_bloch_fft_bands(model,kScan,blochCfg)
%TL_BLOCH_FFT_BANDS Reconstruct bulk Floquet bands from stroboscopic data.
%
% Every k is propagated independently from a complete energy-scaled state
% basis.  Observables are sampled at t=p*T and transformed with an inverse
% FFT so exp(-i*omega*p*T) peaks at positive omega under the repository's
% exp(i*k*x-i*omega*t) convention.  This reconstructs Re(omega) in the
% first temporal Floquet zone; it is not an eigenvalue solver.

if nargin ~= 3 || ~isstruct(model) || ~isscalar(model) || ...
        ~isstruct(blochCfg) || ~isscalar(blochCfg)
    error('Use tl_bloch_fft_bands(model,kScan,blochCfg).');
end
if ~isnumeric(kScan) || isempty(kScan) || ~isvector(kScan) || ...
        ~isreal(kScan) || any(~isfinite(kScan))
    error('kScan must be a nonempty finite real vector.');
end
kScan = kScan(:).';
if any(abs(kScan) > pi/model.cell.a*(1+100*eps))
    error('Bulk FFT kScan must lie in the spatial first Brillouin zone.');
end

stepsPerPeriod = tl_option('integer',blochCfg,'stepsPerPeriod',256,4);
temporalCellCount = tl_option('integer',blochCfg,'temporalCellCount',64,4);
zeroPaddingFactor = tl_option('integer',blochCfg,'zeroPaddingFactor',4,1);
dynamicRangeDb = tl_option('positive',blochCfg,'dynamicRangeDb',60);
returnTemporalSamples = tl_option('logical', ...
    blochCfg,'returnTemporalSamples',false);
T = model.modulation.period;
Omega = model.modulation.OmegaRadPerSec;
nState = model.bulk.stateDimension;
nK = numel(kScan);

[segmentMidpoints,segmentDurations,integrationName] = ...
    tl_time_layers(model,stepsPerPeriod);

sampleCount = temporalCellCount;
observableSamples = complex(zeros(sampleCount,nState,nState,nK));
dominantGrowthPerPeriod = zeros(1,nK);
monodromy = complex(zeros(nState,nState,nK));
sampleIndex = (0:sampleCount-1).';

for ik = 1:nK
    U = eye(nState);
    for layer = 1:numel(segmentDurations)
        A = model.functions.bulkStateMatrix( ...
            kScan(ik),segmentMidpoints(layer));
        U = expm(A*segmentDurations(layer))*U;
    end
    monodromy(:,:,ik) = U;
    state = energy_scaled_basis(model,kScan(ik));
    O = model.functions.observableMatrix(kScan(ik),0);
    observableSamples(1,:,:,ik) = reshape(O*state,[1 nState nState]);
    amplitude = zeros(sampleCount,1);
    amplitude(1) = norm(O*state,'fro');
    for sample = 2:sampleCount
        state = U*state;
        values = O*state;
        if any(~isfinite(values),'all')
            error(['Bulk stroboscopic states overflowed at k=%.16g. ' ...
                'Reduce temporalCellCount or narrow the scan.'],kScan(ik));
        end
        observableSamples(sample,:,:,ik) = ...
            reshape(values,[1 nState nState]);
        amplitude(sample) = norm(values,'fro');
    end
    fitRows = floor(sampleCount/2)+1:sampleCount;
    coordinate = sampleIndex(fitRows)-mean(sampleIndex(fitRows));
    denominator = sum(coordinate.^2);
    dominantGrowthPerPeriod(ik) = coordinate.'* ...
        log(max(amplitude(fitRows),realmin))/denominator;
end

window = periodic_hann(sampleCount);
nFft = zeroPaddingFactor*sampleCount;
windowed = observableSamples.*reshape(window,[],1,1,1);
complexAmplitudeScale = nFft/sum(window);
complexSpectrum = ifft(windowed,nFft,1)*complexAmplitudeScale;
powerByChannel = abs(complexSpectrum).^2;
rawPower = reshape(sum(sum(powerByChannel,2),3),nFft,nK);

frequencyOrders = (-floor(nFft/2):ceil(nFft/2)-1).';
centeredIndices = mod(frequencyOrders,nFft)+1;
rawPower = rawPower(centeredIndices,:);
omega = frequencyOrders*Omega/nFft;
columnMaximum = max(rawPower,[],1);
if any(~isfinite(columnMaximum)) || any(columnMaximum <= 0)
    error('Bulk FFT returned an empty or non-finite spectrum column.');
end
globalMaximum = max(columnMaximum);
globalPower = rawPower/globalMaximum;
columnPower = rawPower./columnMaximum;
floorPower = 10^(-dynamicRangeDb/10);
spectrumDb = 10*log10(max(columnPower,floorPower));
[ridgePower,ridgeIndex] = max(rawPower,[],1);

result.k = kScan;
result.kaOverPi = kScan*model.cell.a/pi;
result.omega = omega;
result.omegaOverOmega = omega/Omega;
result.rawPower = rawPower;
result.globalNormalizedPower = globalPower;
result.columnNormalizedPower = columnPower;
result.spectrumDb = spectrumDb;
result.ridgeOmega = omega(ridgeIndex).';
result.ridgePower = ridgePower;
result.monodromy = monodromy;
result.dominantGrowthPerPeriod = dominantGrowthPerPeriod;
result.dominantImagOmega = dominantGrowthPerPeriod/T;
result.windowName = 'periodic-hann';
result.windowCoherentGain = mean(window);
result.windowEnergyGain = mean(window.^2);
result.windowEnbwBins = sampleCount*sum(window.^2)/sum(window)^2;
result.spectralNormalization = struct( ...
    'name','window-coherent-amplitude-squared', ...
    'complexAmplitudeScale',complexAmplitudeScale, ...
    'powerMeaning','squared coherent amplitude sampled on the FFT grid', ...
    'frequencyIntegralRequiresBinWidth',true);
result.nativeOmegaResolution = Omega/sampleCount;
result.zeroPaddedOmegaSpacing = Omega/nFft;
result.temporalCellCount = temporalCellCount;
result.stepsPerPeriod = stepsPerPeriod;
result.nFft = nFft;
result.integrationName = integrationName;
result.segmentMidpoints = segmentMidpoints;
result.segmentDurations = segmentDurations;
result.modelSnapshot = model.snapshot;
result.blochConfig = blochCfg;
if returnTemporalSamples
    result.observableSamples = observableSamples;
    result.sampleCellIndices = sampleIndex;
end
end

% -------------------------------------------------------------------------
function basis = energy_scaled_basis(model,k)
C = model.functions.capacitanceBulk(0);
[~,~,Lk] = model.functions.bulkMatrices(k);
if strcmp(model.kind,'sspp')
    scales = [sqrt(C),sqrt(Lk)];
elseif isinf(model.resonator.Cblock)
    scales = [sqrt(C),sqrt(Lk),sqrt(model.resonator.L0)];
else
    scales = [sqrt(C),sqrt(Lk),sqrt(model.resonator.L0), ...
        sqrt(model.resonator.Cblock)];
end
basis = diag(scales);
end

% -------------------------------------------------------------------------
function window = periodic_hann(count)
index = (0:count-1).';
window = 0.5-0.5*cos(2*pi*index/count);
end
