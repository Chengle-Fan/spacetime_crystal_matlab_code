function result = tl_xt_fft_bands(observation,referenceObservation,fftCfg)
%TL_XT_FFT_BANDS Reconstruct finite-chain bands with a signed x-t FFT.
%
%   result = tl_xt_fft_bands(observation,referenceObservation,fftCfg)
%
% observation and referenceObservation are independent simulations or
% measurements with fields:
%   data, x, t, cellIndex, channelId, channelOffset, observableName
% data is Nt-by-Nsample.  Samples belonging to each channel must cover the
% same consecutive cells once.  The reference is used only to establish
% source support; this routine never asks PWE/TMM where a peak should lie.
%
% Required fftCfg fields: cellPeriod, temporalPeriod.
% maximumPhysicalOmega optionally limits which raw physical-frequency bins
% are accumulated during Floquet folding; the complete raw spectrum is
% still returned for audit.

if nargin ~= 3 || ~isstruct(fftCfg) || ~isscalar(fftCfg)
    error(['Use tl_xt_fft_bands(observation,referenceObservation,' ...
        'fftCfg).']);
end
a = read_positive(fftCfg,'cellPeriod',[]);
T = read_positive(fftCfg,'temporalPeriod',[]);
zeroPaddingTime = read_integer(fftCfg,'zeroPaddingTime',2,1);
zeroPaddingSpace = read_integer(fftCfg,'zeroPaddingSpace',2,1);
dynamicRangeDb = read_positive(fftCfg,'dynamicRangeDb',60);
referenceThresholdDb = read_positive(fftCfg,'referenceThresholdDb',40);
ridgeThresholdDb = read_positive(fftCfg,'ridgeThresholdDb',30);
ridgeCount = read_integer(fftCfg,'ridgeCount',2,1);
ridgeJumpBins = read_integer(fftCfg,'ridgeJumpBins',6,1);
ridgeExclusionBins = read_integer(fftCfg,'ridgeExclusionBins',3,0);
removeTemporalMean = read_logical(fftCfg,'removeTemporalMean',false);

obs = validate_observation(observation,a,'observation');
ref = validate_observation(referenceObservation,a,'referenceObservation');
assert_compatible(obs,ref);
dt = obs.t(2)-obs.t(1);
nTime = numel(obs.t);
nCell = numel(obs.cells);
Omega = 2*pi/T;
periodCountFloat = nTime*dt/T;
periodCount = round(periodCountFloat);
if periodCount < 1 || abs(periodCountFloat-periodCount) > ...
        200*eps(max(1,abs(periodCountFloat)))
    error(['The endpoint-free time record Nt*dt must contain an integer ' ...
        'number of modulation periods.']);
end
if isfield(fftCfg,'maximumPhysicalOmega') && ...
        ~isempty(fftCfg.maximumPhysicalOmega)
    maximumPhysicalOmega = read_positive( ...
        fftCfg,'maximumPhysicalOmega',[]);
    if maximumPhysicalOmega >= pi/dt
        error(['fftCfg.maximumPhysicalOmega must be below the recorded ' ...
            'Nyquist angular frequency pi/dt.']);
    end
else
    maximumPhysicalOmega = pi/dt;
end

nFftTime = zeroPaddingTime*nTime;
nFftSpace = zeroPaddingSpace*nCell;
timeWindow = periodic_hann(nTime);
spaceWindow = periodic_hann(nCell).';
[rawChannelPower,timeOrders,kOrders] = transform_channels( ...
    obs.dataGrid,timeWindow,spaceWindow,nFftTime,nFftSpace, ...
    removeTemporalMean);
[referenceRawChannelPower,~,~] = transform_channels( ...
    ref.dataGrid,timeWindow,spaceWindow,nFftTime,nFftSpace, ...
    removeTemporalMean);

rawOmega = 2*pi*timeOrders/(nFftTime*dt);
rawK = 2*pi*kOrders/(nFftSpace*a);
foldingPhysicalFrequencyMask = abs(rawOmega) <= ...
    maximumPhysicalOmega+128*eps(maximumPhysicalOmega);
foldingChannelPower = rawChannelPower;
foldingReferenceChannelPower = referenceRawChannelPower;
foldingChannelPower(~foldingPhysicalFrequencyMask,:,:) = 0;
foldingReferenceChannelPower(~foldingPhysicalFrequencyMask,:,:) = 0;
foldBinCount = round(Omega/(2*pi/(nFftTime*dt)));
if foldBinCount < 1 || abs(foldBinCount- ...
        zeroPaddingTime*periodCount) > 0
    error('The time grid is not commensurate with the Floquet folding grid.');
end
[foldedChannelPower,foldedOrders] = fold_temporal_power( ...
    foldingChannelPower,timeOrders,foldBinCount);
[referenceFoldedChannelPower,referenceOrders] = fold_temporal_power( ...
    foldingReferenceChannelPower,timeOrders,foldBinCount);
if ~isequal(foldedOrders,referenceOrders)
    error('Internal reference folding mismatch.');
end
foldedOmega = foldedOrders*Omega/foldBinCount;

rawPower = sum(rawChannelPower,3);
referenceRawPower = sum(referenceRawChannelPower,3);
foldedPower = sum(foldedChannelPower,3);
referenceFoldedPower = sum(referenceFoldedChannelPower,3);
referenceMaximum = max(referenceFoldedPower,[],'all');
if ~isfinite(referenceMaximum) || referenceMaximum <= 0
    error('referenceObservation produced no finite nonzero spectral power.');
end
sourceSupportMask = referenceFoldedPower >= ...
    referenceMaximum*10^(-referenceThresholdDb/10);

globalMaximum = max(foldedPower,[],'all');
if ~isfinite(globalMaximum) || globalMaximum <= 0
    error('observation produced no finite nonzero spectral power.');
end
globalNormalizedPower = foldedPower/globalMaximum;
columnMaximum = max(foldedPower.*sourceSupportMask,[],1);
columnNormalizedPower = zeros(size(foldedPower));
supportedColumns = columnMaximum > 0;
columnNormalizedPower(:,supportedColumns) = ...
    foldedPower(:,supportedColumns)./columnMaximum(supportedColumns);
displayFloor = 10^(-dynamicRangeDb/10);
spectrumDb = 10*log10(max(columnNormalizedPower,displayFloor));
spectrumDb(:,~supportedColumns) = -dynamicRangeDb;

significanceMask = false(size(foldedPower));
significanceMask(:,supportedColumns) = ...
    foldedPower(:,supportedColumns) >= columnMaximum(supportedColumns)* ...
    10^(-ridgeThresholdDb/10);
ridgeMask = sourceSupportMask & significanceMask;
[ridgeIndices,ridgePower] = extract_ridges(foldedPower,ridgeMask, ...
    ridgeCount,ridgeJumpBins,ridgeExclusionBins);
ridgeOmega = nan(size(ridgeIndices));
validRidge = ~isnan(ridgeIndices);
ridgeOmega(validRidge) = foldedOmega(ridgeIndices(validRidge));

result.k = rawK;
result.kaOverPi = rawK*a/pi;
result.rawOmega = rawOmega;
result.rawPower = rawPower;
result.rawChannelPower = rawChannelPower;
result.referenceRawPower = referenceRawPower;
result.referenceRawChannelPower = referenceRawChannelPower;
result.foldingPhysicalFrequencyMask = foldingPhysicalFrequencyMask;
result.foldedOmega = foldedOmega;
result.omegaOverOmega = foldedOmega/Omega;
result.foldedPower = foldedPower;
result.foldedChannelPower = foldedChannelPower;
result.referenceFoldedPower = referenceFoldedPower;
result.referenceFoldedChannelPower = referenceFoldedChannelPower;
result.sourceSupportMask = sourceSupportMask;
result.globalNormalizedPower = globalNormalizedPower;
result.columnNormalizedPower = columnNormalizedPower;
result.spectrumDb = spectrumDb;
result.ridgeIndices = ridgeIndices;
result.ridgeOmega = ridgeOmega;
result.ridgePower = ridgePower;
result.channelIds = obs.channelIds;
result.observableName = obs.observableName;
result.nativeOmegaResolution = 2*pi/(nTime*dt);
result.nativeKResolution = 2*pi/(nCell*a);
result.zeroPaddedOmegaSpacing = 2*pi/(nFftTime*dt);
result.zeroPaddedKSpacing = 2*pi/(nFftSpace*a);
result.maximumPhysicalOmega = maximumPhysicalOmega;
result.nyquistOmega = pi/dt;
result.periodCount = periodCount;
result.window = struct( ...
    'timeName','periodic-hann','spaceName','periodic-hann', ...
    'timeCoherentGain',mean(timeWindow), ...
    'timeEnergyGain',mean(timeWindow.^2), ...
    'timeEnbwBins',nTime*sum(timeWindow.^2)/sum(timeWindow)^2, ...
    'spaceCoherentGain',mean(spaceWindow), ...
    'spaceEnergyGain',mean(spaceWindow.^2), ...
    'spaceEnbwBins',nCell*sum(spaceWindow.^2)/sum(spaceWindow)^2);
result.spectralNormalization = struct( ...
    'name','window-coherent-amplitude-squared', ...
    'timeComplexAmplitudeScale',nFftTime/sum(timeWindow), ...
    'spaceComplexAmplitudeScale',1/sum(spaceWindow), ...
    'powerMeaning','squared coherent amplitude sampled on the k-omega grid', ...
    'spectralIntegralRequiresBinArea',true);
result.observationMetadata = obs.metadata;
result.referenceMetadata = ref.metadata;
result.fftConfig = fftCfg;
result.interpretation = ['The x-t ridge estimates Re(omega) only.  Its ' ...
    'linewidth is not Im(omega), and finite-chain growth includes port ' ...
    'escape, beating, and multimode effects.'];
end

% -------------------------------------------------------------------------
function obs = validate_observation(input,a,label)
if ~isstruct(input) || ~isscalar(input)
    error('%s must be a scalar struct.',label);
end
required = {'data','x','t','cellIndex','channelId','channelOffset', ...
    'observableName'};
missing = required(~isfield(input,required));
if ~isempty(missing)
    error('%s is missing field(s): %s.',label,strjoin(missing,', '));
end
data = input.data;
t = input.t(:);
x = input.x(:);
cellIndex = input.cellIndex(:);
channelOffset = input.channelOffset(:);
if ~isnumeric(data) || ~ismatrix(data) || isempty(data) || ...
        any(~isfinite(data),'all')
    error('%s.data must be a finite nonempty numeric matrix.',label);
end
if ~isnumeric(t) || ~isreal(t) || numel(t) < 4 || any(~isfinite(t)) || ...
        any(diff(t) <= 0)
    error('%s.t must contain at least four increasing finite samples.',label);
end
dt = diff(t);
timeTolerance = max(1e-9*dt(1),128*eps(max(abs(t))));
if max(abs(dt-dt(1))) > timeTolerance
    error('%s.t must be uniformly spaced.',label);
end
if size(data,1) ~= numel(t)
    error('%s.data must have one row per time sample.',label);
end
nSample = size(data,2);
if ~isnumeric(x) || ~isreal(x) || numel(x) ~= nSample || ...
        any(~isfinite(x))
    error('%s.x must have one finite real value per data column.',label);
end
if ~isnumeric(cellIndex) || ~isreal(cellIndex) || ...
        numel(cellIndex) ~= nSample || any(~isfinite(cellIndex)) || ...
        any(cellIndex ~= round(cellIndex))
    error('%s.cellIndex must contain one finite integer per sample.',label);
end
if ~isnumeric(channelOffset) || ~isreal(channelOffset) || ...
        numel(channelOffset) ~= nSample || any(~isfinite(channelOffset))
    error('%s.channelOffset must contain one finite real value per sample.',label);
end
channelId = input.channelId;
if isscalar(channelId) && nSample > 1
    if isnumeric(channelId) || isstring(channelId) || ischar(channelId)
        channelId = repmat_scalar_id(channelId,nSample);
    end
end
if numel(channelId) ~= nSample
    error('%s.channelId must contain one identifier per sample.',label);
end
[channelIds,groups] = group_ids(channelId);
nChannel = numel(channelIds);
cells = unique(cellIndex).';
if numel(cells) < 4 || any(diff(cells) ~= 1)
    error('%s must cover at least four consecutive cell indices.',label);
end
nCell = numel(cells);
dataGrid = complex(zeros(numel(t),nCell,nChannel));
for channel = 1:nChannel
    members = find(groups == channel);
    [channelCells,order] = sort(cellIndex(members));
    members = members(order);
    if ~isequal(channelCells(:).',cells)
        error('%s channel %d must sample every ROI cell exactly once.', ...
            label,channel);
    end
    offsets = channelOffset(members);
    if max(abs(offsets-offsets(1))) > 100*eps(max(1,max(abs(offsets))))
        error('%s.channelOffset must be constant within each channel.',label);
    end
    residual = x(members)-a*channelCells-offsets;
    if max(abs(residual-residual(1))) > 1e-9*max(a,1e-12)
        error('%s.x is inconsistent with cellIndex, channelOffset, and a.',label);
    end
    dataGrid(:,:,channel) = data(:,members);
end
observableName = input.observableName;
if isstring(observableName) && isscalar(observableName)
    observableName = char(observableName);
end
if ~ischar(observableName) || size(observableName,1) ~= 1
    error('%s.observableName must be a text scalar.',label);
end
metadata = input;
metadata = rmfield(metadata,'data');
obs = struct('dataGrid',dataGrid,'t',t,'x',x,'cells',cells, ...
    'channelIds',{channelIds},'observableName',observableName, ...
    'metadata',metadata);
end

% -------------------------------------------------------------------------
function output = repmat_scalar_id(value,count)
if ischar(value)
    output = repmat({value},count,1);
elseif isstring(value)
    output = repmat(value,count,1);
else
    output = repmat(value,count,1);
end
end

% -------------------------------------------------------------------------
function [ids,groups] = group_ids(values)
if ischar(values)
    values = cellstr(values);
end
if isnumeric(values) || islogical(values) || isstring(values) || iscellstr(values)
    [ids,~,groups] = unique(values,'stable');
else
    error('channelId must be numeric, logical, string, char, or cellstr.');
end
if isstring(ids), ids = cellstr(ids); end
if isnumeric(ids) || islogical(ids), ids = num2cell(ids); end
ids = ids(:).';
groups = groups(:);
end

% -------------------------------------------------------------------------
function assert_compatible(obs,ref)
if ~isequal(obs.t,ref.t) || ~isequal(obs.cells,ref.cells) || ...
        size(obs.dataGrid,3) ~= size(ref.dataGrid,3) || ...
        ~isequal(obs.channelIds,ref.channelIds)
    error(['observation and referenceObservation must have identical time, ' ...
        'cell, and channel grids.']);
end
end

% -------------------------------------------------------------------------
function [power,timeOrders,kOrders] = transform_channels( ...
        data,timeWindow,spaceWindow,nFftTime,nFftSpace,removeMean)
if removeMean
    data = data-mean(data,1);
end
windowed = data.*reshape(timeWindow,[],1,1).* ...
    reshape(spaceWindow,1,[],1);
spectrum = ifft(windowed,nFftTime,1)* ...
    (nFftTime/sum(timeWindow));
spectrum = fft(spectrum,nFftSpace,2)/sum(spaceWindow);
timeOrders = (-floor(nFftTime/2):ceil(nFftTime/2)-1).';
kOrders = -floor(nFftSpace/2):ceil(nFftSpace/2)-1;
timeIndices = mod(timeOrders,nFftTime)+1;
kIndices = mod(kOrders,nFftSpace)+1;
spectrum = spectrum(timeIndices,kIndices,:);
power = abs(spectrum).^2;
end

% -------------------------------------------------------------------------
function [folded,foldedOrders] = fold_temporal_power( ...
        raw,timeOrders,foldCount)
foldedNatural = zeros(foldCount,size(raw,2),size(raw,3));
for row = 1:numel(timeOrders)
    target = mod(timeOrders(row),foldCount)+1;
    foldedNatural(target,:,:) = foldedNatural(target,:,:)+raw(row,:,:);
end
foldedOrders = (-floor(foldCount/2):ceil(foldCount/2)-1).';
indices = mod(foldedOrders,foldCount)+1;
folded = foldedNatural(indices,:,:);
end

% -------------------------------------------------------------------------
function [indices,powers] = extract_ridges( ...
        power,support,count,jumpBins,exclusionBins)
[nOmega,nK] = size(power);
indices = nan(count,nK);
powers = nan(count,nK);
available = support;
for ridge = 1:count
    activeColumns = any(available,1);
    transitions = diff([false activeColumns false]);
    starts = find(transitions == 1);
    stops = find(transitions == -1)-1;
    for segment = 1:numel(starts)
        columns = starts(segment):stops(segment);
        [path,pathPower] = ridge_segment( ...
            power(:,columns),available(:,columns),jumpBins);
        indices(ridge,columns) = path;
        powers(ridge,columns) = pathPower;
    end
    for column = 1:nK
        if isnan(indices(ridge,column)), continue; end
        center = indices(ridge,column);
        rows = mod((center-exclusionBins:center+exclusionBins)-1,nOmega)+1;
        available(rows,column) = false;
    end
end
end

% -------------------------------------------------------------------------
function [path,pathPower] = ridge_segment(power,support,jumpBins)
[nOmega,nK] = size(power);
score = -inf(nOmega,nK);
back = zeros(nOmega,nK);
local = log(max(power,realmin));
local(~support) = -inf;
score(:,1) = local(:,1);
for column = 2:nK
    for row = 1:nOmega
        if ~support(row,column), continue; end
        candidates = mod((row-jumpBins:row+jumpBins)-1,nOmega)+1;
        distance = min(abs(candidates-row),nOmega-abs(candidates-row));
        transitionScore = score(candidates,column-1)- ...
            (distance(:).^2)/(2*max(jumpBins,1)^2);
        [best,bestIndex] = max(transitionScore);
        if isfinite(best)
            score(row,column) = local(row,column)+best;
            back(row,column) = candidates(bestIndex);
        end
    end
end
[best,last] = max(score(:,end));
path = nan(1,nK);
pathPower = nan(1,nK);
if ~isfinite(best), return; end
path(end) = last;
pathPower(end) = power(last,end);
for column = nK:-1:2
    last = back(last,column);
    if last == 0, return; end
    path(column-1) = last;
    pathPower(column-1) = power(last,column-1);
end
end

% -------------------------------------------------------------------------
function window = periodic_hann(count)
index = (0:count-1).';
window = 0.5-0.5*cos(2*pi*index/count);
end

% -------------------------------------------------------------------------
function value = read_positive(cfg,name,defaultValue)
value = read_scalar(cfg,name,defaultValue);
if value <= 0, error('fftCfg.%s must be positive.',name); end
end

function value = read_integer(cfg,name,defaultValue,minimum)
value = read_scalar(cfg,name,defaultValue);
if value ~= round(value) || value < minimum
    error('fftCfg.%s must be an integer not smaller than %d.',name,minimum);
end
end

function value = read_logical(cfg,name,defaultValue)
if isfield(cfg,name) && ~isempty(cfg.(name)), value = cfg.(name); else, value = defaultValue; end
if ~islogical(value) || ~isscalar(value)
    error('fftCfg.%s must be a logical scalar.',name);
end
end

function value = read_scalar(cfg,name,defaultValue)
if isfield(cfg,name) && ~isempty(cfg.(name)), value = cfg.(name); else, value = defaultValue; end
if isempty(value) || ~isnumeric(value) || ~isscalar(value) || ...
        ~isreal(value) || ~isfinite(value)
    error('fftCfg.%s must be a finite real scalar.',name);
end
end
