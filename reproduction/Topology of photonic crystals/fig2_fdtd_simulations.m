function fig2_fdtd_simulations(forceRecompute)
%FIG2_FDTD_SIMULATIONS Yee-FDTD reproduction of Lustig et al. Fig. 2.
%
% This implementation uses a real-space D/B Yee grid and leapfrog finite
% differences. It also performs an independent broadband FDTD run and
% obtains the photonic-time-crystal bands from the simulated D(z,t) field:
% spatial FFT -> temporal FFT -> folding into the first temporal Brillouin
% zone. The exact TMM dispersion is plotted only as an audit overlay.
%
%   fig2_fdtd_simulations()       reuse validated cache files
%   fig2_fdtd_simulations(true)   recompute wave packets and FFT bands

if nargin < 1
    forceRecompute = false;
end
validateattributes(forceRecompute,{'logical','numeric'},{'scalar'});
forceRecompute = logical(forceRecompute);

rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(rootDir,'startup_stm.m'));

outputDir = fullfile(fileparts(mfilename('fullpath')),'output');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

% Physical values stated in the published article.
p.c0 = 0.299792458;                 % um/fs
p.eps1 = 3;
p.eps2 = 1;
p.epsBackground = 1;
p.mu = 1;
p.T = 2;                            % fs
p.tStart = 220;                     % fs
p.tEnd = 340;                       % fs
p.nPeriods = round((p.tEnd-p.tStart)/p.T);
p.interfaceTimes = p.tStart:p.T/2:(p.tEnd-p.T/2);
p.nTemporalInterfaces = numel(p.interfaceTimes);
p.publishedFwhm = 45;               % fs, stated in the article

% Reproduction, raster-inference, and numerical choices (not reported FDTD
% settings in the article, which does not publish its grid or time step).
p.branchAuditTime = 400;
p.inferredRasterFwhm = 35;          % fs, fit to the embedded raster
p.fwhm = p.inferredRasterFwhm;
p.sigmaX = p.fwhm/(2*sqrt(log(2)));
p.zInitial = 0;                     % um
p.zMin = -70;                       % um, periodic domain padding
p.zMax = 230;                       % um
p.zPlot = [0 135];                  % um, estimated from the paper raster
p.dxTarget = 0.04;                  % um
p.fdtdDt = p.T/64;                  % fs; 32 steps per half-period
p.recordDt = 0.5;                   % fs

% Broadband FDTD used only for extracting the band structure. Coordinates
% are q=z/c0, so q and t both have units of fs and the solver uses c0=1.
p.band.nPeriods = 24;
p.band.domainPeriods = 40;          % Lq/T; gives Delta(k/k0)=1/40
p.band.Nx = 4096;
p.band.dt = p.T/128;
p.band.recordEvery = 8;             % 16 recorded fields per period
p.band.sourceSigma = p.T/20;        % broadband localized D impulse
p.band.maxKNormalized = 5;
p.band.zeroPaddingFactor = 4;
p.band.displayFloorDb = -45;
p.modelVersion = 'fig2-realspace-db-yee-fdtd-and-field-fft-v3';

validate_parameters(p);

caseBand = struct('name','inband','label','Band propagation', ...
    'lambda',1.4,'clim',[-1.5 0.3],'tStop',465, ...
    'expectedGap',false);
caseGap = struct('name','ingap','label','Gap propagation', ...
    'lambda',0.93,'clim',[-1.5 11],'tStop',520, ...
    'expectedGap',true);

dataBand = run_or_load_wavepacket(caseBand,p,outputDir,forceRecompute);
dataGap = run_or_load_wavepacket(caseGap,p,outputDir,forceRecompute);
bandData = run_or_load_fdtd_bands(p,outputDir,forceRecompute);

fig = figure('Color','w','Position',[40 40 1020 410]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
plot_case(nexttile(tl),dataBand,caseBand,p,'(a)');
plot_case(nexttile(tl),dataGap,caseGap,p,'(b)');

outputFile = fullfile(outputDir,'fig2_fdtd_simulations.png');
save_figure(fig,outputFile,250);
pathTreeFile = fullfile(outputDir,'fig2_temporal_path_tree.png');
plot_path_tree(dataBand.pathTree,p,pathTreeFile);
bandFigureFile = fullfile(outputDir,'fig2_fdtd_fft_band_structure.png');
plot_fdtd_bands(bandData,p,bandFigureFile);

% The paper simultaneously states 60 periods and approximately 20,000
% amplification. Exact Maxwell/TMM propagation exposes their mismatch.
kGap = 2*pi*p.c0/caseGap.lambda;
response20 = temporal_finite_crystal_response(kGap,[p.eps1 p.eps2], ...
    [1 1],[p.T/2 p.T/2],20,p.epsBackground,p.mu,[1;0]);
response60 = temporal_finite_crystal_response(kGap,[p.eps1 p.eps2], ...
    [1 1],[p.T/2 p.T/2],60,p.epsBackground,p.mu,[1;0]);
coherent20 = abs(response20.forward)+abs(response20.backward);
coherent60 = abs(response60.forward)+abs(response60.backward);

audit = struct();
audit.wavepacketMethod = dataBand.method;
audit.bandMethod = bandData.method;
audit.wavepacketCourant = dataBand.fdtd.courant;
audit.bandCourant = bandData.fdtd.courant;
audit.wavepacketTemporalInterfaceCount = ...
    dataBand.fdtd.temporalInterfaceCount;
audit.wavepacketEvolutionUsesSpectralFilter = ...
    dataBand.fdtd.spectralFilterApplied;
audit.coherentTmm20 = coherent20;
audit.coherentTmm60 = coherent60;
audit.lnCoherentTmm20 = log(coherent20);
audit.lnCoherentTmm60 = log(coherent60);
audit.log10CoherentTmm20 = log10(coherent20);
audit.log10CoherentTmm60 = log10(coherent60);
audit.paperClaimedGain = 2e4;
audit.fdtdMaximumGain60 = dataGap.maximumGain;
audit.fdtdToTmm60Ratio = dataGap.maximumGain/coherent60;
audit.clippedPixelFraction = dataGap.clippedPixelFraction;
audit.temporalInterfaceCount = dataBand.interfaceAudit.count;
audit.pathCountsFirstSixInterfaces = dataBand.pathTree.pathCounts;
audit.maximumPathTreeResidual = dataBand.pathTree.maximumRelativeResidual;
audit.bandFftMedianErrorOmegaT = bandData.ridgeAudit.medianErrorOmegaT;
audit.bandFftMaximumErrorOmegaT = bandData.ridgeAudit.maximumErrorOmegaT;
audit.bandFftPassbandPointCount = bandData.ridgeAudit.passbandPointCount;
audit.bandFftNativeResolutionOmegaT = ...
    bandData.ridgeAudit.nativeFrequencyResolutionOmegaT;
audit.bandMinimumInitialSpectrumRelative = ...
    bandData.ridgeAudit.minimumInitialSpectrumRelative;
audit.note = ['The article''s stated 60-period duration and 20000-fold ' ...
    'gain are mutually inconsistent for its stated material parameters.'];

if dataBand.fdtd.spectralFilterApplied || dataGap.fdtd.spectralFilterApplied
    error('The wave-packet evolution must not use a spectral filter.');
end
if dataBand.fdtd.temporalInterfaceCount ~= p.nTemporalInterfaces || ...
        dataGap.fdtd.temporalInterfaceCount ~= p.nTemporalInterfaces
    error('The FDTD run did not process all temporal interfaces.');
end
if ~isfinite(audit.fdtdToTmm60Ratio) || ...
        audit.fdtdToTmm60Ratio < 0.65 || audit.fdtdToTmm60Ratio > 1.35
    error('The FDTD wave-packet gain is inconsistent with the carrier TMM.');
end
if bandData.ridgeAudit.medianErrorOmegaT > 0.10 || ...
        bandData.ridgeAudit.maximumErrorOmegaT > 0.30
    error('The FDTD-FFT band ridges failed the TMM cross-check.');
end

save(fullfile(outputDir,'fig2_audit.mat'), ...
    'p','caseBand','caseGap','audit','-v7.3');

fprintf(['Fig. 2 gain audit at lambda=0.93 um: N=20 gives %.5g, ' ...
    'N=60 gives %.5g, FDTD wave-packet maximum %.5g.\n'], ...
    coherent20,coherent60,dataGap.maximumGain);
fprintf(['FDTD-FFT band audit: median |Delta(Omega T)|=%.4f, ' ...
    'maximum=%.4f over %d pass-band k samples.\n'], ...
    bandData.ridgeAudit.medianErrorOmegaT, ...
    bandData.ridgeAudit.maximumErrorOmegaT, ...
    bandData.ridgeAudit.passbandPointCount);
fprintf('Saved: %s\n',outputFile);
fprintf('Saved: %s\n',pathTreeFile);
fprintf('Saved: %s\n',bandFigureFile);
end

function validate_parameters(p)
if abs(p.nPeriods*p.T-(p.tEnd-p.tStart)) > 10*eps(p.tEnd)
    error('The PTC interval is not an integer number of periods.');
end
if p.nTemporalInterfaces ~= 2*p.nPeriods
    error('Expected two material changes per PTC period.');
end
wavepacketTimes = [p.T/2 p.tStart p.tEnd p.recordDt];
if any(abs(wavepacketTimes/p.fdtdDt- ...
        round(wavepacketTimes/p.fdtdDt)) > 1e-12)
    error('p.fdtdDt must align with interfaces and record times.');
end
bandRecordDt = p.band.recordEvery*p.band.dt;
bandTimes = [p.T/4 p.T/2 p.band.nPeriods*p.T bandRecordDt];
if any(abs(bandTimes/p.band.dt-round(bandTimes/p.band.dt)) > 1e-12)
    error('The broadband FDTD grid does not align with the time crystal.');
end
if mod(round(p.T/bandRecordDt),2) ~= 0
    error('The broadband run needs an even number of records per period.');
end
end

function data = run_or_load_wavepacket(caseDef,p,outputDir,forceRecompute)
dataFile = fullfile(outputDir,['fig2_' caseDef.name '_data.mat']);
if ~forceRecompute && exist(dataFile,'file')
    cached = load(dataFile,'data');
    if isfield(cached,'data')
        candidate = cached.data;
    else
        candidate = struct();
    end
    required = {'modelVersion','method','D','logAmplitude','lambda', ...
        't','z','maximumGain','clippedPixelFraction', ...
        'connectedZoneIsGap','parameters', ...
        'fdtd','interfaceAudit','pathTree'};
    cacheValid = all(isfield(candidate,required));
    if cacheValid
        cacheValid = isstruct(candidate.fdtd) && ...
            all(isfield(candidate.fdtd,{'spectralFilterApplied', ...
            'temporalInterfaceCount','courant'})) && ...
            isstruct(candidate.interfaceAudit) && ...
            isfield(candidate.interfaceAudit,'count') && ...
            isstruct(candidate.pathTree);
    end
    if cacheValid && strcmp(caseDef.name,'inband')
        cacheValid = all(isfield(candidate.pathTree,{'pathCounts', ...
            'maximumRelativeResidual','segments','interfaceTimes'})) && ...
            isstruct(candidate.pathTree.segments) && ...
            all(isfield(candidate.pathTree.segments,{'t0','t1','z0', ...
            'z1','amplitude','direction'}));
    end
    if cacheValid
        cacheValid = strcmp(candidate.modelVersion,p.modelVersion) && ...
            isequaln(candidate.parameters,p) && ...
            strcmp(candidate.method,'real-space D/B Yee FDTD') && ...
            ~isempty(candidate.t) && ~isempty(candidate.z) && ...
            isequal(size(candidate.D), ...
            [numel(candidate.t) numel(candidate.z)]) && ...
            isequal(size(candidate.logAmplitude),size(candidate.D)) && ...
            all(isfinite(candidate.D(:))) && ...
            all(isfinite(candidate.logAmplitude(:))) && ...
            all(isfinite(candidate.t(:))) && all(isfinite(candidate.z(:))) && ...
            isscalar(candidate.maximumGain) && ...
            isfinite(candidate.maximumGain) && ...
            isscalar(candidate.clippedPixelFraction) && ...
            isfinite(candidate.clippedPixelFraction) && ...
            isscalar(candidate.lambda) && isfinite(candidate.lambda) && ...
            abs(candidate.lambda-caseDef.lambda) < 1e-12 && ...
            isscalar(candidate.connectedZoneIsGap) && ...
            candidate.connectedZoneIsGap == caseDef.expectedGap && ...
            abs(candidate.t(end)-caseDef.tStop) < 1e-12 && ...
            isscalar(candidate.fdtd.spectralFilterApplied) && ...
            ~candidate.fdtd.spectralFilterApplied && ...
            isscalar(candidate.fdtd.courant) && ...
            isfinite(candidate.fdtd.courant) && ...
            candidate.fdtd.courant < 1 && ...
            isscalar(candidate.fdtd.temporalInterfaceCount) && ...
            candidate.fdtd.temporalInterfaceCount == ...
            p.nTemporalInterfaces && ...
            isscalar(candidate.interfaceAudit.count) && ...
            candidate.interfaceAudit.count == p.nTemporalInterfaces;
    end
    if cacheValid
        data = candidate;
        fprintf('Loaded validated Yee-FDTD cache: %s\n',dataFile);
        return;
    end
end

fprintf('Running real-space Yee-FDTD %s case (lambda=%.2f um)...\n', ...
    caseDef.name,caseDef.lambda);

domainLengthZ = p.zMax-p.zMin;
Nx = 2^nextpow2(domainLengthZ/p.dxTarget);
domainLengthQ = domainLengthZ/p.c0;
dq = domainLengthQ/Nx;
q = p.zMin/p.c0+(0:Nx-1)*dq;
z = p.c0*q;

orders = [0:Nx/2-1 -Nx/2:-1];
kQ = (2*pi/domainLengthQ)*orders;       % rad/fs in q=z/c0
kPhysical = kQ/p.c0;                   % rad/um
kCarrier = 2*pi/caseDef.lambda;
kCarrierQ = p.c0*kCarrier;

normalizedPosition = (z-p.zInitial)/(p.c0*p.sigmaX);
requestedInitialD = exp(-normalizedPosition.^2).* ...
    exp(1i*kCarrier*(z-p.zInitial));
Dk = fft(requestedInitialD);

% Project the source once. The subsequent evolution remains pure real-space
% FDTD and applies no per-step or per-interface FFT filter.
supportTolerance = 1e-13;
[connectedZone,zoneIsGap,zoneBounds] = carrier_connected_zone( ...
    kQ,kPhysical,kCarrier,p);
if zoneIsGap ~= caseDef.expectedGap
    error('The carrier is not in the expected Floquet band/gap.');
end
physicalSupport = kPhysical > 0 & connectedZone & ...
    abs(Dk) > supportTolerance*max(abs(Dk));
Dk(~physicalSupport) = 0;
initialD = ifft(Dk);
sourceProjectionResidual = norm(initialD-requestedInitialD)/ ...
    norm(requestedInitialD);
if sourceProjectionResidual > 1e-3
    error(['The selected Floquet zone changes the requested source by ' ...
        '%.3g; use a narrower source.'],sourceProjectionResidual);
end

% Initialize H at the Yee location q+dq/2 and time -dt/2. omegaYee is
% the numerical vacuum dispersion, so this is a one-way Yee eigenpacket.
courant = p.fdtdDt/dq;
omegaYee = (2/p.fdtdDt)*asin(courant*sin(kQ*dq/2));
Hhalf0 = ifft(Dk.*exp(1i*(kQ*dq/2+omegaYee*p.fdtdDt/2)));

nSteps = round(caseDef.tStop/p.fdtdDt);
recordEvery = round(p.recordDt/p.fdtdDt);
cfg = struct();
cfg.x = q;
cfg.dt = p.fdtdDt;
cfg.nSteps = nSteps;
cfg.epsFun = @(xq,tq) wavepacket_epsilon_grid(xq,tq,p);
cfg.muFun = @(xq,tq) ones(size(xq));
cfg.E0 = initialD;                   % D=E in the epsilon=1 background
cfg.Hhalf0 = Hhalf0;
cfg.boundary = 'periodic';
cfg.recordEvery = recordEvery;
cfg.storeFields = false;
cfg.storeD = true;
cfg.temporalInterfaces = p.interfaceTimes;
cfg.maxWaveSpeed = 1;               % include the epsilon=1 half-period
cfg.precision = 'double';
cfg.progressBar = false;
out = fdtd1d_db(cfg);

if out.sampledCourant >= 1
    error('The wave-packet FDTD violates the CFL condition.');
end
if ~isempty(out.spectralFilterMask)
    error('No spectral filter is permitted during FDTD evolution.');
end

crop = z >= p.zPlot(1) & z <= p.zPlot(2);
zCrop = z(crop);
Dcrop = single(out.D(:,crop));
t = out.t;
initialPeak = max(abs(Dcrop(1,:)));
logAmplitude = log10(max(abs(Dcrop)/single(initialPeak),single(1e-15)));

[~,idStart] = min(abs(t-p.tStart));
[~,idEnd] = min(abs(t-p.tEnd));
[~,idBranchAudit] = min(abs(t-p.branchAuditTime));
[~,idFinal] = min(abs(t-caseDef.tStop));
[~,peakStartId] = max(abs(Dcrop(idStart,:)));
zPeakAtStart = zCrop(peakStartId);
expectedStart = p.zInitial+p.c0*p.tStart;
if abs(zPeakAtStart-expectedStart) > 2
    error('The FDTD onset ridge is inconsistent with vacuum propagation.');
end

data = struct();
data.modelVersion = p.modelVersion;
data.method = 'real-space D/B Yee FDTD';
data.name = caseDef.name;
data.lambda = caseDef.lambda;
data.kCarrierPhysical = kCarrier;
data.kCarrierNormalizedCoordinate = kCarrierQ;
data.kCarrierNumerical = kCarrierQ; % compatibility alias: q=z/c0 coordinate
data.z = zCrop;
data.t = t(:);
data.D = Dcrop;
data.logAmplitude = logAmplitude;
data.logBase = 10;
data.displayQuantity = 'log10(|D|/max|D(t=0)|)';
data.initialPeak = double(initialPeak);
data.maximumGain = double(max(abs(Dcrop(:)))/initialPeak);
data.zPeakAtStart = zPeakAtStart;
data.expectedZPeakAtStart = expectedStart;
data.dxMicrometre = p.c0*out.dx;
data.recordDt = p.recordDt;
data.spectralSupportCount = nnz(physicalSupport);
data.spectralSupportTolerance = supportTolerance;
data.sourceProjectionResidual = sourceProjectionResidual;
data.connectedZoneIsGap = zoneIsGap;
data.connectedZoneBoundsRadPerMicrometre = zoneBounds;
data.clippedPixelFraction = mean(logAmplitude(:) >= caseDef.clim(2));
data.parameters = p;
data.fdtd = struct( ...
    'scheme','D/B Yee leapfrog finite differences', ...
    'Nx',Nx, ...
    'dxNormalizedTime',out.dx, ...
    'dxMicrometre',p.c0*out.dx, ...
    'dtFs',out.dt, ...
    'courant',out.sampledCourant, ...
    'pointsPerCarrierWavelength',caseDef.lambda/(p.c0*out.dx), ...
    'boundary',out.boundary, ...
    'temporalInterfaceCount',out.processedTemporalInterfaceCount, ...
    'configuredTemporalInterfaceCount',numel(out.temporalInterfaces), ...
    'temporalInterfaces',out.temporalInterfaces, ...
    'maximumInterfaceGridError', ...
    max(abs(out.temporalInterfaces-p.interfaceTimes)), ...
    'spectralFilterApplied',~isempty(out.spectralFilterMask), ...
    'sourceProjectionOnly',true);
data.interfaceAudit = struct( ...
    'count',out.processedTemporalInterfaceCount, ...
    'times',out.temporalInterfaces, ...
    'maximumGridAlignmentError', ...
    max(abs(out.temporalInterfaces-p.interfaceTimes)), ...
    'continuityVariables','D and B', ...
    'crossingRule','centred average of D/eps before and after interface');

if strcmp(caseDef.name,'inband')
    groupSpeed = floquet_group_speed(kCarrierQ,p);
    expectedEnd = expectedStart+[-1 1]*p.c0*groupSpeed* ...
        (p.tEnd-p.tStart);
    expectedFinal = reshape(expectedEnd(:)+ ...
        [-1 1]*p.c0*(caseDef.tStop-p.tEnd),1,[]);
    expectedBranchAudit = sort(reshape(expectedEnd(:)+ ...
        [-1 1]*p.c0*(p.branchAuditTime-p.tEnd),1,[]));
else
    groupSpeed = 0;
    expectedEnd = [expectedStart expectedStart];
    expectedFinal = expectedStart+[-1 1]*p.c0*(caseDef.tStop-p.tEnd);
    expectedBranchAudit = [];
end
data.floquetGroupSpeedOverC = groupSpeed;
data.expectedPeaksAtEnd = expectedEnd;
data.measuredPeaksAtEnd = window_peaks( ...
    zCrop,abs(Dcrop(idEnd,:)),expectedEnd,12);
visible = expectedFinal >= p.zPlot(1) & expectedFinal <= p.zPlot(2);
data.expectedVisiblePeaksAtFinal = expectedFinal(visible);
data.measuredVisiblePeaksAtFinal = window_peaks( ...
    zCrop,abs(Dcrop(idFinal,:)),expectedFinal(visible),12);
if strcmp(caseDef.name,'inband')
    data.expectedFourBranchesAtAudit = expectedBranchAudit;
    data.measuredFourBranchesAtAudit = window_peaks( ...
        zCrop,abs(Dcrop(idBranchAudit,:)),expectedBranchAudit,8);
    data.pathTree = build_path_tree(kCarrierQ,p,6);
else
    data.expectedFourBranchesAtAudit = [];
    data.measuredFourBranchesAtAudit = [];
    data.pathTree = struct();
end

trajectoryFailure = any(abs(data.measuredPeaksAtEnd(:)- ...
    expectedEnd(:)) > 8);
if ~strcmp(caseDef.name,'inband')
    trajectoryFailure = trajectoryFailure || ...
        any(abs(data.measuredVisiblePeaksAtFinal(:)- ...
        expectedFinal(visible).') > 8);
end
if trajectoryFailure
    error('The Yee-FDTD ridges failed the Floquet trajectory audit.');
end
if strcmp(caseDef.name,'inband') && ...
        (any(~isfinite(data.measuredFourBranchesAtAudit)) || ...
        any(abs(data.measuredFourBranchesAtAudit- ...
        expectedBranchAudit) > 8))
    error('The four FDTD output branches failed their position audit.');
end
if strcmp(caseDef.name,'inband') && ...
        (data.pathTree.maximumRelativeResidual > 1e-11 || ...
        ~isequal(data.pathTree.pathCounts,2.^(1:6)))
    error('The temporal-interface path-tree audit failed.');
end
if strcmp(caseDef.name,'inband') && data.maximumGain > 20
    error('The in-band FDTD run developed unexpected gain.');
end

save(dataFile,'data','-v7.3');
fprintf(['  Nx=%d, dz=%.5f um, dt=%.4f fs, CFL=%.4f, ' ...
    'points/lambda=%.1f, source projection=%.3e, max gain=%.5g.\n'], ...
    Nx,data.fdtd.dxMicrometre,data.fdtd.dtFs,data.fdtd.courant, ...
    data.fdtd.pointsPerCarrierWavelength,sourceProjectionResidual, ...
    data.maximumGain);
fprintf('  Processed %d interfaces; evolution spectral filter: %d.\n', ...
    data.fdtd.temporalInterfaceCount,data.fdtd.spectralFilterApplied);
fprintf('  Saved: %s\n',dataFile);
end

function bandData = run_or_load_fdtd_bands(p,outputDir,forceRecompute)
dataFile = fullfile(outputDir,'fig2_fdtd_fft_band_data.mat');
if ~forceRecompute && exist(dataFile,'file')
    cached = load(dataFile,'bandData');
    if isfield(cached,'bandData')
        candidate = cached.bandData;
    else
        candidate = struct();
    end
    required = {'modelVersion','method','D','q','z','t', ...
        'kNormalized','omegaT','spectralDb','ridgePositiveOmegaT', ...
        'ridgeNegativeOmegaT', ...
        'exactTmmPositiveOmegaT','inBand','ridgeAudit','parameters','fdtd'};
    cacheValid = all(isfield(candidate,required));
    if cacheValid
        cacheValid = isstruct(candidate.fdtd) && ...
            all(isfield(candidate.fdtd,{'courant','spectralFilterApplied'})) && ...
            isstruct(candidate.ridgeAudit) && ...
            all(isfield(candidate.ridgeAudit,{'medianErrorOmegaT', ...
            'maximumErrorOmegaT','passbandPointCount', ...
            'minimumInitialSpectrumRelative', ...
            'nativeFrequencyResolutionOmegaT'}));
    end
    if cacheValid
        cacheValid = strcmp(candidate.modelVersion,p.modelVersion) && ...
            isequaln(candidate.parameters,p) && ...
            strcmp(candidate.method, ...
            'broadband real-space Yee-FDTD plus 2-D field FFT') && ...
            ~isempty(candidate.t) && ~isempty(candidate.q) && ...
            ~isempty(candidate.kNormalized) && ~isempty(candidate.omegaT) && ...
            numel(candidate.z) == numel(candidate.q) && ...
            isequal(size(candidate.D), ...
            [numel(candidate.t) numel(candidate.q)]) && ...
            isequal(size(candidate.spectralDb), ...
            [numel(candidate.omegaT) numel(candidate.kNormalized)]) && ...
            numel(candidate.ridgePositiveOmegaT) == ...
            numel(candidate.kNormalized) && ...
            numel(candidate.ridgeNegativeOmegaT) == ...
            numel(candidate.kNormalized) && ...
            numel(candidate.exactTmmPositiveOmegaT) == ...
            numel(candidate.kNormalized) && ...
            numel(candidate.inBand) == numel(candidate.kNormalized) && ...
            all(isfinite(candidate.D(:))) && ...
            all(isfinite(candidate.q(:))) && all(isfinite(candidate.z(:))) && ...
            all(isfinite(candidate.t(:))) && ...
            all(isfinite(candidate.kNormalized(:))) && ...
            all(isfinite(candidate.omegaT(:))) && ...
            all(isfinite(candidate.spectralDb(:))) && ...
            all(isfinite(candidate.ridgePositiveOmegaT(:))) && ...
            all(isfinite(candidate.ridgeNegativeOmegaT(:))) && ...
            all(isfinite(candidate.exactTmmPositiveOmegaT(:))) && ...
            isscalar(candidate.fdtd.courant) && ...
            isfinite(candidate.fdtd.courant) && ...
            candidate.fdtd.courant < 1 && ...
            isscalar(candidate.fdtd.spectralFilterApplied) && ...
            ~candidate.fdtd.spectralFilterApplied && ...
            isscalar(candidate.ridgeAudit.medianErrorOmegaT) && ...
            isfinite(candidate.ridgeAudit.medianErrorOmegaT) && ...
            isscalar(candidate.ridgeAudit.maximumErrorOmegaT) && ...
            isfinite(candidate.ridgeAudit.maximumErrorOmegaT) && ...
            isscalar(candidate.ridgeAudit.passbandPointCount) && ...
            isfinite(candidate.ridgeAudit.passbandPointCount) && ...
            isscalar(candidate.ridgeAudit.minimumInitialSpectrumRelative) && ...
            candidate.ridgeAudit.minimumInitialSpectrumRelative > 1e-3 && ...
            isscalar(candidate.ridgeAudit.nativeFrequencyResolutionOmegaT) && ...
            isfinite(candidate.ridgeAudit.nativeFrequencyResolutionOmegaT);
    end
    if cacheValid
        bandData = candidate;
        fprintf('Loaded validated FDTD-FFT band cache: %s\n',dataFile);
        return;
    end
end

fprintf('Running broadband Yee-FDTD for field-FFT band extraction...\n');
Lq = p.band.domainPeriods*p.T;
Nx = p.band.Nx;
dq = Lq/Nx;
q = (-Nx/2:Nx/2-1)*dq;
tStop = p.band.nPeriods*p.T;
nSteps = round(tStop/p.band.dt);
% Use a low-epsilon-centred symmetric unit cell. The initial and final
% states are then sampled away from a material jump at the same phase.
interfaceTimes = (p.T/4):(p.T/2):(tStop-p.T/4);

initialD = exp(-(q/p.band.sourceSigma).^2);
cfg = struct();
cfg.x = q;
cfg.dt = p.band.dt;
cfg.nSteps = nSteps;
cfg.epsFun = @(xq,tq) infinite_ptc_epsilon_grid(xq,tq,p);
cfg.muFun = @(xq,tq) ones(size(xq));
cfg.E0 = initialD/p.eps2;           % make the initialized D equal initialD
cfg.Hhalf0 = zeros(size(q));        % excites both +/- Floquet branches
cfg.boundary = 'periodic';
cfg.recordEvery = p.band.recordEvery;
cfg.storeFields = false;
cfg.storeD = true;
cfg.temporalInterfaces = interfaceTimes;
cfg.maxWaveSpeed = 1;               % default sampling can miss eps=1
cfg.precision = 'double';
cfg.progressBar = false;
out = fdtd1d_db(cfg);

if out.sampledCourant >= 1 || ~isempty(out.spectralFilterMask)
    error('The broadband FDTD CFL/filter audit failed.');
end

[kNormalized,omegaT,spectralDb,ridgePositive,ridgeNegative,exactOmegaT, ...
    inBand,ridgeAudit] = field_fft_bands(out.D,out.t,q,p);

bandData = struct();
bandData.modelVersion = p.modelVersion;
bandData.method = 'broadband real-space Yee-FDTD plus 2-D field FFT';
bandData.fieldQuantity = 'electric displacement D(q,t)';
bandData.fftPipeline = ['FFT along q, FFT along t, incoherent sum of ' ...
    'Floquet harmonics folded modulo 2*pi/T'];
bandData.q = q;
bandData.z = p.c0*q;
bandData.t = out.t;
bandData.D = out.D;
bandData.kNormalized = kNormalized;
bandData.omegaT = omegaT;
bandData.spectralDb = spectralDb;
bandData.ridgePositiveOmegaT = ridgePositive;
bandData.ridgeNegativeOmegaT = ridgeNegative;
bandData.exactTmmPositiveOmegaT = exactOmegaT;
bandData.inBand = inBand;
bandData.ridgeAudit = ridgeAudit;
bandData.initialD = single(initialD);
bandData.maximumFieldGain = double(max(abs(out.D(:)))/max(abs(out.D(1,:))));
bandData.parameters = p;
bandData.fdtd = struct( ...
    'scheme','D/B Yee leapfrog finite differences', ...
    'Nx',Nx, ...
    'dxNormalizedTime',out.dx, ...
    'dxMicrometre',p.c0*out.dx, ...
    'dtFs',out.dt, ...
    'recordDtFs',p.band.recordEvery*out.dt, ...
    'courant',out.sampledCourant, ...
    'boundary',out.boundary, ...
    'periodCount',p.band.nPeriods, ...
    'recordsPerPeriod',round(p.T/(p.band.recordEvery*out.dt)), ...
    'temporalInterfaceCount',out.processedTemporalInterfaceCount, ...
    'configuredTemporalInterfaceCount',numel(out.temporalInterfaces), ...
    'spectralFilterApplied',~isempty(out.spectralFilterMask));

save(dataFile,'bandData','-v7.3');
fprintf(['  Nx=%d, dq=%.6f fs, dt=%.6f fs, CFL=%.3f, ' ...
    '%d periods, %d k samples.\n'],Nx,dq,out.dt,out.sampledCourant, ...
    p.band.nPeriods,numel(kNormalized));
fprintf('  Saved raw FDTD field and FFT bands: %s\n',dataFile);
end

function [kNormalized,omegaT,spectralDb,ridgePositive,ridgeNegative, ...
    exactOmegaT, ...
    inBand,audit] = field_fft_bands(Dhistory,t,q,p)
% Convert a broadband FDTD field history into Floquet quasifrequency bands.
if size(Dhistory,1) ~= numel(t) || size(Dhistory,2) ~= numel(q)
    error('The FDTD field history dimensions are inconsistent.');
end
dtRecord = mean(diff(t));
Danalysis = Dhistory(1:end-1,:);
nTime = size(Danalysis,1);
Nx = size(Danalysis,2);
dq = q(2)-q(1);

orders = 0:Nx/2-1;
kNormalizedAll = orders*p.T/(Nx*dq);
selected = find(kNormalizedAll <= p.band.maxKNormalized+1e-12);
kNormalized = kNormalizedAll(selected);

initialSpectrum = abs(fft(Dhistory(1,:),[],2));
initialSpectrumRelative = initialSpectrum(selected)/max(initialSpectrum);
minimumInitialSpectrumRelative = min(initialSpectrumRelative);
if minimumInitialSpectrumRelative <= 1e-3
    error(['The broadband source does not reliably excite every selected ' ...
        'k column; minimum relative source amplitude is %.3g.'], ...
        minimumInitialSpectrumRelative);
end

% Actual spatial FFT of the complete FDTD field history.
DkAll = fft(Danalysis,[],2);
Dk = DkAll(:,selected);
clear DkAll Danalysis;

% Per-k display normalization prevents unstable gaps from hiding pass bands.
columnNorm = sqrt(sum(abs(Dk).^2,1));
columnNorm(columnNorm == 0) = 1;
Dk = Dk./columnNorm;
window = 0.5-0.5*cos(2*pi*(0:nTime-1)'/nTime);
if isa(Dk,'single'), window = single(window); end
Dk = Dk.*window;

nFft = p.band.zeroPaddingFactor*nTime;
complexSpectrum = fftshift(ifft(Dk,nFft,1),1);
powerFull = abs(complexSpectrum).^2;
clear complexSpectrum Dk;

% Floquet harmonics differ by 2*pi in omega*T and align to FFT bins.
nFoldFloat = nFft*dtRecord/p.T;
nFold = round(nFoldFloat);
if abs(nFold-nFoldFloat) > 1e-10 || mod(nFold,2) ~= 0 || ...
        mod(nFft,nFold) ~= 0
    error('FFT sampling is not commensurate with the modulation period.');
end
nHarmonicBlocks = nFft/nFold;
powerFold = sum(reshape(powerFull, ...
    nFold,nHarmonicBlocks,numel(selected)),2);
powerFold = reshape(powerFold,nFold,numel(selected));
powerFold = fftshift(powerFold,1);
omegaT = (-nFold/2:nFold/2-1)'*(2*pi/nFold);

columnMaximum = max(powerFold,[],1);
columnMaximum(columnMaximum == 0) = 1;
powerFold = powerFold./columnMaximum;
powerFloor = 1e-8;
if isa(powerFold,'single'), powerFloor = single(powerFloor); end
spectralDb = 10*log10(max(powerFold,powerFloor));

positiveRows = find(omegaT >= 0);
[~,localPeak] = max(powerFold(positiveRows,:),[],1);
ridgePositive = omegaT(positiveRows(localPeak)).';
negativeRows = find(omegaT <= 0);
[~,localPeak] = max(powerFold(negativeRows,:),[],1);
ridgeNegative = omegaT(negativeRows(localPeak)).';

kQ = kNormalized*(2*pi/p.T);
bands = temporal_crystal_bands(kQ,[p.eps1 p.eps2], ...
    [p.mu p.mu],[p.T/2 p.T/2]);
halfTrace = real(bands.halfTrace);
inBand = abs(halfTrace) <= 1+5e-11;
exactOmegaT = acos(min(1,max(-1,halfTrace)));

positiveErrors = abs(ridgePositive(inBand)-exactOmegaT(inBand));
negativeErrors = abs(ridgeNegative(inBand)+exactOmegaT(inBand));
errors = [positiveErrors negativeErrors];
sortedErrors = sort(errors);
audit = struct();
audit.passbandPointCount = nnz(inBand);
audit.medianErrorOmegaT = median(errors);
audit.maximumErrorOmegaT = max(errors);
audit.positiveBranchMedianErrorOmegaT = median(positiveErrors);
audit.negativeBranchMedianErrorOmegaT = median(negativeErrors);
audit.percentile95ErrorOmegaT = sortedErrors( ...
    max(1,ceil(0.95*numel(sortedErrors))));
nativePeriodCount = nTime*dtRecord/p.T;
audit.nativeFrequencyResolutionOmegaT = 2*pi/nativePeriodCount;
audit.zeroPaddedGridSpacingOmegaT = 2*pi/nFold;
audit.spatialBinSpacingKOverK0 = p.T/(Nx*dq);
audit.zeroPaddingFactor = p.band.zeroPaddingFactor;
audit.minimumInitialSpectrumRelative = minimumInitialSpectrumRelative;
end

function epsGrid = wavepacket_epsilon_grid(x,t,p)
epsGrid = ptc_epsilon_scalar(t,p)*ones(size(x));
end

function epsGrid = infinite_ptc_epsilon_grid(x,t,p)
tau = mod(t,p.T);
if tau < p.T/4 || tau >= 3*p.T/4
    epsValue = p.eps2;
else
    epsValue = p.eps1;
end
epsGrid = epsValue*ones(size(x));
end

function epsValue = ptc_epsilon_scalar(t,p)
if t < p.tStart || t >= p.tEnd
    epsValue = p.epsBackground;
elseif mod(t-p.tStart,p.T) < p.T/2
    epsValue = p.eps1;
else
    epsValue = p.eps2;
end
end

function plot_fdtd_bands(data,p,outputFile)
fig = figure('Color','w','Position',[70 70 820 520]);
ax = axes(fig);
imagesc(ax,data.kNormalized,data.omegaT,data.spectralDb);
axis(ax,'xy');
colormap(ax,parula(256));
clim(ax,[p.band.displayFloorDb 0]);
hold(ax,'on');

[gapStarts,gapEnds] = logical_segments(~data.inBand);
for g = 1:numel(gapStarts)
    x1 = data.kNormalized(gapStarts(g));
    x2 = data.kNormalized(gapEnds(g));
    patch(ax,[x1 x2 x2 x1],[-pi -pi pi pi],[0.65 0.65 0.65], ...
        'EdgeColor','none','FaceAlpha',0.12,'HandleVisibility','off');
end

[bandStarts,bandEnds] = logical_segments(data.inBand);
for b = 1:numel(bandStarts)
    ids = bandStarts(b):bandEnds(b);
    plot(ax,data.kNormalized(ids),data.exactTmmPositiveOmegaT(ids), ...
        'k-','LineWidth',1.0,'HandleVisibility','off');
    plot(ax,data.kNormalized(ids),-data.exactTmmPositiveOmegaT(ids), ...
        'k-','LineWidth',1.0,'HandleVisibility','off');
end
hTheory = plot(ax,nan,nan,'k-','LineWidth',1.0, ...
    'DisplayName','TMM audit');
passbandIds = find(data.inBand);
markerIds = passbandIds(1:2:end);
hFDTD = scatter(ax,data.kNormalized(markerIds), ...
    data.ridgePositiveOmegaT(markerIds),10,'w','filled', ...
    'DisplayName','FDTD-FFT pass-band ridge');
scatter(ax,data.kNormalized(markerIds), ...
    data.ridgeNegativeOmegaT(markerIds),10,'w','filled', ...
    'HandleVisibility','off');

xlabel(ax,'k/k_0,  k_0=2\pi/(Tc)');
ylabel(ax,'\Omega_F T');
xlim(ax,[0 p.band.maxKNormalized]);
ylim(ax,[-pi pi]);
yticks(ax,[-pi 0 pi]);
yticklabels(ax,{'-\pi','0','\pi'});
title(ax,{'Photonic-time-crystal bands from the Yee-FDTD field', ...
    sprintf(['2-D FFT, Floquet folded; median ridge error ' ...
    '|\Delta\Omega T|=%.3f'],data.ridgeAudit.medianErrorOmegaT)}, ...
    'FontWeight','normal');
cb = colorbar(ax);
cb.Label.String = 'Folded FFT power [dB; normalized at each k]';
legend(ax,[hFDTD hTheory],'Location','southoutside', ...
    'Orientation','horizontal','Box','off');
box(ax,'on');
set(ax,'FontSize',10,'Layer','top');
save_figure(fig,outputFile,250);
end

function plot_case(ax,data,caseDef,p,panelLabel)
imagesc(ax,data.z,data.t,data.logAmplitude);
axis(ax,'xy');
colormap(ax,jet(256));
clim(ax,caseDef.clim);
xlim(ax,p.zPlot);
ylim(ax,[0 caseDef.tStop]);
xticks(ax,[0 50 100]);
yticks(ax,[0 200 400]);
xlabel(ax,'z [\mum]');
ylabel(ax,'t [fs]');
title(ax,caseDef.label,'FontWeight','normal');
text(ax,-0.12,1.08,panelLabel,'Units','normalized', ...
    'FontSize',13,'Clipping','off');
cb = colorbar(ax);
cb.Label.String = data.displayQuantity;
cb.Label.Interpreter = 'none';
if strcmp(caseDef.name,'inband')
    cb.Ticks = [-1.5 -1 -0.5 0];
else
    cb.Ticks = [0 5 10];
end
pbaspect(ax,[1.24 1 1]);
set(ax,'FontSize',9,'Layer','top');
box(ax,'on');
end

function tree = build_path_tree(kQ,p,nLevels)
% Analytic diagnostic only: enumerate coherent temporal-interface paths.
validateattributes(nLevels,{'numeric'},{'scalar','integer','positive'});
if nLevels > numel(p.interfaceTimes)
    error('The path-tree depth exceeds the number of PTC interfaces.');
end

directions = 1;
amplitudes = 1;
positions = p.zInitial+p.c0*p.tStart;
directState = [1;0];
epsBefore = p.epsBackground;

segmentT0 = [];
segmentT1 = [];
segmentZ0 = [];
segmentZ1 = [];
segmentAmplitude = [];
segmentDirection = [];
segmentLevel = [];
pathCounts = zeros(1,nLevels);
relativeResidual = zeros(1,nLevels);

for level = 1:nLevels
    tInterface = p.interfaceTimes(level);
    epsAfter = ptc_epsilon_scalar(tInterface+1e-9,p);
    J = temporal_interface_matrix_jump( ...
        epsBefore,p.mu,epsAfter,p.mu,'DB');

    nParents = numel(directions);
    childDirections = zeros(1,2*nParents);
    childAmplitudes = complex(zeros(1,2*nParents));
    childPositions = zeros(1,2*nParents);
    childId = 0;
    for parent = 1:nParents
        inputChannel = 1+(directions(parent) < 0);
        for outputChannel = 1:2
            childId = childId+1;
            childDirections(childId) = 3-2*outputChannel;
            childAmplitudes(childId) = ...
                J(outputChannel,inputChannel)*amplitudes(parent);
            childPositions(childId) = positions(parent);
        end
    end

    dtLayer = p.T/2;
    omega = kQ/sqrt(epsAfter*p.mu);
    endPositions = childPositions+childDirections* ...
        (p.c0/sqrt(epsAfter*p.mu))*dtLayer;
    endAmplitudes = childAmplitudes.* ...
        exp(-1i*childDirections*omega*dtLayer);

    nChildren = numel(childDirections);
    segmentT0 = [segmentT0 repmat(tInterface,1,nChildren)]; %#ok<AGROW>
    segmentT1 = [segmentT1 repmat(tInterface+dtLayer,1,nChildren)]; %#ok<AGROW>
    segmentZ0 = [segmentZ0 childPositions]; %#ok<AGROW>
    segmentZ1 = [segmentZ1 endPositions]; %#ok<AGROW>
    segmentAmplitude = [segmentAmplitude abs(childAmplitudes)]; %#ok<AGROW>
    segmentDirection = [segmentDirection childDirections]; %#ok<AGROW>
    segmentLevel = [segmentLevel repmat(level,1,nChildren)]; %#ok<AGROW>

    phase = [exp(-1i*omega*dtLayer);exp(1i*omega*dtLayer)];
    directState = phase.*(J*directState);
    coherentSum = [sum(endAmplitudes(childDirections > 0)); ...
        sum(endAmplitudes(childDirections < 0))];
    relativeResidual(level) = norm(coherentSum-directState)/ ...
        max(norm(directState),eps);

    directions = childDirections;
    amplitudes = endAmplitudes;
    positions = endPositions;
    epsBefore = epsAfter;
    pathCounts(level) = nChildren;
end

tree = struct();
tree.depth = nLevels;
tree.interfaceTimes = p.interfaceTimes(1:nLevels);
tree.pathCounts = pathCounts;
tree.relativeResidual = relativeResidual;
tree.maximumRelativeResidual = max(relativeResidual);
tree.maximumSpreadMicrometre = max(positions)-min(positions);
tree.pulseSpatialFwhmMicrometre = p.c0*p.fwhm;
tree.leafDirections = directions;
tree.leafAmplitudes = amplitudes;
tree.leafPositions = positions;
tree.segments = struct('t0',segmentT0,'t1',segmentT1, ...
    'z0',segmentZ0,'z1',segmentZ1,'amplitude',segmentAmplitude, ...
    'direction',segmentDirection,'level',segmentLevel);
end

function plot_path_tree(tree,p,outputFile)
fig = figure('Color','w','Position',[60 60 690 430]);
ax = axes(fig);
hold(ax,'on');
segments = tree.segments;
maximumAmplitude = max(segments.amplitude);
forwardColor = [0.85 0.20 0.15];
reverseColor = [0.10 0.35 0.80];
for j = 1:numel(segments.t0)
    if segments.direction(j) > 0
        color = forwardColor;
    else
        color = reverseColor;
    end
    relativeAmplitude = segments.amplitude(j)/maximumAmplitude;
    lineWidth = 0.35+1.8*sqrt(relativeAmplitude);
    plot(ax,[segments.z0(j) segments.z1(j)], ...
        [segments.t0(j) segments.t1(j)],'Color',color, ...
        'LineWidth',lineWidth,'HandleVisibility','off');
end
for tInterface = tree.interfaceTimes
    yline(ax,tInterface,':','Color',[0.65 0.65 0.65], ...
        'LineWidth',0.6,'HandleVisibility','off');
end
hForward = plot(ax,nan,nan,'Color',forwardColor,'LineWidth',1.8, ...
    'DisplayName','time-refracted (+)');
hReverse = plot(ax,nan,nan,'Color',reverseColor,'LineWidth',1.8, ...
    'DisplayName','time-reflected (-)');

allZ = [segments.z0 segments.z1];
zPad = max(0.5,0.12*(max(allZ)-min(allZ)));
xlim(ax,[min(allZ)-zPad max(allZ)+zPad]);
ylim(ax,[p.tStart p.tStart+tree.depth*p.T/2]);
xlabel(ax,'z [\mum]');
ylabel(ax,'t [fs]');
title(ax,{sprintf('First %d temporal interfaces: coherent path expansion', ...
    tree.depth),sprintf(['2, 4, 8, 16, 32, 64 paths; spread %.2f um ' ...
    '< pulse FWHM %.2f um'],tree.maximumSpreadMicrometre, ...
    tree.pulseSpatialFwhmMicrometre)},'FontWeight','normal', ...
    'FontSize',10);
legend(ax,[hForward hReverse],'Location','eastoutside');
box(ax,'on');
grid(ax,'on');
set(ax,'FontSize',9,'Layer','top');
save_figure(fig,outputFile,250);
close(fig);
end

function [zoneMask,isGap,bounds] = carrier_connected_zone( ...
    kQ,kPhysical,kCarrier,p)
q1 = sqrt(p.eps1/p.mu);
q2 = sqrt(p.eps2/p.mu);
theta1 = kQ*(p.T/2)/sqrt(p.eps1*p.mu);
theta2 = kQ*(p.T/2)/sqrt(p.eps2*p.mu);
halfTrace = cos(theta1).*cos(theta2)- ...
    0.5*(q1/q2+q2/q1).*sin(theta1).*sin(theta2);
gapMask = abs(real(halfTrace)) > 1+1e-12;

positiveIds = find(kPhysical > 0);
[~,localCarrierId] = min(abs(kPhysical(positiveIds)-kCarrier));
carrierId = positiveIds(localCarrierId);
isGap = gapMask(carrierId);
sameZone = (gapMask == isGap) & kPhysical > 0;

firstId = carrierId;
while firstId > positiveIds(1) && sameZone(firstId-1)
    firstId = firstId-1;
end
lastId = carrierId;
while lastId < positiveIds(end) && sameZone(lastId+1)
    lastId = lastId+1;
end
zoneMask = false(size(kPhysical));
zoneMask(firstId:lastId) = true;
bounds = kPhysical([firstId lastId]);
end

function speed = floquet_group_speed(kQ,p)
dk = 1e-5*max(1,abs(kQ));
b = temporal_crystal_bands([kQ-dk kQ+dk], ...
    [p.eps1 p.eps2],[p.mu p.mu],[p.T/2 p.T/2]);
h = min(1,max(-1,real(b.halfTrace)));
omega = acos(h)/p.T;
speed = abs(diff(omega)/(2*dk));
end

function positions = window_peaks(z,amplitude,expected,halfWidth)
positions = nan(size(expected));
for j = 1:numel(expected)
    ids = find(abs(z-expected(j)) <= halfWidth);
    if isempty(ids), continue; end
    [~,localId] = max(amplitude(ids));
    positions(j) = z(ids(localId));
end
end

function [starts,ends] = logical_segments(mask)
d = diff([false mask(:).' false]);
starts = find(d == 1);
ends = find(d == -1)-1;
end

function save_figure(fig,outputFile,resolution)
try
    exportgraphics(fig,outputFile,'Resolution',resolution);
catch
    print(fig,outputFile,'-dpng',sprintf('-r%d',resolution));
end
end
