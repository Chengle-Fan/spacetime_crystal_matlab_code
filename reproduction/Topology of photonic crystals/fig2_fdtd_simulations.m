function fig2_fdtd_simulations(forceRecompute)
%FIG2_FDTD_SIMULATIONS Reproduce Fig. 2 of Lustig et al. (2018).
%
% The medium is homogeneous in space, hence every spatial Fourier
% component evolves independently.  This implementation propagates the
% exact Maxwell state [D;B] in k space and reconstructs D(z,t) by FFT.  It
% is the spectral (zero-spatial-discretization-error) counterpart of the
% paper's FDTD calculation and, unlike a long real-space FDTD run, cannot
% seed unstable momentum gaps with roundoff noise outside the pulse band.

if nargin < 1
    forceRecompute = false;
end
validateattributes(forceRecompute,{'logical','numeric'},{'scalar'});
forceRecompute = logical(forceRecompute);

rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(rootDir,'startup_stm.m'));

outputDir = fullfile(fileparts(mfilename('fullpath')),'output');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

% Parameters stated in the published article.
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
p.branchAuditTime = 400;            % all four band outputs are in view
p.publishedFwhm = 45;               % fs, value stated in the article text
p.inferredRasterFwhm = 35;          % fs, fit to the embedded figure in log10
p.fwhm = p.inferredRasterFwhm;      % reproduce the plotted pulse width
p.sigmaX = p.fwhm/(2*sqrt(log(2)));
p.zInitial = 0;                     % um
p.zMin = -70;                       % periodic FFT domain is padded
p.zMax = 230;
p.zPlot = [0 135];                  % axes calibrated from the paper raster
p.dxTarget = 0.04;                  % um; rounded to an FFT-friendly grid
p.recordDt = 0.5;                   % fs
p.modelVersion = 'fig2-exact-k-connected-zone-v22';

if abs(p.nPeriods*p.T-(p.tEnd-p.tStart)) > 10*eps(p.tEnd)
    error('The PTC interval is not an integer number of periods.');
end
if p.nTemporalInterfaces ~= 2*p.nPeriods
    error('Expected two temporal interfaces per PTC period.');
end
gridTimes = [p.T/2 p.tStart p.tEnd];
if any(abs(gridTimes/p.recordDt-round(gridTimes/p.recordDt)) > 1e-12)
    error('recordDt must resolve every temporal interface exactly.');
end

% The white brace is a display inset reused at the same normalized height
% in the two published panels, even though their time axes differ.  These
% rasterMarkerTimes affect only that annotation; the simulated PTC remains
% exactly p.tStart..p.tEnd in both cases.
caseBand = struct('name','inband','label','Band propagation', ...
    'lambda',1.4,'clim',[-1.5 0.3],'tStop',465, ...
    'expectedGap',false, ...
    'rasterMarkerTimes',[p.tStart p.tEnd]);
caseGap = struct('name','ingap','label','Gap propagation', ...
    'lambda',0.93,'clim',[-1.5 11],'tStop',520, ...
    'expectedGap',true, ...
    'rasterMarkerTimes',520/465*[p.tStart p.tEnd]);

dataBand = run_or_load_case(caseBand,p,outputDir,forceRecompute);
dataGap = run_or_load_case(caseGap,p,outputDir,forceRecompute);

% The text simultaneously specifies 60 periods and approximately 20,000
% amplification.  Exact Maxwell propagation makes the inconsistency
% explicit: the latter occurs after about 20 periods, not after 60.
kGap = 2*pi*p.c0/caseGap.lambda;
response20 = temporal_finite_crystal_response(kGap,[p.eps1 p.eps2], ...
    [1 1],[p.T/2 p.T/2],20,p.epsBackground,p.mu,[1;0]);
response60 = temporal_finite_crystal_response(kGap,[p.eps1 p.eps2], ...
    [1 1],[p.T/2 p.T/2],60,p.epsBackground,p.mu,[1;0]);
coherent20 = abs(response20.forward)+abs(response20.backward);
coherent60 = abs(response60.forward)+abs(response60.backward);

fprintf(['Fig. 2 parameter audit at lambda=0.93 um: N=20 gives ' ...
    '%.5g (log10 %.3f), while the published N=60 gives %.5g ' ...
    '(log10 %.3f).  The image keeps the stated 60-period dynamics and ' ...
    'the raster-calibrated base-10 color limits; unclipped fields are ' ...
    'saved.\n'],coherent20,log10(coherent20), ...
    coherent60,log10(coherent60));

fig = figure('Color','w','Position',[40 40 1020 410]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
plot_case(nexttile(tl),dataBand,caseBand,p,'(a)');
plot_case(nexttile(tl),dataGap,caseGap,p,'(b)');

outputFile = fullfile(outputDir,'fig2_fdtd_simulations.png');
exportgraphics(fig,outputFile,'Resolution',250);
pathTreeFile = fullfile(outputDir,'fig2_temporal_path_tree.png');
plot_path_tree(dataBand.pathTree,p,pathTreeFile);

audit = struct();
audit.coherentTmm20 = coherent20;
audit.coherentTmm60 = coherent60;
audit.lnCoherentTmm20 = log(coherent20);
audit.lnCoherentTmm60 = log(coherent60);
audit.log10CoherentTmm20 = log10(coherent20);
audit.log10CoherentTmm60 = log10(coherent60);
audit.paperClaimedGain = 2e4;
audit.exactSpectralMaximumGain60 = dataGap.maximumGain;
audit.spectralToTmm60Ratio = dataGap.maximumGain/coherent60;
audit.clippedPixelFraction = dataGap.clippedPixelFraction;
audit.temporalInterfaceCount = dataBand.interfaceAudit.count;
audit.maximumInterfaceMixingResidual = ...
    dataBand.interfaceAudit.maximumRelativeResidual;
audit.pathCountsFirstSixInterfaces = dataBand.pathTree.pathCounts;
audit.maximumPathTreeResidual = dataBand.pathTree.maximumRelativeResidual;
audit.note = ['The stated 60-period duration and 20000-fold gain are ' ...
    'mutually inconsistent for the stated material parameters.'];
if ~isfinite(audit.spectralToTmm60Ratio) || ...
        audit.spectralToTmm60Ratio < 0.65 || ...
        audit.spectralToTmm60Ratio > 1.35
    error('The exact wavepacket peak is inconsistent with the carrier TMM.');
end
save(fullfile(outputDir,'fig2_audit.mat'),'p','caseBand','caseGap','audit');
fprintf('Saved: %s\n',outputFile);
fprintf('Saved: %s\n',pathTreeFile);
end

function data = run_or_load_case(caseDef,p,outputDir,forceRecompute)
dataFile = fullfile(outputDir,['fig2_' caseDef.name '_data.mat']);
if ~forceRecompute && exist(dataFile,'file')
    cached = load(dataFile,'data');
    required = {'modelVersion','D','logAmplitude','logBase','lambda', ...
        't','z','dxMicrometre','connectedZoneIsGap', ...
        'sourceProjectionResidual','parameters'};
    if isfield(cached,'data')
        candidate = cached.data;
    else
        candidate = struct();
    end
    cacheValid = all(isfield(candidate,required));
    if cacheValid
        cacheValid = strcmp(candidate.modelVersion,p.modelVersion) && ...
            isequaln(candidate.parameters,p) && ...
            all(isfinite(candidate.D(:))) && ...
            all(isfinite(candidate.logAmplitude(:))) && ...
            candidate.logBase == 10 && ...
            abs(candidate.lambda-caseDef.lambda) < 1e-12 && ...
            candidate.connectedZoneIsGap == caseDef.expectedGap && ...
            abs(candidate.t(end)-caseDef.tStop) < 1e-12 && ...
            abs(candidate.z(1)-p.zPlot(1)) <= ...
            2*candidate.dxMicrometre && ...
            abs(candidate.z(end)-p.zPlot(2)) <= ...
            2*candidate.dxMicrometre;
    end
    if cacheValid
        data = candidate;
        fprintf('Loaded validated Fig. 2 cache: %s\n',dataFile);
        return;
    end
end

fprintf('Running exact-k Fig. 2 %s case (lambda=%.2f um)...\n', ...
    caseDef.name,caseDef.lambda);

domainLength = p.zMax-p.zMin;
Nx = 2^nextpow2(domainLength/p.dxTarget);
dxZ = domainLength/Nx;
z = p.zMin+(0:Nx-1)*dxZ;

if mod(Nx,2) ~= 0
    error('The exact-k reconstruction requires an even FFT grid.');
end
orders = [0:Nx/2-1 -Nx/2:-1];
kPhysical = (2*pi/domainLength)*orders;                 % rad/um
kNumerical = p.c0*kPhysical;                            % rad/fs

kCarrier = 2*pi/caseDef.lambda;
normalizedPosition = (z-p.zInitial)/(p.c0*p.sigmaX);
requestedInitialD = exp(-normalizedPosition.^2).* ...
    exp(1i*kCarrier*(z-p.zInitial));
Dk = fft(requestedInitialD);

% Analytic initial signal: retain the positive-k part of the connected
% Floquet band (or gap) that contains the carrier.  The paper explicitly
% states that the pulse momentum components reside in the selected band or
% gap.  Removing the Gaussian's tiny tails in neighbouring zones is
% essential here: after 60 periods an initially 1e-6 tail in another gap
% can become brighter than the intended in-band pulse.  This mask is
% applied once to the source, never again at temporal interfaces.
supportTolerance = 1e-13;
[connectedZone,zoneIsGap,zoneBounds] = carrier_connected_zone( ...
    kNumerical,kPhysical,kCarrier,p);
if zoneIsGap ~= caseDef.expectedGap
    error('The carrier is not in the expected Floquet band/gap.');
end
physicalSupport = kPhysical > 0 & connectedZone & abs(Dk) > ...
    supportTolerance*max(abs(Dk));
Dk(~physicalSupport) = 0;
initialD = ifft(Dk);
sourceProjectionResidual = norm(initialD-requestedInitialD)/ ...
    norm(requestedInitialD);
if sourceProjectionResidual > 1e-3
    error(['The selected connected Floquet zone changes the requested ' ...
        'input pulse by %.3g; use a spectrally narrower source.'], ...
        sourceProjectionResidual);
end
Bk = Dk;                              % free-space forward wave: B=D

t = (0:round(caseDef.tStop/p.recordDt))*p.recordDt;
if abs(t(end)-caseDef.tStop) > 10*eps(caseDef.tStop)
    error('The panel stop time must align with p.recordDt.');
end
nTime = numel(t);
crop = z >= p.zPlot(1) & z <= p.zPlot(2);
zCrop = z(crop);
Dcrop = complex(zeros(nTime,nnz(crop),'single'));
Dcrop(1,:) = single(initialD(crop));

stateD = Dk;
stateB = Bk;
epsPrevious = p.epsBackground;
interfaceCount = 0;
interfaceTimesObserved = nan(1,p.nTemporalInterfaces);
interfaceEpsBefore = nan(1,p.nTemporalInterfaces);
interfaceEpsAfter = nan(1,p.nTemporalInterfaces);
interfaceTau = nan(1,p.nTemporalInterfaces);
interfaceRho = nan(1,p.nTemporalInterfaces);
interfaceResidual = nan(1,p.nTemporalInterfaces);
for it = 2:nTime
    tMid = 0.5*(t(it-1)+t(it));
    epsStep = ptc_epsilon_scalar(tMid,p);
    dt = t(it)-t(it-1);
    if epsStep ~= epsPrevious
        interfaceCount = interfaceCount+1;
        if interfaceCount > p.nTemporalInterfaces
            error('More temporal interfaces were encountered than expected.');
        end

        % Make the otherwise implicit split explicit.  Both columns of J
        % contribute: each incoming +/- frequency channel generates a
        % refracted and a time-reflected outgoing channel.  Histories with
        % the same outgoing channel must be added coherently.
        state = [stateD;stateB];
        before = temporal_db_to_directional( ...
            state,epsPrevious,p.mu);
        afterFromContinuity = temporal_db_to_directional( ...
            state,epsStep,p.mu);
        [J,tau,rho] = temporal_interface_matrix_jump( ...
            epsPrevious,p.mu,epsStep,p.mu,'DB');
        afterFromSplitting = J*before;
        residual = norm(afterFromContinuity(:)- ...
            afterFromSplitting(:))/max(norm(afterFromContinuity(:)),eps);

        interfaceTimesObserved(interfaceCount) = t(it-1);
        interfaceEpsBefore(interfaceCount) = epsPrevious;
        interfaceEpsAfter(interfaceCount) = epsStep;
        interfaceTau(interfaceCount) = tau;
        interfaceRho(interfaceCount) = rho;
        interfaceResidual(interfaceCount) = residual;
        epsPrevious = epsStep;
    end

    % D and B stay continuous at the interface.  Propagation in the new
    % segment then advances all coherently recombined histories exactly.
    [stateD,stateB] = propagate_spectrum( ...
        stateD,stateB,kNumerical,epsStep,p.mu,dt,physicalSupport);
    field = ifft(stateD);
    Dcrop(it,:) = single(field(crop));
end
if interfaceCount ~= p.nTemporalInterfaces
    error('Encountered %d temporal interfaces; expected %d.', ...
        interfaceCount,p.nTemporalInterfaces);
end
if max(interfaceResidual) > 5e-13 || any(abs(interfaceRho) < 1e-12)
    error('The explicit temporal refraction/reflection audit failed.');
end

initialPeak = max(abs(Dcrop(1,:)));
% The caption says only "log scale".  Both the four branch colours and
% the unsaturated gap-growth slope in the embedded raster identify base
% ten.  Natural log would hide the weakest fourth band output entirely.
logAmplitude = log10(max(abs(Dcrop)/initialPeak,1e-15));

[~,idStart] = min(abs(t-p.tStart));
[~,idEnd] = min(abs(t-p.tEnd));
[~,idBranchAudit] = min(abs(t-p.branchAuditTime));
[~,idFinal] = min(abs(t-caseDef.tStop));
[~,peakStartId] = max(abs(Dcrop(idStart,:)));
zPeakAtStart = zCrop(peakStartId);
expectedStart = p.zInitial+p.c0*p.tStart;
if abs(zPeakAtStart-expectedStart) > 2
    error('Exact-k onset ridge is inconsistent with free-space propagation.');
end

data = struct();
data.modelVersion = p.modelVersion;
data.method = 'exact independent-k Maxwell propagation and FFT';
data.name = caseDef.name;
data.lambda = caseDef.lambda;
data.kCarrierPhysical = kCarrier;
data.kCarrierNumerical = p.c0*kCarrier;
data.z = zCrop;
data.t = t(:);
data.D = Dcrop;
data.logAmplitude = logAmplitude;
data.logBase = 10;
data.displayQuantity = 'log10(|D|/max|D(t=0)|)';
data.initialPeak = initialPeak;
data.maximumGain = max(abs(Dcrop(:)))/initialPeak;
data.zPeakAtStart = zPeakAtStart;
data.expectedZPeakAtStart = expectedStart;
data.dxMicrometre = dxZ;
data.recordDt = p.recordDt;
data.spectralSupportCount = nnz(physicalSupport);
data.spectralSupportTolerance = supportTolerance;
data.sourceProjectionResidual = sourceProjectionResidual;
data.connectedZoneIsGap = zoneIsGap;
data.connectedZoneBoundsRadPerMicrometre = zoneBounds;
data.clippedPixelFraction = mean(logAmplitude(:) >= caseDef.clim(2));
data.parameters = p;
data.interfaceAudit = struct( ...
    'count',interfaceCount, ...
    'times',interfaceTimesObserved(1:interfaceCount), ...
    'epsBefore',interfaceEpsBefore(1:interfaceCount), ...
    'epsAfter',interfaceEpsAfter(1:interfaceCount), ...
    'tau',interfaceTau(1:interfaceCount), ...
    'rho',interfaceRho(1:interfaceCount), ...
    'relativeResidual',interfaceResidual(1:interfaceCount), ...
    'maximumRelativeResidual',max(interfaceResidual));

if strcmp(caseDef.name,'inband')
    groupSpeed = floquet_group_speed(p.c0*kCarrier,p);
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
    data.pathTree = build_path_tree(p.c0*kCarrier,p,6);
else
    data.expectedFourBranchesAtAudit = [];
    data.measuredFourBranchesAtAudit = [];
    data.pathTree = struct();
end

trajectoryFailure = any(abs(data.measuredPeaksAtEnd(:)- ...
    expectedEnd(:)) > 8);
% At the paper's band-panel endpoint the two inner outputs overlap and
% interfere, so they are not two independently locatable final peaks.  The
% separated four-branch audit at 400 fs below is the meaningful check.
if ~strcmp(caseDef.name,'inband')
    trajectoryFailure = trajectoryFailure || ...
        any(abs(data.measuredVisiblePeaksAtFinal(:)- ...
        expectedFinal(visible).') > 8);
end
if trajectoryFailure
    error('The reconstructed ridges fail the Floquet trajectory audit.');
end
if strcmp(caseDef.name,'inband') && ...
        (any(~isfinite(data.measuredFourBranchesAtAudit)) || ...
        any(abs(data.measuredFourBranchesAtAudit- ...
        expectedBranchAudit) > 8))
    error('The four refracted/reflected output branches failed audit.');
end
if strcmp(caseDef.name,'inband') && ...
        (data.pathTree.maximumRelativeResidual > 1e-11 || ...
        ~isequal(data.pathTree.pathCounts,2.^(1:6)))
    error('The repeated temporal-interface path-tree audit failed.');
end
if strcmp(caseDef.name,'inband') && data.maximumGain > 20
    error('The in-band exact-k run developed unexpected gain.');
end

save(dataFile,'data','-v7.3');
fprintf(['  Nx=%d, active k=%d, dz=%.5f um, onset z=%.3f um, ' ...
    'source projection=%.3e, maximum |D|/|D0|=%.5g, ' ...
    'clipped pixels=%.3f%%.\n'], ...
    Nx,data.spectralSupportCount,dxZ,zPeakAtStart, ...
    sourceProjectionResidual,data.maximumGain, ...
    100*data.clippedPixelFraction);
fprintf(['  %d explicit interface mixes, max directional-basis residual ' ...
    '%.3e.\n'],interfaceCount,max(interfaceResidual));
if strcmp(caseDef.name,'inband')
    fprintf('  Four branches at t=%.0f fs: expected [%s], measured [%s] um.\n', ...
        p.branchAuditTime,num2str(expectedBranchAudit,'%.2f '), ...
        num2str(data.measuredFourBranchesAtAudit,'%.2f '));
    fprintf(['  First six interfaces: path counts [%s], max coherent-' ...
        'sum residual %.3e, path spread %.3f um (pulse FWHM %.3f um).\n'], ...
        num2str(data.pathTree.pathCounts,'%d '), ...
        data.pathTree.maximumRelativeResidual, ...
        data.pathTree.maximumSpreadMicrometre, ...
        data.pathTree.pulseSpatialFwhmMicrometre);
end
fprintf('  Saved: %s\n',dataFile);
end

function [Dnew,Bnew] = propagate_spectrum( ...
    D,B,k,epsValue,muValue,dt,support)
rootRatio = sqrt(epsValue/muValue);
phase = k*dt/sqrt(epsValue*muValue);
c = cos(phase);
s = sin(phase);
Dnew = zeros(size(D));
Bnew = zeros(size(B));
Dnew(support) = c(support).*D(support) - ...
    1i*rootRatio*s(support).*B(support);
Bnew(support) = -1i*s(support)/rootRatio.*D(support) + ...
    c(support).*B(support);
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

function tree = build_path_tree(kNumerical,p,nLevels)
% Enumerate histories without treating them as incoherent intensities.
% This diagnostic is deliberately limited to a few interfaces: after m
% jumps there are 2^m histories, although only two coherent +/- channels.
validateattributes(nLevels,{'numeric'},{'scalar','integer','positive'});
if nLevels > numel(p.interfaceTimes)
    error('The path-tree depth exceeds the number of PTC interfaces.');
end

directions = 1;                     % +1: refracted, -1: time reversed
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

    % Propagate all children to the next half-period boundary.  Histories
    % retain their individual centres and complex phases for the audit.
    dtLayer = p.T/2;
    omega = kNumerical/sqrt(epsAfter*p.mu);
    endPositions = childPositions + childDirections* ...
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
% Auxiliary diagnostic; the paper panel remains the coherent total field.
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
    tree.depth),sprintf(['2, 4, 8, 16, 32, 64 paths; spread %.2f μm ' ...
    '< pulse FWHM %.2f μm'],tree.maximumSpreadMicrometre, ...
    tree.pulseSpatialFwhmMicrometre)},'FontWeight','normal', ...
    'FontSize',10);
legend(ax,[hForward hReverse],'Location','eastoutside');
box(ax,'on');
grid(ax,'on');
set(ax,'FontSize',9,'Layer','top');
exportgraphics(fig,outputFile,'Resolution',250);
close(fig);
end

function plot_case(ax,data,caseDef,p,panelLabel)
imagesc(ax,data.z,data.t,data.logAmplitude);
axis(ax,'xy');
colormap(ax,jet(256));
clim(ax,caseDef.clim);
hold(ax,'on');
zMarker = p.zPlot(2)-7;
plot(ax,[zMarker zMarker],caseDef.rasterMarkerTimes,'w-','LineWidth',1.2);
plot(ax,zMarker,caseDef.rasterMarkerTimes(1),'w_', ...
    'MarkerSize',8,'LineWidth',1.2);
plot(ax,zMarker,caseDef.rasterMarkerTimes(2),'w_', ...
    'MarkerSize',8,'LineWidth',1.2);
text(ax,zMarker-3,mean(caseDef.rasterMarkerTimes),'\epsilon(t)', ...
    'Color','w','HorizontalAlignment','right','FontWeight','bold');
text(ax,zMarker,caseDef.rasterMarkerTimes(2)+0.045*caseDef.tStop,'T', ...
    'Color','w','HorizontalAlignment','center','FontWeight','bold');
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
if strcmp(caseDef.name,'inband')
    cb.Ticks = [-1.5 -1 -0.5 0];
else
    cb.Ticks = [0 5 10];
end
pbaspect(ax,[1.24 1 1]);
set(ax,'FontSize',9,'Layer','top');
box(ax,'on');
end

function [zoneMask,isGap,bounds] = carrier_connected_zone( ...
    kNumerical,kPhysical,kCarrier,p)
% Return the one connected pass-band/gap interval containing the carrier.
q1 = sqrt(p.eps1/p.mu);
q2 = sqrt(p.eps2/p.mu);
theta1 = kNumerical*(p.T/2)/sqrt(p.eps1*p.mu);
theta2 = kNumerical*(p.T/2)/sqrt(p.eps2*p.mu);
halfTrace = cos(theta1).*cos(theta2) - ...
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

function speed = floquet_group_speed(k,p)
dk = 1e-5*max(1,abs(k));
b = temporal_crystal_bands([k-dk k+dk], ...
    [p.eps1 p.eps2],[1 1],[p.T/2 p.T/2]);
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
