function result = tl_fdtd_gaussian_k_scan(model,scanCfg)
%TL_FDTD_GAUSSIAN_K_SCAN Run finite-chain FDTD for Gaussian packets.
%
%   result = tl_fdtd_gaussian_k_scan(model,scanCfg)
%
% Each kScan entry is the centre of a finite-width complex voltage packet.
% The packet is embedded in a complete Q/Phi finite-chain state using an
% instantaneous local bulk eigenvector at the packet centre, then advanced
% by tl_fdtd1d.  The requested V envelope is imposed exactly even when the
% finite-chain capacitance is spatially nonuniform.  Only fixed
% node-voltage probes are retained for every k; one representative run can
% additionally retain a full spatial voltage/current history.
%
% Required scanCfg fields:
%   kScan, dt, nSteps, pulseIntensityFwhmCells, voltageAmplitude,
%   centerNode, probeNodeIndices
%
% Important optional fields:
%   targetInitialFrequencyHz   selects the positive-frequency CROW branch
%   representativeK           nearest k whose full field is retained
%   recordEvery               default 1
%   precision                 'single' (default) or 'double'
%   boundaryType              'open' (default), 'matched', or 'short'
%   modulationEnabled         default true
%   modulationStart/End       passed to tl_fdtd1d
%   initialTailTolerance      default 1e-4
%   requireNoBoundaryArrival  default true
%   zeroKPropagationDirection default +1 for a zero-frequency acoustic mode
%
% pulseIntensityFwhmCells is the FWHM of |V|^2, not of V.  The returned
% horizontal coordinate remains the source centre k_c; every column has a
% finite circular Bloch-spectrum width and is not an exact bulk eigenvalue.

if nargin ~= 2 || ~isstruct(model) || ~isscalar(model) || ...
        ~isstruct(scanCfg) || ~isscalar(scanCfg)
    error('Use tl_fdtd_gaussian_k_scan(model,scanCfg).');
end
requiredModel = {'kind','cell','series','shunt','resonator','finite', ...
    'modulation','functions','snapshot','derived'};
missingModel = requiredModel(~isfield(model,requiredModel));
if ~isempty(missingModel)
    error('model is missing field(s): %s.',strjoin(missingModel,', '));
end
required = {'kScan','dt','nSteps','pulseIntensityFwhmCells', ...
    'voltageAmplitude','centerNode','probeNodeIndices'};
missing = required(~isfield(scanCfg,required));
if ~isempty(missing)
    error('scanCfg is missing field(s): %s.',strjoin(missing,', '));
end

kScan = validate_k_scan(scanCfg.kScan,model.cell.a);
nK = numel(kScan);
dt = read_positive(scanCfg,'dt',[]);
nSteps = read_integer(scanCfg,'nSteps',[],1);
recordEvery = read_integer(scanCfg,'recordEvery',1,1);
if mod(nSteps,recordEvery) ~= 0
    error('scanCfg.recordEvery must divide scanCfg.nSteps.');
end
fwhmCells = read_positive(scanCfg,'pulseIntensityFwhmCells',[]);
voltageAmplitude = read_positive(scanCfg,'voltageAmplitude',[]);
nNode = model.finite.cellCount;
centerNode = read_integer(scanCfg,'centerNode',[],1);
if centerNode > nNode
    error('scanCfg.centerNode exceeds the finite node count.');
end
probeNodeIndices = validate_indices( ...
    scanCfg.probeNodeIndices,nNode,'scanCfg.probeNodeIndices');
precision = read_text(scanCfg,'precision','single',{'single','double'});
boundaryType = read_text(scanCfg,'boundaryType','open', ...
    {'open','matched','short'});
modulationEnabled = read_logical(scanCfg,'modulationEnabled',true);
modulationStart = read_nonnegative(scanCfg,'modulationStart',0);
modulationEnd = read_positive_or_inf(scanCfg,'modulationEnd',Inf);
initialTailTolerance = read_fraction( ...
    scanCfg,'initialTailTolerance',1e-4);
requireNoBoundaryArrival = read_logical( ...
    scanCfg,'requireNoBoundaryArrival',true);
progressEvery = read_integer( ...
    scanCfg,'progressEvery',max(1,ceil(nK/10)),1);
returnProbeSignals = read_logical(scanCfg,'returnProbeSignals',true);
zeroKPropagationDirection = read_real( ...
    scanCfg,'zeroKPropagationDirection',1);
if ~ismember(zeroKPropagationDirection,[-1 1])
    error('scanCfg.zeroKPropagationDirection must equal -1 or +1.');
end

if isfield(scanCfg,'targetInitialFrequencyHz') && ...
        ~isempty(scanCfg.targetInitialFrequencyHz)
    targetOmega = 2*pi*read_nonnegative( ...
        scanCfg,'targetInitialFrequencyHz',[]);
else
    targetOmega = model.modulation.OmegaRadPerSec/2;
end
if isfield(scanCfg,'representativeK') && ~isempty(scanCfg.representativeK)
    representativeK = read_real(scanCfg,'representativeK',[]);
else
    representativeK = kScan(ceil(nK/2));
end
[~,representativeKIndex] = min(abs(kScan-representativeK));

a = model.cell.a;
xNode = (0:nNode-1)*a;
xCenter = xNode(centerNode);
fwhm = fwhmCells*a;
amplitudeWidth = fwhm/sqrt(2*log(2));
distanceToBoundary = min(xCenter-xNode(1),xNode(end)-xCenter);
tailRadius = amplitudeWidth*sqrt(log(1/initialTailTolerance));
minimumSeriesInductance = model.series.Ls-2*abs(model.series.mutualS);
minimumCapacitance = model.shunt.totalC0-model.shunt.deltaC* ...
    max(model.finite.amplitudeScaleByCell);
velocityScale = a/sqrt(minimumSeriesInductance*minimumCapacitance);
travelDistanceScale = velocityScale*nSteps*dt;
edgeAmplitudeBound = exp(-(distanceToBoundary/amplitudeWidth)^2);
if edgeAmplitudeBound > initialTailTolerance
    error(['The initial voltage-packet amplitude at a chain endpoint exceeds ' ...
        'scanCfg.initialTailTolerance. Increase cellCount or reduce FWHM.']);
end
if requireNoBoundaryArrival && ...
        tailRadius+travelDistanceScale >= distanceToBoundary
    error(['The finite chain is too short for the requested record under ' ...
        'the registered coupling-velocity scale. Increase cellCount, ' ...
        'shorten the analysis, or explicitly set requireNoBoundaryArrival ' ...
        'to false and perform a boundary convergence study.']);
end

if strcmp(boundaryType,'short') && any(probeNodeIndices == [1 nNode])
    error('Voltage probes at short-circuited terminal nodes are not useful.');
end
nBranch = nNode-1;
xBranchStart = (0:nBranch-1)*a;
xBranchMidpoint = xBranchStart+a/2;
nodeEnvelope = exp(-((xNode-xCenter)/amplitudeWidth).^2);

nRecord = nSteps/recordEvery+1;
if returnProbeSignals
    probeSignals = complex(zeros( ...
        nRecord,numel(probeNodeIndices),nK,precision));
else
    probeSignals = complex(zeros(0,0,0,precision));
end
time = [];
representativeField = struct();
initialOmega = zeros(1,nK);
initialModeIndex = zeros(1,nK);
initialVoltageWeight = zeros(1,nK);
initialModeReferenceK = zeros(1,nK);
carrierGroupVelocity = zeros(1,nK);
carrierOmegaLeapfrog = complex(zeros(1,nK));
initialDirectionalPower = zeros(1,nK);
initialVoltageEnvelopeRelativeError = zeros(1,nK);
stabilityAudit = [];
if modulationEnabled && modulationStart == 0
    initialBulkCapacitance = model.functions.capacitanceBulk(0);
    initialFiniteCapacitance = model.functions.capacitanceFinite( ...
        (1:nNode).',0);
else
    initialBulkCapacitance = model.shunt.totalC0;
    initialFiniteCapacitance = ...
        model.shunt.totalC0*ones(nNode,1);
end
initialCenterCapacitance = initialFiniteCapacitance(centerNode);
initialModeCapacitance = initialCenterCapacitance;
capacitanceUniformTolerance = ...
    128*eps(max(abs(initialFiniteCapacitance)));
initialFiniteCapacitanceUniform = ...
    max(abs(initialFiniteCapacitance-initialCenterCapacitance)) <= ...
    capacitanceUniformTolerance;
if ~initialFiniteCapacitanceUniform
    warning('tl_fdtd_gaussian_k_scan:NonuniformInitialCapacitance', ...
        ['The initial finite-chain capacitance is spatially nonuniform. ' ...
        'The voltage Gaussian is imposed exactly, while the companion flux/' ...
        'resonator state is a local centre-cell narrow-band approximation.']);
end
lossyModeTrackingApproximation = model.series.Rs > 0 || ...
    model.shunt.Gp > 0 || ...
    (strcmp(model.kind,'crow') && model.resonator.R0 > 0);
if lossyModeTrackingApproximation
    warning('tl_fdtd_gaussian_k_scan:LossyModalTrackingApproximation', ...
        ['Neighbouring-k tracking uses a right-eigenvector energy-overlap ' ...
        'heuristic for this lossy model.  Near an exceptional point, use ' ...
        'nonzero-k port excitation or a biorthogonal modal analysis.']);
end

for kIndex = 1:nK
    kCenter = kScan(kIndex);
    [modeState,omegaMode,modeIndex,voltageWeight] = ...
        select_initial_mode(model,kCenter,targetOmega, ...
        initialModeCapacitance);
    modeReferenceK = kCenter;
    zeroKTolerance = 64*eps(pi/a);
    useSignedZeroLimit = abs(kCenter) <= zeroKTolerance && ...
        strcmp(model.kind,'sspp');
    if useSignedZeroLimit
        % SSPP 的声学支在 k=0 简并；含损耗时还可能成为纯衰减模。
        % 无条件使用用户指定一侧的 k->0 极限补全阻抗关系，避免
        % eig 的零点排序决定方向。载波空间中心仍严格保持 k_c=0。
        modeReferenceK = zeroKPropagationDirection*1e-6*pi/a;
        [modeState,omegaMode,modeIndex,voltageWeight] = ...
            select_initial_mode(model,modeReferenceK,targetOmega, ...
            initialModeCapacitance);
    end
    omegaLeapfrog = discrete_carrier_omega(omegaMode,dt);
    groupVelocity = estimate_group_velocity(model,kCenter,modeReferenceK, ...
        modeState,omegaMode,initialModeCapacitance, ...
        zeroKPropagationDirection,useSignedZeroLimit,dt);
    initialOmega(kIndex) = omegaMode;
    initialModeIndex(kIndex) = modeIndex;
    initialVoltageWeight(kIndex) = voltageWeight;
    initialModeReferenceK(kIndex) = modeReferenceK;
    carrierGroupVelocity(kIndex) = groupVelocity;
    carrierOmegaLeapfrog(kIndex) = omegaLeapfrog;

    modeState = modeState*(voltageAmplitude/voltageWeight);
    [directionalPower,powerScale] = modal_series_power( ...
        model,modeReferenceK,modeState,initialModeCapacitance);
    initialDirectionalPower(kIndex) = directionalPower;
    if useSignedZeroLimit
        frequencyResolution = ...
            1e-9*model.derived.omegaMaximumEstimate;
        powerResolution = 1e-9*max(powerScale,realmin);
        if abs(real(omegaMode)) <= frequencyResolution || ...
                zeroKPropagationDirection*directionalPower <= powerResolution
            error(['The requested SSPP k_c=0 signed limit has no resolvable ' ...
                'propagating frequency/power direction.  Use a nonzero k_c ' ...
                'or a physical port excitation for this lossy/overdamped model.']);
        end
    end
    nodeCarrier = exp(1i*kCenter*(xNode-xCenter));
    branchCarrier = exp(1i*kCenter*(xBranchStart-xCenter));
    halfStepFactor = exp(1i*omegaLeapfrog*dt/2);
    centerAtMinusHalfStep = xCenter-groupVelocity*dt/2;
    branchEnvelopeHalf = exp(-((xBranchMidpoint- ...
        centerAtMinusHalfStep)/amplitudeWidth).^2);
    nodeEnvelopeHalf = exp(-((xNode-centerAtMinusHalfStep)/ ...
        amplitudeWidth).^2);

    targetVoltage0 = voltageAmplitude*nodeEnvelope.*nodeCarrier;
    q0 = initialFiniteCapacitance.'.*targetVoltage0;
    initialVoltageEnvelopeRelativeError(kIndex) = max(abs( ...
        q0./initialFiniteCapacitance.'-targetVoltage0))/voltageAmplitude;
    phiHalf0 = modeState(2)*halfStepFactor* ...
        branchEnvelopeHalf.*branchCarrier;

    fdtdCfg = struct();
    fdtdCfg.dt = dt;
    fdtdCfg.nSteps = nSteps;
    fdtdCfg.recordEvery = recordEvery;
    fdtdCfg.boundaryType = boundaryType;
    fdtdCfg.modulationEnabled = modulationEnabled;
    fdtdCfg.modulationStart = modulationStart;
    fdtdCfg.modulationEnd = modulationEnd;
    fdtdCfg.precision = precision;
    fdtdCfg.Q0 = q0(:);
    fdtdCfg.PhiHalf0 = phiHalf0(:);
    fdtdCfg.source = struct('type','none');
    if strcmp(model.kind,'crow')
        fdtdCfg.I0Half0 = modeState(3)*halfStepFactor/model.resonator.L0* ...
            (nodeEnvelopeHalf.*nodeCarrier).';
        if ~isinf(model.resonator.Cblock)
            fdtdCfg.Qblock0 = (modeState(4)/modeState(1))*q0(:);
        end
    end
    if ~isempty(stabilityAudit)
        fdtdCfg.stabilityAudit = stabilityAudit;
    end

    if kIndex == representativeKIndex
        fdtdCfg.nodeIndices = 1:nNode;
        fdtdCfg.branchIndices = 1:nBranch;
    else
        fdtdCfg.nodeIndices = probeNodeIndices;
        fdtdCfg.branchIndices = nearest_branches(probeNodeIndices,nBranch);
    end
    field = tl_fdtd1d(model,fdtdCfg);
    if isempty(stabilityAudit)
        stabilityAudit = field.grid.stabilityAudit;
    end
    if isempty(time)
        time = field.node.t(:);
    elseif ~isequal(time,field.node.t(:))
        error('The recorded time grids differ between k scans.');
    end

    if returnProbeSignals
        [found,columns] = ismember(probeNodeIndices,field.node.indices);
        if ~all(found)
            error('Internal error: a requested voltage probe was not recorded.');
        end
        probeSignals(:,:,kIndex) = field.node.V(:,columns);
    end
    if kIndex == representativeKIndex
        representativeField = field;
        representativeField.k = kCenter;
        representativeField.kIndex = kIndex;
    end
    if mod(kIndex,progressEvery) == 0 || kIndex == 1 || kIndex == nK
        fprintf('Finite-chain Gaussian k scan: %d/%d, k_c a/pi=%+.6g\n', ...
            kIndex,nK,kCenter*a/pi);
    end
end

result.k = kScan;
result.kaOverPi = kScan*a/pi;
result.probeSignals = probeSignals;
result.time = time;
result.probeNodeIndices = probeNodeIndices;
result.probePositions = xNode(probeNodeIndices);
result.representativeField = representativeField;
result.representativeKIndex = representativeKIndex;
result.initialOmega = initialOmega;
result.initialModeIndex = initialModeIndex;
result.initialVoltageWeightBeforeScaling = initialVoltageWeight;
result.initialModeReferenceK = initialModeReferenceK;
result.carrierGroupVelocity = carrierGroupVelocity;
result.carrierOmegaLeapfrog = carrierOmegaLeapfrog;
result.initialDirectionalPower = initialDirectionalPower;
result.zeroKPropagationDirection = zeroKPropagationDirection;
result.initialBulkCapacitance = initialBulkCapacitance;
result.initialCenterCapacitance = initialCenterCapacitance;
result.initialModeCapacitance = initialModeCapacitance;
result.initialFiniteCapacitance = initialFiniteCapacitance;
result.initialFiniteCapacitanceUniform = initialFiniteCapacitanceUniform;
result.lossyModeTrackingApproximation = ...
    lossyModeTrackingApproximation;
result.initialVoltageEnvelopeRelativeError = ...
    initialVoltageEnvelopeRelativeError;
result.pulseIntensityFwhm = fwhm;
result.pulseIntensityFwhmCells = fwhmCells;
result.approximateIntensityKFwhm = 4*log(2)/fwhm;
result.initialAmplitudeWidth = amplitudeWidth;
result.initialTailTolerance = initialTailTolerance;
result.initialEdgeAmplitudeBound = edgeAmplitudeBound;
result.distanceToBoundary = distanceToBoundary;
result.couplingVelocityScale = velocityScale;
result.travelDistanceScale = travelDistanceScale;
result.requireNoBoundaryArrival = requireNoBoundaryArrival;
result.recordEvery = recordEvery;
result.recordTimeStep = recordEvery*dt;
result.stabilityAudit = stabilityAudit;
result.observableName = 'node voltage V (electric-field circuit proxy)';
result.modelSnapshot = model.snapshot;
result.scanConfig = scanCfg;
result.interpretation = ['Each column is a finite-chain response weighted ' ...
    'by a finite-width Gaussian voltage source centred at k_c.  For a ' ...
    'spatially nonuniform initial capacitance, its companion state uses a ' ...
    'local centre-cell narrow-band mode approximation.'];
end

% -------------------------------------------------------------------------
function [power,powerScale] = modal_series_power(model,k,state,capacitance)
physicalState = midpoint_gauge(state,k,model.cell.a);
seriesInductance = model.series.Ls+ ...
    2*model.series.mutualS*cos(k*model.cell.a);
voltage = physicalState(1)/capacitance;
current = physicalState(2)/seriesInductance;
power = 0.5*real(voltage*conj(current));
powerScale = 0.5*abs(voltage)*abs(current);
if ~isfinite(power) || ~isfinite(powerScale)
    error('Selected initial mode has non-finite series power.');
end
end

% -------------------------------------------------------------------------
function velocity = estimate_group_velocity( ...
        model,k,referenceK,referenceState,referenceOmega,modeCapacitance, ...
        zeroDirection,useSignedZeroLimit,dt)
a = model.cell.a;
kLimit = pi/a;
step = 1e-4*kLimit;
if useSignedZeroLimit
    kLeft = zeroDirection*0.5*step;
    kRight = zeroDirection*1.5*step;
else
    kLeft = max(-kLimit,k-step);
    kRight = min(kLimit,k+step);
end
if kRight == kLeft
    velocity = 0;
    return;
end
[~,omegaLeft] = select_tracked_mode(model,kLeft,modeCapacitance, ...
    referenceState,referenceOmega,referenceK);
[~,omegaRight] = select_tracked_mode(model,kRight,modeCapacitance, ...
    referenceState,referenceOmega,referenceK);
omegaLeft = discrete_carrier_omega(omegaLeft,dt);
omegaRight = discrete_carrier_omega(omegaRight,dt);
velocity = (real(omegaRight)-real(omegaLeft))/(kRight-kLeft);
if ~isfinite(velocity)
    error('Unable to estimate a finite carrier group velocity.');
end
end

% -------------------------------------------------------------------------
function omegaDiscrete = discrete_carrier_omega(omegaContinuous,dt)
% 对无损交错 leapfrog 振子，sin(omega_d*dt/2)=omega*dt/2。
% 有损扫描仍把它作为窄带初值近似；实际演化始终由完整有限链内核决定。
omegaDiscrete = (2/dt)*asin(omegaContinuous*dt/2);
if ~isfinite(omegaDiscrete)
    error('Unable to construct a finite leapfrog carrier frequency.');
end
end

% -------------------------------------------------------------------------
function [state,omegaSelected,index,voltageWeight] = ...
        select_initial_mode(model,k,targetOmega,modeCapacitance)
[Aconstant,AinverseC] = model.functions.bulkMatrices(k);
A = Aconstant+AinverseC/modeCapacitance;
[vectors,values] = eig(A,'vector');
omega = 1i*values;
voltage = vectors(1,:)/modeCapacitance;
tolerance = 1e-10*max([1;abs(omega)]);
positive = find(real(omega) >= -tolerance & ...
    abs(voltage(:)) > 1e-12*max(abs(voltage)));
if isempty(positive)
    error('No voltage-observable nonnegative-frequency initial mode at k=%.16g.',k);
end
if strcmp(model.kind,'sspp')
    [~,local] = max(real(omega(positive)));
else
    distance = abs(real(omega(positive))-targetOmega);
    distance = distance+0.01*targetOmega* ...
        (1-abs(voltage(positive)).'/max(abs(voltage(positive))));
    [~,local] = min(distance);
end
index = positive(local);
state = vectors(:,index);
omegaSelected = omega(index);
voltageWeight = voltage(index);
if ~isfinite(voltageWeight) || abs(voltageWeight) <= realmin
    error('Selected initial mode has zero or non-finite voltage weight.');
end
end

% -------------------------------------------------------------------------
function [state,omegaSelected,index] = select_tracked_mode( ...
        model,k,modeCapacitance,referenceState,referenceOmega,referenceK)
[Aconstant,AinverseC] = model.functions.bulkMatrices(k);
A = Aconstant+AinverseC/modeCapacitance;
[vectors,values] = eig(A,'vector');
omega = 1i*values;
voltage = vectors(1,:)/modeCapacitance;
tolerance = 1e-10*max([1;abs(omega);abs(referenceOmega)]);
candidates = find(real(omega) >= -tolerance & ...
    abs(voltage(:)) > 1e-12*max(abs(voltage)));
if isempty(candidates)
    error('No trackable voltage-observable mode at k=%.16g.',k);
end

% 用正定电路能量度量比较相邻 k 的右本征向量，避免不同量纲的 Q/Phi
% 分量直接做欧氏内积，也避免 CROW 近交叉处按 targetOmega 重新跳支。
% 体模型的 Phi gauge 锚在支路起点；比较前转到物理支路中点。
weights = state_energy_weights( ...
    model,referenceK,modeCapacitance,numel(referenceState));
referencePhysical = midpoint_gauge(referenceState,referenceK,model.cell.a);
referenceNorm = sqrt(real( ...
    referencePhysical'*(weights.*referencePhysical)));
overlap = zeros(size(candidates));
for candidateIndex = 1:numel(candidates)
    vector = midpoint_gauge( ...
        vectors(:,candidates(candidateIndex)),k,model.cell.a);
    vectorNorm = sqrt(real(vector'*(weights.*vector)));
    overlap(candidateIndex) = abs( ...
        referencePhysical'*(weights.*vector))/ ...
        max(referenceNorm*vectorNorm,realmin);
end
bestOverlap = max(overlap);
nearBest = find(overlap >= bestOverlap-128*eps(max(1,bestOverlap)));
if numel(nearBest) > 1
    [~,nearIndex] = min(abs(omega(candidates(nearBest))-referenceOmega));
    local = nearBest(nearIndex);
else
    local = nearBest;
end
index = candidates(local);
state = vectors(:,index);
omegaSelected = omega(index);
end

% -------------------------------------------------------------------------
function state = midpoint_gauge(state,k,a)
state(2) = state(2)*exp(-1i*k*a/2);
end

% -------------------------------------------------------------------------
function weights = state_energy_weights(model,k,capacitance,stateCount)
seriesInductance = model.series.Ls+ ...
    2*model.series.mutualS*cos(k*model.cell.a);
weights = [1/capacitance;1/seriesInductance];
if strcmp(model.kind,'crow')
    weights(end+1,1) = 1/model.resonator.L0;
    if stateCount == 4
        weights(end+1,1) = 1/model.resonator.Cblock;
    end
end
if numel(weights) ~= stateCount || any(~isfinite(weights)) || ...
        any(weights <= 0)
    error('Unable to construct a positive finite modal-energy metric.');
end
end

% -------------------------------------------------------------------------
function indices = nearest_branches(nodes,nBranch)
indices = unique(max(1,min(nBranch,nodes)));
end

function k = validate_k_scan(value,a)
if ~isnumeric(value) || isempty(value) || ~isvector(value) || ...
        ~isreal(value) || any(~isfinite(value))
    error('scanCfg.kScan must be a nonempty finite real vector.');
end
k = value(:).';
if any(abs(k) > pi/a*(1+100*eps))
    error('scanCfg.kScan must lie in the first spatial Brillouin zone.');
end
end

function indices = validate_indices(value,count,label)
indices = value(:).';
if ~isnumeric(value) || isempty(value) || ~isvector(value) || ...
        ~isreal(value) || any(~isfinite(value)) || ...
        any(value ~= round(value)) || any(value < 1) || any(value > count) || ...
        numel(unique(value)) ~= numel(value)
    error('%s must contain unique valid integer indices.',label);
end
end

function value = read_text(cfg,name,defaultValue,allowed)
if isfield(cfg,name) && ~isempty(cfg.(name)), value = cfg.(name); else, value = defaultValue; end
if isstring(value) && isscalar(value), value = char(value); end
if ~ischar(value) || size(value,1) ~= 1
    error('scanCfg.%s must be a text scalar.',name);
end
match = strcmpi(value,allowed);
if ~any(match)
    error('scanCfg.%s must be one of: %s.',name,strjoin(allowed,', '));
end
value = allowed{find(match,1)};
end

function value = read_logical(cfg,name,defaultValue)
if isfield(cfg,name) && ~isempty(cfg.(name)), value = cfg.(name); else, value = defaultValue; end
if ~islogical(value) || ~isscalar(value)
    error('scanCfg.%s must be a logical scalar.',name);
end
end

function value = read_positive(cfg,name,defaultValue)
value = read_real(cfg,name,defaultValue);
if value <= 0, error('scanCfg.%s must be positive.',name); end
end

function value = read_nonnegative(cfg,name,defaultValue)
value = read_real(cfg,name,defaultValue);
if value < 0, error('scanCfg.%s must be nonnegative.',name); end
end

function value = read_integer(cfg,name,defaultValue,minimum)
value = read_real(cfg,name,defaultValue);
if value ~= round(value) || value < minimum
    error('scanCfg.%s must be an integer not smaller than %d.',name,minimum);
end
end

function value = read_real(cfg,name,defaultValue)
if isfield(cfg,name) && ~isempty(cfg.(name)), value = cfg.(name); else, value = defaultValue; end
if isempty(value) || ~isnumeric(value) || ~isscalar(value) || ...
        ~isreal(value) || ~isfinite(value)
    error('scanCfg.%s must be a finite real scalar.',name);
end
end

function value = read_positive_or_inf(cfg,name,defaultValue)
if isfield(cfg,name) && ~isempty(cfg.(name)), value = cfg.(name); else, value = defaultValue; end
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        isnan(value) || value <= 0
    error('scanCfg.%s must be positive and may be Inf.',name);
end
end

function value = read_fraction(cfg,name,defaultValue)
value = read_real(cfg,name,defaultValue);
if value <= 0 || value >= 1
    error('scanCfg.%s must lie strictly between zero and one.',name);
end
end
