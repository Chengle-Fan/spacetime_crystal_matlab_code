function result = fdtd_gaussian_k_scan(cfg)
%FDTD_GAUSSIAN_K_SCAN  逐中心 k 运行有限空间高斯波包 FDTD 并记录 E 探针。
% 本函数只负责可复用的“构造 E/H 初值 -> fdtd1d -> 探针/代表场记录”数据
% 流。有限样品几何、材料函数、边界安全距离和 FFT 窗由调用入口显式定义。

if nargin ~= 1 || ~isstruct(cfg) || ~isscalar(cfg)
    error('必须以标量结构体 cfg 调用 fdtd_gaussian_k_scan。');
end
required = {'x','dt','nSteps','epsFun','muFun','kScan','xCenter', ...
    'pulseIntensityFwhm','fieldAmplitude','initialEpsilon', ...
    'muRelative','probeIndices'};
missing = required(~isfield(cfg,required));
if ~isempty(missing)
    error('cfg 缺少字段：%s。',strjoin(missing,', '));
end

%% ===================== 网格、激发与记录配置检查 =====================

if ~isnumeric(cfg.x) || ~isvector(cfg.x) || numel(cfg.x) < 5 || ...
        ~isreal(cfg.x) || any(~isfinite(cfg.x))
    error('cfg.x 必须是至少含 5 点的有限实向量。');
end
x = cfg.x(:).';
dxVector = diff(x);
if any(dxVector <= 0)
    error('cfg.x 必须严格递增。');
end
dx = mean(dxVector);
if max(abs(dxVector-dx)) > 1e-10*max(1,abs(dx))
    error('cfg.x 必须是均匀网格。');
end
Nx = numel(x);

dt = read_positive_scalar(cfg.dt,'cfg.dt');
nSteps = read_integer(cfg.nSteps,'cfg.nSteps',0);
if ~isa(cfg.epsFun,'function_handle') || ~isa(cfg.muFun,'function_handle')
    error('cfg.epsFun 和 cfg.muFun 必须是函数句柄。');
end
if ~isnumeric(cfg.kScan) || ~isvector(cfg.kScan) || isempty(cfg.kScan) || ...
        ~isreal(cfg.kScan) || any(~isfinite(cfg.kScan))
    error('cfg.kScan 必须是非空有限实向量。');
end
kScan = cfg.kScan(:).';
nK = numel(kScan);

xCenter = read_finite_real_scalar(cfg.xCenter,'cfg.xCenter');
if xCenter <= x(1) || xCenter >= x(end)
    error('cfg.xCenter 必须严格位于空间网格内部。');
end
pulseIntensityFwhm = read_positive_scalar( ...
    cfg.pulseIntensityFwhm,'cfg.pulseIntensityFwhm');
fieldAmplitude = read_positive_scalar(cfg.fieldAmplitude,'cfg.fieldAmplitude');
initialEpsilon = read_positive_scalar( ...
    cfg.initialEpsilon,'cfg.initialEpsilon');
muRelative = read_positive_scalar(cfg.muRelative,'cfg.muRelative');

probeIndices = validate_indices(cfg.probeIndices,Nx,'cfg.probeIndices');
if isfield(cfg,'representativeSpatialIndices') && ...
        ~isempty(cfg.representativeSpatialIndices)
    representativeSpatialIndices = validate_indices( ...
        cfg.representativeSpatialIndices,Nx, ...
        'cfg.representativeSpatialIndices');
else
    representativeSpatialIndices = probeIndices;
end
% 代表性全场记录中始终包含探针，便于从同一次 FDTD 数据交叉核查。
representativeSpatialIndices = unique( ...
    [representativeSpatialIndices probeIndices]);

if isfield(cfg,'representativeK') && ~isempty(cfg.representativeK)
    representativeK = read_finite_real_scalar( ...
        cfg.representativeK,'cfg.representativeK');
else
    representativeK = kScan(ceil(nK/2));
end
[~,representativeKIndex] = min(abs(kScan-representativeK));

recordEvery = read_optional_integer(cfg,'recordEvery',1,1);
recordPrecision = read_optional_choice(cfg,'recordPrecision','single', ...
    {'double','single'});
if isfield(cfg,'zeroKPropagationDirection') && ...
        ~isempty(cfg.zeroKPropagationDirection)
    zeroKPropagationDirection = read_finite_real_scalar( ...
        cfg.zeroKPropagationDirection,'cfg.zeroKPropagationDirection');
    if ~ismember(zeroKPropagationDirection,[-1 1])
        error('cfg.zeroKPropagationDirection 必须等于 -1 或 +1。');
    end
else
    zeroKPropagationDirection = 1;
end
progressEvery = read_optional_integer(cfg,'progressEvery',max(1,ceil(nK/10)),1);
if isfield(cfg,'progressLabel') && ~isempty(cfg.progressLabel)
    progressLabel = cfg.progressLabel;
    if isstring(progressLabel) && isscalar(progressLabel)
        progressLabel = char(progressLabel);
    elseif ~ischar(progressLabel) || size(progressLabel,1) ~= 1
        error('cfg.progressLabel 必须是文本标量。');
    end
else
    progressLabel = 'FDTD k 扫描';
end
if isfield(cfg,'progressKScale') && ~isempty(cfg.progressKScale)
    progressKScale = read_positive_scalar( ...
        cfg.progressKScale,'cfg.progressKScale');
else
    progressKScale = 1;
end
if isfield(cfg,'progressKLabel') && ~isempty(cfg.progressKLabel)
    progressKLabel = cfg.progressKLabel;
    if isstring(progressKLabel) && isscalar(progressKLabel)
        progressKLabel = char(progressKLabel);
    elseif ~ischar(progressKLabel) || size(progressKLabel,1) ~= 1
        error('cfg.progressKLabel 必须是文本标量。');
    end
else
    progressKLabel = 'k';
end

if mod(nSteps,recordEvery) ~= 0
    error('cfg.recordEvery 必须整除 cfg.nSteps，保证所有扫描具有相同终止记录。');
end

amplitudeGaussianWidth = pulseIntensityFwhm/sqrt(2*log(2));
if any(abs(kScan)+6/amplitudeGaussianWidth >= pi/dx)
    error('某些高斯波包的主要空间谱接近 Nyquist 边界；请细化网格或增大 FWHM。');
end

%% ===================== 公共 fdtd1d 配置 =====================

fdtdBase = struct();
fdtdBase.x = x;
fdtdBase.dt = dt;
fdtdBase.nSteps = nSteps;
fdtdBase.epsFun = cfg.epsFun;
fdtdBase.muFun = cfg.muFun;
fdtdBase.recordEvery = recordEvery;
fdtdBase.recordField = 'E';
fdtdBase.recordPrecision = recordPrecision;

copyFields = {'boundaryType','spongeCells','spongeStrength', ...
    'temporalInterfaces','certifiedMinimumEpsilon','certifiedMinimumMu'};
for fieldIndex = 1:numel(copyFields)
    fieldName = copyFields{fieldIndex};
    if isfield(cfg,fieldName) && ~isempty(cfg.(fieldName))
        fdtdBase.(fieldName) = cfg.(fieldName);
    end
end

%% ===================== 逐 k 高斯初值与完整空间推进 =====================

nRecordedTimes = nSteps/recordEvery+1;
probeSignals = complex(zeros( ...
    nRecordedTimes,numel(probeIndices),nK,recordPrecision));
probeTime = [];
representativeField = struct();
xH = x(1:end-1)+dx/2;
initialSpeed = 1/sqrt(initialEpsilon*muRelative);
electricEnvelope = exp(-2*log(2)* ...
    ((x-xCenter)/pulseIntensityFwhm).^2);
electricCarrierCoordinate = x-xCenter;
magneticCarrierCoordinate = xH-xCenter;
propagationDirections = sign(kScan);
propagationDirections(propagationDirections == 0) = ...
    zeroKPropagationDirection;
carrierOmegaYee = zeros(1,nK);
carrierGroupSpeedYee = zeros(1,nK);

for kIndex = 1:nK
    kCenter = kScan(kIndex);
    E0 = fieldAmplitude*electricEnvelope.* ...
        exp(1i*kCenter*electricCarrierCoordinate);

    yeeArgument = initialSpeed*dt/dx*sin(abs(kCenter)*dx/2);
    if abs(yeeArgument) >= 1
        error('扫描索引 %d 的载波 k=%.16g 不满足 Yee 离散色散条件。', ...
            kIndex,kCenter);
    end
    omegaYee = (2/dt)*asin(yeeArgument);
    propagationSign = propagationDirections(kIndex);
    % 包络半步位移使用与载波一致的 Yee 群速度，而不是连续
    % 介质相速度。这仍是有限带宽包络的窄带初值，不是精确本征模。
    groupSpeedYee = initialSpeed*cos(abs(kCenter)*dx/2) / ...
        cos(omegaYee*dt/2);
    centerAtMinusHalfStep = xCenter- ...
        propagationSign*groupSpeedYee*dt/2;
    magneticEnvelope = exp(-2*log(2)* ...
        ((xH-centerAtMinusHalfStep)/pulseIntensityFwhm).^2);
    Hhalf0 = propagationSign*sqrt(initialEpsilon/muRelative)* ...
        fieldAmplitude.*magneticEnvelope.* ...
        exp(1i*kCenter*magneticCarrierCoordinate) .* ...
        exp(1i*omegaYee*dt/2);
    carrierOmegaYee(kIndex) = omegaYee;
    carrierGroupSpeedYee(kIndex) = propagationSign*groupSpeedYee;

    if kIndex == representativeKIndex
        recordSpatialIndices = representativeSpatialIndices;
    else
        recordSpatialIndices = probeIndices;
    end

    fdtdCfg = fdtdBase;
    fdtdCfg.E0 = E0;
    fdtdCfg.Hhalf0 = Hhalf0;
    fdtdCfg.recordSpatialIndices = recordSpatialIndices;
    fieldResult = fdtd1d(fdtdCfg);

    if isempty(probeTime)
        probeTime = fieldResult.t;
    elseif ~isequal(probeTime,fieldResult.t)
        error('不同 k 仿真的探针时间网格不一致。');
    end
    [probeIsRecorded,probeColumns] = ismember(probeIndices,recordSpatialIndices);
    if ~all(probeIsRecorded)
        error('内部错误：某些探针没有包含在当前 FDTD 记录中。');
    end
    probeSignals(:,:,kIndex) = fieldResult.E(:,probeColumns);

    if kIndex == representativeKIndex
        representativeField.x = fieldResult.x;
        representativeField.t = fieldResult.t;
        representativeField.E = fieldResult.E;
        representativeField.k = kCenter;
        representativeField.kIndex = kIndex;
        representativeField.courant = fieldResult.courant;
        representativeField.boundaryType = fieldResult.boundaryType;
    end

    if mod(kIndex,progressEvery) == 0 || kIndex == 1 || kIndex == nK
        fprintf('%s：%d/%d，%s=%+.7g\n', ...
            progressLabel,kIndex,nK,progressKLabel,kCenter/progressKScale);
    end
end

result.kScan = kScan;
result.probeSignals = probeSignals;
result.time = probeTime;
result.probeIndices = probeIndices;
result.probePositions = x(probeIndices);
result.representativeField = representativeField;
result.representativeKIndex = representativeKIndex;
result.dx = dx;
result.dt = dt;
result.recordEvery = recordEvery;
result.recordPrecision = recordPrecision;
result.pulseIntensityFwhm = pulseIntensityFwhm;
result.amplitudeGaussianWidth = amplitudeGaussianWidth;
result.initialEpsilon = initialEpsilon;
result.muRelative = muRelative;
result.propagationDirections = propagationDirections;
result.zeroKPropagationDirection = zeroKPropagationDirection;
result.carrierOmegaYee = carrierOmegaYee;
result.carrierGroupSpeedYee = carrierGroupSpeedYee;
result.courant = representativeField.courant;
end

% -------------------------------------------------------------------------
function indices = validate_indices(value,Nx,label)
indices = value(:).';
if ~isnumeric(indices) || isempty(indices) || ~isreal(indices) || ...
        any(~isfinite(indices)) || any(indices ~= round(indices)) || ...
        any(indices < 1) || any(indices > Nx) || any(diff(indices) <= 0)
    error('%s 必须是 [1,Nx] 内严格递增、互不重复的有限整数向量。',label);
end
end

% -------------------------------------------------------------------------
function value = read_positive_scalar(value,label)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value <= 0
    error('%s 必须是有限正实标量。',label);
end
end

% -------------------------------------------------------------------------
function value = read_finite_real_scalar(value,label)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value)
    error('%s 必须是有限实标量。',label);
end
end

% -------------------------------------------------------------------------
function value = read_integer(value,label,minimumValue)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value ~= round(value) || value < minimumValue
    error('%s 必须是不小于 %d 的有限整数。',label,minimumValue);
end
end

% -------------------------------------------------------------------------
function value = read_optional_integer(cfg,fieldName,defaultValue,minimumValue)
if isfield(cfg,fieldName) && ~isempty(cfg.(fieldName))
    value = read_integer(cfg.(fieldName),['cfg.' fieldName],minimumValue);
else
    value = defaultValue;
end
end

% -------------------------------------------------------------------------
function value = read_optional_choice(cfg,fieldName,defaultValue,allowedValues)
if isfield(cfg,fieldName) && ~isempty(cfg.(fieldName))
    value = cfg.(fieldName);
else
    value = defaultValue;
end
if isstring(value) && isscalar(value)
    value = char(value);
elseif ~ischar(value) || size(value,1) ~= 1
    error('cfg.%s 必须是文本标量。',fieldName);
end
matched = strcmpi(value,allowedValues);
if ~any(matched)
    error('cfg.%s 必须是：%s。',fieldName,strjoin(allowedValues,', '));
end
value = allowedValues{find(matched,1)};
end
