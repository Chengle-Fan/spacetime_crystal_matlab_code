function result = fdtd1d(cfg)
%FDTD1D  有限一维介质的 D/B-Yee 时域有限差分求解器。
% 仅处理初值激发与截断空间域；可记录复电场 E 或电位移 D。

if nargin ~= 1 || ~isstruct(cfg) || ~isscalar(cfg)
    error('必须以标量结构体 cfg 调用 fdtd1d。');
end
required = {'x','dt','nSteps','epsFun','E0','Hhalf0'};
% 必需输入分别给出 E 网格、时间步数、材料函数，以及相差半个时间步的
% E(t=0) 和 H(t=-dt/2) 初值。muFun、记录选项、边界与材料下界均为可选。
missing = required(~isfield(cfg,required));
if ~isempty(missing)
    error('cfg 缺少字段：%s。',strjoin(missing,', '));
end

%% ===================== 网格和步数检查 =====================

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
xH = x(1:end-1)+dx/2;

dt = cfg.dt;
nSteps = cfg.nSteps;
if ~isnumeric(dt) || ~isscalar(dt) || ~isreal(dt) || ~isfinite(dt) || dt <= 0
    error('cfg.dt 必须是有限正实标量。');
end
if ~isnumeric(nSteps) || ~isscalar(nSteps) || ~isreal(nSteps) || ...
        ~isfinite(nSteps) || nSteps < 0 || nSteps ~= round(nSteps)
    error('cfg.nSteps 必须是非负整数。');
end
if ~isa(cfg.epsFun,'function_handle')
    error('cfg.epsFun 必须是函数句柄。');
end
if isfield(cfg,'muFun') && ~isempty(cfg.muFun)
    muFun = cfg.muFun;
    if ~isa(muFun,'function_handle')
        error('cfg.muFun 必须是函数句柄。');
    end
else
    muFun = @(position,time) ones(size(position));
end

recordEvery = read_integer_option(cfg,'recordEvery',1,1);
recordField = read_text_option(cfg,'recordField','E',{'E','D'});
recordPrecision = read_text_option(cfg,'recordPrecision','double', ...
    {'double','single'});

% 默认记录整个 E 网格；大规模场图可只给出严格递增的整数下标 ROI，
% 内部推进仍使用全空间双精度场，不因减少输出而改变数值解。
if isfield(cfg,'recordSpatialIndices') && ...
        ~isempty(cfg.recordSpatialIndices)
    recordSpatialIndices = cfg.recordSpatialIndices(:).';
    if ~isnumeric(recordSpatialIndices) || ~isreal(recordSpatialIndices) || ...
            any(~isfinite(recordSpatialIndices)) || ...
            any(recordSpatialIndices ~= round(recordSpatialIndices)) || ...
            any(recordSpatialIndices < 1) || any(recordSpatialIndices > Nx) || ...
            any(diff(recordSpatialIndices) <= 0)
        error(['cfg.recordSpatialIndices 必须是 [1,Nx] 内严格递增、' ...
            '互不重复的有限整数向量。']);
    end
else
    recordSpatialIndices = 1:Nx;
end

% fixed-distant 是固定截断端点，只有在场于仿真结束前到不了边界时才可
% 近似无反射；sponge 是简化阻尼层而非 PML。显式指定 fixed-distant 时
% 默认完全关闭海绵，防止配置名称与实际执行的边界不一致。
boundaryWasSpecified = isfield(cfg,'boundaryType') && ...
    ~isempty(cfg.boundaryType);
if boundaryWasSpecified
    boundaryType = read_text_option(cfg,'boundaryType','fixed-distant', ...
        {'fixed-distant','sponge'});
else
    boundaryType = '';
end
spongeCellsWasSpecified = isfield(cfg,'spongeCells') && ...
    ~isempty(cfg.spongeCells);
spongeStrengthWasSpecified = isfield(cfg,'spongeStrength') && ...
    ~isempty(cfg.spongeStrength);
legacySpongeWasRequested = ~boundaryWasSpecified && ...
    (spongeCellsWasSpecified || spongeStrengthWasSpecified);
if strcmp(boundaryType,'fixed-distant') || ...
        (~boundaryWasSpecified && ~legacySpongeWasRequested)
    defaultSpongeCells = 0;
    defaultSpongeStrength = 0;
else
    defaultSpongeCells = floor(Nx/10);
    defaultSpongeStrength = 8;
end
spongeCells = read_integer_option(cfg,'spongeCells',defaultSpongeCells,0);
if 2*spongeCells >= Nx
    error('两侧 spongeCells 不能重叠。');
end
if isfield(cfg,'spongeStrength') && ~isempty(cfg.spongeStrength)
    spongeStrength = cfg.spongeStrength;
else
    spongeStrength = defaultSpongeStrength;
end
if ~isnumeric(spongeStrength) || ~isscalar(spongeStrength) || ...
        ~isreal(spongeStrength) || ~isfinite(spongeStrength) || spongeStrength < 0
    error('cfg.spongeStrength 必须是非负有限实标量。');
end
if ~boundaryWasSpecified
    % 旧调用常用单独的零值字段关闭海绵；把另一个未显式字段也归零，
    % 同时仍拒绝调用者亲自给出的“零宽度、正强度”等矛盾组合。
    if spongeCellsWasSpecified && spongeCells == 0 && ...
            ~spongeStrengthWasSpecified
        spongeStrength = 0;
    elseif spongeStrengthWasSpecified && spongeStrength == 0 && ...
            ~spongeCellsWasSpecified
        spongeCells = 0;
    end
    if spongeCells > 0 && spongeStrength > 0
        boundaryType = 'sponge';
    else
        boundaryType = 'fixed-distant';
    end
end
if strcmp(boundaryType,'fixed-distant') && ...
        (spongeCells ~= 0 || spongeStrength ~= 0)
    error(['boundaryType=''fixed-distant'' 时 spongeCells 和 ' ...
        'spongeStrength 必须同时为 0。']);
elseif strcmp(boundaryType,'sponge') && ...
        (spongeCells == 0 || spongeStrength == 0)
    error(['boundaryType=''sponge'' 时 spongeCells 和 spongeStrength ' ...
        '必须同时为正。']);
end

% 对解析上已知材料取值范围的算例，调用者可以提供严格正的全局下界，
% 跳过昂贵的全时空预扫描。两个下界必须同时提供；推进期间每次实际材料
% 求值仍会检查是否违反声明，以尽早发现错误的“认证下界”。
hasCertifiedEpsilon = isfield(cfg,'certifiedMinimumEpsilon') && ...
    ~isempty(cfg.certifiedMinimumEpsilon);
hasCertifiedMu = isfield(cfg,'certifiedMinimumMu') && ...
    ~isempty(cfg.certifiedMinimumMu);
if xor(hasCertifiedEpsilon,hasCertifiedMu)
    error(['cfg.certifiedMinimumEpsilon 和 cfg.certifiedMinimumMu ' ...
        '必须同时提供或同时省略。']);
end
if hasCertifiedEpsilon
    epsilonMinimum = read_positive_option(cfg.certifiedMinimumEpsilon, ...
        'cfg.certifiedMinimumEpsilon');
    muMinimum = read_positive_option(cfg.certifiedMinimumMu, ...
        'cfg.certifiedMinimumMu');
    cflAuditSource = 'certified-material-bounds';
else
    epsilonMinimum = inf;
    muMinimum = inf;
    cflAuditSource = 'full-grid-time-audit';
end

E = read_initial_field(cfg.E0,Nx,'cfg.E0');
H = read_initial_field(cfg.Hhalf0,Nx-1,'cfg.Hhalf0');

%% ===================== 时间突变节点检查 =====================

if isfield(cfg,'temporalInterfaces') && ~isempty(cfg.temporalInterfaces)
    temporalInterfaces = cfg.temporalInterfaces(:).';
    if ~isnumeric(temporalInterfaces) || ~isreal(temporalInterfaces) || ...
            any(~isfinite(temporalInterfaces)) || any(diff(temporalInterfaces) <= 0)
        error('cfg.temporalInterfaces 必须是严格递增的有限实向量。');
    end
    interfaceNodes = round(temporalInterfaces/dt);
    alignmentError = abs(temporalInterfaces-interfaceNodes*dt);
    alignmentTolerance = min(dt/4, ...
        max(128*eps(max(1,nSteps*dt)),1e-10*dt));
    if any(alignmentError > alignmentTolerance)
        error('每个时间突变必须严格对齐到整数时间节点 n*dt。');
    end
    if any(diff(interfaceNodes) <= 0)
        error('时间突变舍入到 FDTD 节点后必须仍然严格递增且互不重复。');
    end
    temporalInterfaces = interfaceNodes*dt;
    if any(interfaceNodes < 0) || any(interfaceNodes >= nSteps)
        error('时间突变必须位于 [0,nSteps*dt) 内。');
    end
else
    temporalInterfaces = [];
    interfaceNodes = [];
end
interfaceIdByStep = zeros(1,nSteps,'uint32');
if ~isempty(interfaceNodes)
    interfaceIdByStep(interfaceNodes+1) = uint32(1:numel(interfaceNodes));
end

%% ===================== 初始 D/B 和严格 CFL 审计 =====================

% 理想时间界面上连续的是 D/B，而 E/H 通常跳变，因此内核始终以
% D=epsilon*E、B=mu*H 初始化和推进，仅在需要计算旋度时恢复 E/H。
epsilonNow = evaluate_material(cfg.epsFun,x,0,'epsFun');
muHalf = evaluate_material(muFun,xH,-dt/2,'muFun');
if hasCertifiedEpsilon
    enforce_certified_bound(epsilonNow,epsilonMinimum,'epsFun',0);
    enforce_certified_bound(muHalf,muMinimum,'muFun',-dt/2);
end
D = epsilonNow.*E;
B = muHalf.*H;

% 未提供认证下界时，epsilon 和 mu 分别在实际 E/H Yee 网格与全部推进
% 时刻审计；否则直接采用解析下界。两种路径都以最小 epsilon、mu 组合
% 得到保守波速上界。
if ~hasCertifiedEpsilon
    epsilonAuditTimes = (0:nSteps)*dt;
    if ~isempty(temporalInterfaces)
        epsilonAuditTimes = [epsilonAuditTimes,temporalInterfaces-dt/4, ...
            temporalInterfaces+dt/4];
    end
    for timeValue = epsilonAuditTimes
        epsilonValue = evaluate_material(cfg.epsFun,x,timeValue,'epsFun');
        epsilonMinimum = min(epsilonMinimum,min(epsilonValue));
    end
    muAuditTimes = (-0.5:nSteps-0.5)*dt;
    for timeValue = muAuditTimes
        muValue = evaluate_material(muFun,xH,timeValue,'muFun');
        muMinimum = min(muMinimum,min(muValue));
    end
end
maximumWaveSpeed = 1/sqrt(epsilonMinimum*muMinimum);
courant = maximumWaveSpeed*dt/dx;
if ~isfinite(courant) || courant >= 1
    error('一维 Courant 数 %.6g 必须严格小于 1；请减小 dt 并重新检查网格收敛。',courant);
end

%% ===================== 截断边界与可选海绵层 =====================

% 两种边界最终都不更新 D 的首尾节点。sponge 只是在固定端点之前附加
% exp(-sigma*dt) 阻尼，不是严格开放边界或 PML；在有 k-gap 增益的时间
% 晶体中，空间非均匀阻尼还可能散射并放大额外 k 分量。
dampingE = ones(1,Nx);
dampingH = ones(1,Nx-1);
if spongeCells > 0 && spongeStrength > 0
    distanceE = min(0:Nx-1,Nx-1:-1:0);
    distanceH = min((0:Nx-2)+0.5,(Nx-2:-1:0)+0.5);
    maskE = distanceE < spongeCells;
    maskH = distanceH < spongeCells;
    depthE = (spongeCells-distanceE(maskE))/spongeCells;
    depthH = (spongeCells-distanceH(maskH))/spongeCells;
    dampingE(maskE) = exp(-spongeStrength*dt*depthE.^3);
    dampingH(maskH) = exp(-spongeStrength*dt*depthH.^3);
end

%% ===================== 记录数组 =====================

recordSteps = 0:recordEvery:nSteps;
fieldHistory = complex(zeros(numel(recordSteps), ...
    numel(recordSpatialIndices),recordPrecision));
tHistory = recordSteps(:)*dt;
if strcmp(recordField,'E')
    fieldHistory(1,:) = E(recordSpatialIndices);
else
    fieldHistory(1,:) = D(recordSpatialIndices);
end
nextRecord = 2;

%% ===================== D/B-Yee 时间推进 =====================

% 每一步先用 Faraday 更新半时间格上的 B，再用 Ampere 更新整数时间格上的
% D。截断边界不做首尾环绕；仅 sponge 模式会在固定端点前逐渐阻尼。
for step = 1:nSteps
    % B 从 t=(n-1/2)dt 推到 (n+1/2)dt。若中心 t=n*dt 是材料突变，
    % 用界面两侧 E 的平均值构造中心时间近似；其他时刻直接使用当前 E。
    interfaceId = interfaceIdByStep(step);
    if interfaceId ~= 0
        interfaceTime = temporalInterfaces(double(interfaceId));
        epsilonBefore = evaluate_material(cfg.epsFun,x,interfaceTime-dt/4,'epsFun');
        epsilonAfter = evaluate_material(cfg.epsFun,x,interfaceTime+dt/4,'epsFun');
        if hasCertifiedEpsilon
            enforce_certified_bound(epsilonBefore,epsilonMinimum, ...
                'epsFun',interfaceTime-dt/4);
            enforce_certified_bound(epsilonAfter,epsilonMinimum, ...
                'epsFun',interfaceTime+dt/4);
        end
        electricForB = 0.5*(D./epsilonBefore+D./epsilonAfter);
    else
        electricForB = E;
    end

    B = B-(dt/dx)*diff(electricForB);
    B = B.*dampingH;
    timeHalf = (step-0.5)*dt;
    muHalf = evaluate_material(muFun,xH,timeHalf,'muFun');
    if hasCertifiedMu
        enforce_certified_bound(muHalf,muMinimum,'muFun',timeHalf);
    end
    H = B./muHalf;

    % 截断端点没有环绕差分，只更新内部 D；海绵模式再附加空间阻尼。
    D(2:end-1) = D(2:end-1)-(dt/dx)*(H(2:end)-H(1:end-1));
    D = D.*dampingE;
    timeNow = step*dt;
    epsilonNow = evaluate_material(cfg.epsFun,x,timeNow,'epsFun');
    if hasCertifiedEpsilon
        enforce_certified_bound(epsilonNow,epsilonMinimum,'epsFun',timeNow);
    end
    E = D./epsilonNow;

    if nextRecord <= numel(recordSteps) && step == recordSteps(nextRecord)
        if strcmp(recordField,'E')
            fieldHistory(nextRecord,:) = E(recordSpatialIndices);
        else
            fieldHistory(nextRecord,:) = D(recordSpatialIndices);
        end
        nextRecord = nextRecord+1;
    end
end

% 只保存用户选择的一个观测量。时间节点上的 E 使用材料回调在该节点的
% 右连续取值；理想时间界面上的 D 则由推进变量直接记录并保持连续。
result.x = x(recordSpatialIndices);
result.t = tHistory;
result.(recordField) = fieldHistory;
result.fieldComponent = recordField;
result.recordPrecision = recordPrecision;
result.recordSpatialIndices = recordSpatialIndices;
result.fullDomain = [x(1) x(end)];
result.fullGridPointCount = Nx;
result.dx = dx;
result.dt = dt;
result.courant = courant;
result.maximumWaveSpeed = maximumWaveSpeed;
result.cflAuditSource = cflAuditSource;
result.boundaryType = boundaryType;
result.spongeCells = spongeCells;
result.spongeStrength = spongeStrength;
result.domainTraversalTime = (x(end)-x(1))/maximumWaveSpeed;
end

% -------------------------------------------------------------------------
function field = read_initial_field(value,expectedLength,label)
% 初始场允许为复数，但必须是长度精确匹配且全部有限的向量。
if ~isnumeric(value) || ~isvector(value) || numel(value) ~= expectedLength || ...
        any(~isfinite(value))
    error('%s 必须是含 %d 个有限数值的向量。',label,expectedLength);
end
field = value(:).';
end

% -------------------------------------------------------------------------
function value = evaluate_material(materialFunction,grid,timeValue,label)
% 每次材料求值都执行尺寸、有限性和正实性检查，禁止伪造 NaN/Inf 场。
value = materialFunction(grid,timeValue);
if ~isnumeric(value) || ~isequal(size(value),size(grid)) || ...
        ~isreal(value) || any(~isfinite(value)) || any(value <= 0)
    error('%s 在 t=%.16g 必须返回与网格同尺寸的有限正实数组。',label,timeValue);
end
value = value(:).';
end

% -------------------------------------------------------------------------
function value = read_integer_option(cfg,fieldName,defaultValue,minimumValue)
% 读取可选整数参数，并统一执行有限性和下界检查。
if isfield(cfg,fieldName) && ~isempty(cfg.(fieldName))
    value = cfg.(fieldName);
else
    value = defaultValue;
end
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value < minimumValue || value ~= round(value)
    error('cfg.%s 必须是不小于 %d 的有限整数。',fieldName,minimumValue);
end
end

% -------------------------------------------------------------------------
function value = read_text_option(cfg,fieldName,defaultValue,allowedValues)
% 接受字符行向量或标量 string，并统一成允许列表中的字符文本。
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

% -------------------------------------------------------------------------
function value = read_positive_option(value,label)
% 认证材料下界必须是真正可用于 CFL 上界的严格正有限实标量。
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value <= 0
    error('%s 必须是有限正实标量。',label);
end
end

% -------------------------------------------------------------------------
function enforce_certified_bound(value,certifiedMinimum,label,timeValue)
% 推进时核查实际材料值没有低于调用者声明的全局下界。
tolerance = 128*eps(max(1,certifiedMinimum));
if min(value) < certifiedMinimum-tolerance
    error(['%s 在 t=%.16g 违反认证下界：实际最小值 %.16g < ' ...
        '声明值 %.16g。'],label,timeValue,min(value),certifiedMinimum);
end
end
