%RUN_FDTD_FIELD  计算有限一维光学时间晶体的场分布。
% 第一模块：由复数 E/H 初值推进，只绘制用户选择的 E 或 D 时空幅度。

clear; clc; close all;

%% ===================== 用户参数：单位与材料 =====================

% 'normalized' 表示 c0=1，长度和时间使用彼此一致的归一化单位；
% 'um-fs' 表示长度用 um、时间用 fs，程序会把物理空间换成光行时坐标
% z/c0 后再调用 fdtd1d，避免把 1.4 um 误当成归一化长度。
unitSystem = 'normalized';

% 时间调制可选 'square'（方波）或 'sinusoidal'（正弦）。
modulationType = 'square';
epsHigh = 4;
epsLow = 1;
muRelative = 1;
backgroundEpsilon = 1;
inactiveSampleEpsilon = 1;  % 调制尚未开启或已经关闭时的样品介电常数

temporalPeriod = 1;
dutyCycle = 0.5;            % 方波高介电状态所占周期比例
squarePhaseFraction = 0;    % 开启时刻在方波周期内的位置，范围 [0,1)
sinePhase = 0;              % 正弦调制开启时的初相，单位 rad

%% ===================== 用户参数：仿真时间与有限调制窗口 =====================

stepsPerPeriod = 100;
simulationTime = 12*temporalPeriod;
recordEvery = 2;

% PTC 只在此时间区间内开启。端点必须与 dt 对齐；结束时刻允许等于
% simulationTime。默认在开启前和关闭后各保留两个周期观察传播。
modulationStartTime = 2*temporalPeriod;
modulationEndTime = 10*temporalPeriod;

% 'sample' 只调制下方定义的有限空间样品；'global' 在整个计算域同时
% 调制，适合验证空间均匀的纯时间晶体。
modulationRegion = 'sample';

%% ===================== 用户参数：有限空间区域和边界 =====================

% 纯时间晶体没有空间晶格。cellSize 只是几何长度单位，样品长度等于
% sampleCellCount*cellSize，不能把它当作 Bloch 晶格常数。
cellSize = 1;
sampleCellCount = 40;
pointsPerCell = 20;
backgroundCellCount = 5;

% 'fixed-distant' 是推荐默认值：关闭阻尼，并要求主场在仿真结束前到不了
% 固定截断端点。'sponge' 会在 boundaryLayerCellCount 范围内启用简化
% 阻尼层；它不是 PML，且可能把极小的额外 k 分量散射进增长带隙。
boundaryType = 'fixed-distant';
boundaryLayerCellCount = 3;
spongeStrength = 8;

%% ===================== 用户参数：初值激发与记录 =====================

% 可选 'normalized-k' 或 'wavelength'。前者设置 k*c0/Omega，后者直接
% 设置与 cellSize 相同长度单位的真空中心波长。
excitationInput = 'normalized-k';
kCenterNormalized = 0.70;
centerWavelength = 1.4;

fieldAmplitude = 1;
pulseIntensityFwhm = 1.8;  % 强度 |E|^2 的空间 FWHM，不是含糊的宽度参数

% 只选择一个观测量和一块记录区域，避免无用的 D/E 双历史或全域历史。
fieldComponentToPlot = 'E';   % 'E' 或 'D'
recordRegion = 'full';        % 'full' 或 'sample'
recordPrecision = 'single';   % 只影响历史数组；内部推进始终为 double

%% ===================== 基本参数和单位检查 =====================

unitSystem = validate_choice(unitSystem,{'normalized','um-fs'},'unitSystem');
if strcmp(unitSystem,'normalized')
    vacuumLightSpeed = 1;
    spaceUnitLabel = '归一化长度';
else
    vacuumLightSpeed = 0.299792458;  % um/fs
    spaceUnitLabel = 'z (\mum)';
end

modulationType = validate_choice(modulationType, ...
    {'square','sinusoidal'},'modulationType');
modulationRegion = validate_choice(modulationRegion, ...
    {'sample','global'},'modulationRegion');
boundaryType = validate_choice(boundaryType, ...
    {'fixed-distant','sponge'},'boundaryType');
excitationInput = validate_choice(excitationInput, ...
    {'normalized-k','wavelength'},'excitationInput');
fieldComponentToPlot = validate_choice(fieldComponentToPlot, ...
    {'E','D'},'fieldComponentToPlot');
recordRegion = validate_choice(recordRegion, ...
    {'full','sample'},'recordRegion');
recordPrecision = validate_choice(recordPrecision, ...
    {'double','single'},'recordPrecision');

positiveParameters = {epsHigh,epsLow,muRelative,backgroundEpsilon, ...
    inactiveSampleEpsilon,temporalPeriod,cellSize,pulseIntensityFwhm, ...
    fieldAmplitude};
validPositive = cellfun(@(value) isnumeric(value) && isscalar(value) && ...
    isreal(value) && isfinite(value) && value > 0,positiveParameters);
if ~all(validPositive) || epsHigh < epsLow
    error(['介电常数、mu、周期、几何长度、脉冲 FWHM 和场幅必须为' ...
        '有限正数，并满足 epsHigh >= epsLow。']);
end

integerParameters = {sampleCellCount,pointsPerCell,backgroundCellCount, ...
    boundaryLayerCellCount,stepsPerPeriod,recordEvery};
validIntegers = cellfun(@(value) isnumeric(value) && isscalar(value) && ...
    isreal(value) && isfinite(value) && value == round(value), ...
    integerParameters);
if ~all(validIntegers) || sampleCellCount < 1 || pointsPerCell < 2 || ...
        backgroundCellCount < 1 || boundaryLayerCellCount < 0 || ...
        stepsPerPeriod < 4 || recordEvery < 1
    error('空间/时间离散数、样品数、背景数和边界层数必须是允许范围内的整数。');
end
if ~isnumeric(simulationTime) || ~isscalar(simulationTime) || ...
        ~isreal(simulationTime) || ~isfinite(simulationTime) || ...
        simulationTime <= 0 || ~isnumeric(modulationStartTime) || ...
        ~isscalar(modulationStartTime) || ~isreal(modulationStartTime) || ...
        ~isfinite(modulationStartTime) || ~isnumeric(modulationEndTime) || ...
        ~isscalar(modulationEndTime) || ~isreal(modulationEndTime) || ...
        ~isfinite(modulationEndTime) || modulationStartTime < 0 || ...
        modulationEndTime <= modulationStartTime || ...
        modulationEndTime > simulationTime
    error('调制窗口必须满足 0 <= start < end <= simulationTime。');
end
if ~isnumeric(spongeStrength) || ~isscalar(spongeStrength) || ...
        ~isreal(spongeStrength) || ~isfinite(spongeStrength) || ...
        spongeStrength <= 0
    error('spongeStrength 必须是有限正实数；fixed-distant 模式会自动停用它。');
end

T = temporalPeriod;
Omega = 2*pi/T;
dt = T/stepsPerPeriod;
nStepsFloat = simulationTime/dt;
if abs(nStepsFloat-round(nStepsFloat)) > 1e-10*max(1,nStepsFloat)
    error('simulationTime 必须是 dt 的整数倍。');
end
nSteps = round(nStepsFloat);
if mod(nSteps,recordEvery) ~= 0
    error('总步数必须能被 recordEvery 整除，确保记录包含终止时刻。');
end
if mod(stepsPerPeriod,recordEvery) ~= 0
    error('recordEvery 必须整除 stepsPerPeriod，保证记录与调制周期可公度。');
end

windowNodes = [modulationStartTime modulationEndTime]/dt;
if any(abs(windowNodes-round(windowNodes)) > 1e-10*max(1,nSteps))
    error('modulationStartTime 和 modulationEndTime 必须严格对齐到整数时间节点。');
end

%% ===================== 构造有限时间材料函数 =====================

switch modulationType
    case 'square'
        highStepsFloat = dutyCycle*stepsPerPeriod;
        phaseStepsFloat = squarePhaseFraction*stepsPerPeriod;
        if ~isnumeric(dutyCycle) || ~isscalar(dutyCycle) || ...
                ~isreal(dutyCycle) || ~isfinite(dutyCycle) || ...
                dutyCycle <= 0 || dutyCycle >= 1 || ...
                abs(highStepsFloat-round(highStepsFloat)) > 1e-10
            error('方波 dutyCycle 位于 (0,1)，且 dutyCycle*stepsPerPeriod 必须为整数。');
        end
        if ~isnumeric(squarePhaseFraction) || ...
                ~isscalar(squarePhaseFraction) || ...
                ~isreal(squarePhaseFraction) || ...
                ~isfinite(squarePhaseFraction) || ...
                squarePhaseFraction < 0 || squarePhaseFraction >= 1 || ...
                abs(phaseStepsFloat-round(phaseStepsFloat)) > 1e-10
            error(['squarePhaseFraction 必须位于 [0,1)，并对应整数个' ...
                ' FDTD 时间步。']);
        end
        if round(highStepsFloat) < 2 || ...
                stepsPerPeriod-round(highStepsFloat) < 2
            error('方波高、低介电时间层均至少需要两个 FDTD 步。');
        end

    case 'sinusoidal'
        if ~isnumeric(sinePhase) || ~isscalar(sinePhase) || ...
                ~isreal(sinePhase) || ~isfinite(sinePhase)
            error('sinePhase 必须是有限实标量。');
        end
end

epsilonTime = @(time) finite_temporal_epsilon(time,modulationType, ...
    epsHigh,epsLow,inactiveSampleEpsilon,T,dutyCycle, ...
    squarePhaseFraction,sinePhase,modulationStartTime,modulationEndTime);

% 自动比较每个整数节点左右各 dt/4 的材料值，从而同时捕获调制开启、
% 关闭和方波内部界面。结束时刻等于最终记录节点时不再参加下一次 B 更新。
candidateNodes = 0:nSteps-1;
candidateTimes = candidateNodes*dt;
epsilonBeforeNode = epsilonTime(candidateTimes-dt/4);
epsilonAfterNode = epsilonTime(candidateTimes+dt/4);
interfaceTolerance = 128*eps(max([1 epsHigh epsLow inactiveSampleEpsilon]));
temporalInterfaceNodes = candidateNodes( ...
    abs(epsilonAfterNode-epsilonBeforeNode) > interfaceTolerance);
temporalInterfaces = temporalInterfaceNodes*dt;

%% ===================== 构造空间网格和材料区域 =====================

dz = cellSize/pointsPerCell;
nSample = sampleCellCount*pointsPerCell;
nBackground = backgroundCellCount*pointsPerCell;
nBoundary = boundaryLayerCellCount*pointsPerCell;
Nx = nSample+2*nBackground+2*nBoundary;
z = (0:Nx-1)*dz;

% fdtd1d 的 Maxwell 方程采用 c0=1，因此内部坐标为 x=z/c0；物理绘图仍
% 使用 z。归一化单位下 vacuumLightSpeed=1，两套坐标自动相同。
x = z/vacuumLightSpeed;
dx = x(2)-x(1);
iSampleFirst = nBoundary+nBackground+1;
iSampleLast = iSampleFirst+nSample-1;
sampleEdgesZ = [z(iSampleFirst)-dz/2,z(iSampleLast)+dz/2];
sampleEdgesX = sampleEdgesZ/vacuumLightSpeed;
inSample = @(position) position >= sampleEdgesX(1) & ...
    position < sampleEdgesX(2);

if strcmp(modulationRegion,'sample')
    epsilonFunction = @(position,time) backgroundEpsilon+ ...
        (epsilonTime(time)-backgroundEpsilon).*double(inSample(position));
else
    epsilonFunction = @(position,time) epsilonTime(time).*ones(size(position));
end
muFunction = @(position,time) muRelative*ones(size(position));

%% ===================== 把用户激发参数转换为一致的 k 和波长 =====================

switch excitationInput
    case 'normalized-k'
        if ~isnumeric(kCenterNormalized) || ~isscalar(kCenterNormalized) || ...
                ~isreal(kCenterNormalized) || ~isfinite(kCenterNormalized)
            error('kCenterNormalized 必须是有限实标量。');
        end
        kPhysical = kCenterNormalized*Omega/vacuumLightSpeed;
        if kPhysical == 0
            centerWavelength = inf;
        else
            centerWavelength = 2*pi/abs(kPhysical);
        end

    case 'wavelength'
        if ~isnumeric(centerWavelength) || ~isscalar(centerWavelength) || ...
                ~isreal(centerWavelength) || ~isfinite(centerWavelength) || ...
                centerWavelength <= 0
            error('centerWavelength 必须是有限正实标量。');
        end
        kPhysical = 2*pi/centerWavelength;
        kCenterNormalized = vacuumLightSpeed*kPhysical/Omega;
end
kSolver = vacuumLightSpeed*kPhysical;
if abs(kSolver) >= pi/dx
    error('中心波数超过 Yee 空间 Nyquist 上限；请减小 cellSize/pointsPerCell。');
end

%% ===================== 构造具有明确强度 FWHM 的单向波包 =====================

zCenter = mean(sampleEdgesZ);
xCenter = zCenter/vacuumLightSpeed;
amplitudeGaussianWidth = pulseIntensityFwhm/sqrt(2*log(2));
if abs(kPhysical)+6/amplitudeGaussianWidth >= pi/dz
    error('高斯波包主谱接近空间 Nyquist 边界；请细化空间网格或增大 FWHM。');
end

% 有限空间样品模式要求初值在两个材料端面前充分衰减，否则“局部均匀
% 阻抗初值”会跨越空间介质界面而不再是单向波包。
if strcmp(modulationRegion,'sample')
    sampleLength = sampleCellCount*cellSize;
    edgeRelativeAmplitude = exp(-(sampleLength/(2*amplitudeGaussianWidth))^2);
    if edgeRelativeAmplitude > 1e-4
        error('初始波包在样品端面仍不可忽略；请增大样品或减小脉冲 FWHM。');
    end
end

electricEnvelope = exp(-2*log(2)*((z-zCenter)/pulseIntensityFwhm).^2);
E0 = fieldAmplitude*electricEnvelope.*exp(1i*kPhysical*(z-zCenter));

xH = x(1:end-1)+dx/2;
zH = vacuumLightSpeed*xH;
epsInitialH = epsilonTime(-dt/2);
initialSpeedSolver = 1/sqrt(epsInitialH*muRelative);
yeeArgument = initialSpeedSolver*dt/dx*sin(abs(kSolver)*dx/2);
if abs(yeeArgument) >= 1
    error('初始载波与网格不满足 Yee 离散色散条件；请减小 dt 或 dx。');
end
omegaYee = (2/dt)*asin(yeeArgument);
propagationSign = sign(kSolver);

% H 位于 t=-dt/2；除载波半时间相位外，高斯中心也按初始群速回退半步。
initialSpeedPhysical = vacuumLightSpeed/sqrt(epsInitialH*muRelative);
centerAtMinusHalfStep = zCenter- ...
    propagationSign*initialSpeedPhysical*dt/2;
magneticEnvelope = exp(-2*log(2)* ...
    ((zH-centerAtMinusHalfStep)/pulseIntensityFwhm).^2);
Hhalf0 = propagationSign*sqrt(epsInitialH/muRelative)*fieldAmplitude .* ...
    magneticEnvelope.*exp(1i*kSolver*(xH-xCenter)) .* ...
    exp(1i*omegaYee*dt/2);

%% ===================== 边界安全检查与记录 ROI =====================

minimumEpsilon = min([epsHigh epsLow backgroundEpsilon inactiveSampleEpsilon]);
maximumSpeedPhysical = vacuumLightSpeed/sqrt(minimumEpsilon*muRelative);
if strcmp(boundaryType,'fixed-distant')
    distanceToBoundary = min(zCenter-z(1),z(end)-zCenter);
    initialTailMargin = 6*amplitudeGaussianWidth;
    if maximumSpeedPhysical*simulationTime+initialTailMargin >= ...
            distanceToBoundary
        error(['fixed-distant 边界距离不足：主场或可分辨高斯尾场可能在' ...
            '仿真结束前到达端点。请增加背景/边界层或缩短仿真时间。']);
    end
    spongeCellsForSolver = 0;
    spongeStrengthForSolver = 0;
else
    if nBoundary < 1
        error('sponge 边界要求 boundaryLayerCellCount 至少为 1。');
    end
    warning(['当前 sponge 不是 PML；在 k-gap 增益问题中必须改变边界' ...
        '距离、宽度和强度检查谱泄漏收敛。']);
    spongeCellsForSolver = nBoundary;
    spongeStrengthForSolver = spongeStrength;
end

if strcmp(recordRegion,'sample')
    recordSpatialIndices = iSampleFirst:iSampleLast;
else
    recordSpatialIndices = 1:Nx;
end

%% ===================== D/B-Yee 推进 =====================

fdtdCfg = struct();
fdtdCfg.x = x;
fdtdCfg.dt = dt;
fdtdCfg.nSteps = nSteps;
fdtdCfg.epsFun = epsilonFunction;
fdtdCfg.muFun = muFunction;
fdtdCfg.E0 = E0;
fdtdCfg.Hhalf0 = Hhalf0;
fdtdCfg.recordEvery = recordEvery;
fdtdCfg.recordField = fieldComponentToPlot;
fdtdCfg.recordPrecision = recordPrecision;
fdtdCfg.recordSpatialIndices = recordSpatialIndices;
fdtdCfg.boundaryType = boundaryType;
fdtdCfg.spongeCells = spongeCellsForSolver;
fdtdCfg.spongeStrength = spongeStrengthForSolver;
fdtdCfg.temporalInterfaces = temporalInterfaces;

% 当前材料的全局下界可由用户参数解析证明。内核仍会在每次实际求值时
% 核查材料没有低于声明值，从而避免为简单方波重复扫描全部时空网格。
fdtdCfg.certifiedMinimumEpsilon = minimumEpsilon;
fdtdCfg.certifiedMinimumMu = muRelative;
fieldResult = fdtd1d(fdtdCfg);

%% ===================== 只绘制所选场分布 =====================

recordedZ = vacuumLightSpeed*fieldResult.x;
fieldHistory = fieldResult.(fieldComponentToPlot);
figure('Color','w','Position',[100 80 980 620]);
imagesc(recordedZ,fieldResult.t/T,abs(fieldHistory));
axis xy;
colormap(parula(256));
fieldColorbar = colorbar;
fieldColorbar.Label.String = sprintf('|%s(z,t)|',fieldComponentToPlot);
hold on;
if strcmp(modulationRegion,'sample') && ...
        sampleEdgesZ(1) >= recordedZ(1) && sampleEdgesZ(2) <= recordedZ(end)
    xline(sampleEdgesZ(1),'w--','LineWidth',1.2);
    xline(sampleEdgesZ(2),'w--','LineWidth',1.2);
end
xlabel(spaceUnitLabel);
ylabel('归一化时间 t/T');
title(sprintf(['有限时间晶体场分布：k c_0/\\Omega=%.3f，' ...
    '\\lambda_0=%.4g'],kCenterNormalized,centerWavelength));
box on;

fprintf(['FDTD 场分布：记录 %s，边界=%s，Courant=%.6g，' ...
    '调制窗口=[%.6g,%.6g]。\n'],fieldComponentToPlot, ...
    fieldResult.boundaryType,fieldResult.courant, ...
    modulationStartTime,modulationEndTime);

% -------------------------------------------------------------------------
function value = validate_choice(value,allowedValues,label)
% 把字符向量或标量 string 统一为允许列表中的标准字符文本。
if isstring(value) && isscalar(value)
    value = char(value);
elseif ~ischar(value) || size(value,1) ~= 1
    error('%s 必须是文本标量。',label);
end
matched = strcmpi(value,allowedValues);
if ~any(matched)
    error('%s 必须是：%s。',label,strjoin(allowedValues,', '));
end
value = allowedValues{find(matched,1)};
end

% -------------------------------------------------------------------------
function epsilon = finite_temporal_epsilon(time,modulationType, ...
        epsHigh,epsLow,inactiveEpsilon,T,dutyCycle,squarePhaseFraction, ...
        sinePhase,startTime,endTime)
% 在有限开启窗口外返回静态介电常数；窗口内保留输入时间数组的尺寸。
epsilon = inactiveEpsilon*ones(size(time));
active = time >= startTime & time < endTime;
if ~any(active,'all')
    return;
end

switch modulationType
    case 'square'
        phase = mod((time(active)-startTime)/T+squarePhaseFraction,1);
        epsilon(active) = epsLow+(epsHigh-epsLow).*double(phase < dutyCycle);
    case 'sinusoidal'
        epsMean = (epsHigh+epsLow)/2;
        epsAmplitude = (epsHigh-epsLow)/2;
        phase = 2*pi*(time(active)-startTime)/T+sinePhase;
        epsilon(active) = epsMean+epsAmplitude*cos(phase);
end
end
