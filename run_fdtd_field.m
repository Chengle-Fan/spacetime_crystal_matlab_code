%RUN_FDTD_FIELD  计算有限光学时间晶体的场分布。
% 第一模块：由复数 E/H 初值计算有限样品，只绘制 |E(x,t)|。

clear; clc; close all;

%% ===================== 用户参数：样品材料 =====================

% 时间调制可选 'square'（方波）或 'sinusoidal'（正弦）。
modulationType = 'square';

% 样品内部介电常数范围和磁导率；所有值必须是有限正实数。
epsHigh = 4;
epsLow = 1;
muRelative = 1;
backgroundEpsilon = 1;

% 时间调制参数。
temporalPeriod = 1;
dutyCycle = 0.5;       % 方波高介电状态的时间占比
sinePhase = 0;         % 正弦调制初相，单位为弧度

%% ===================== 用户参数：有限样品几何 =====================

% 纯时间晶体在样品内部沿 x 均匀，没有真实的空间晶格。这里的 cellSize
% 只是控制样品长度的几何单位，sampleCellCount*cellSize 才是样品长度；
% 后续 FFT 不会把 cellSize 当作晶格常数，也不会据此折叠空间波数。
cellSize = 1;              % 几何长度单元，不是空间晶格常数
sampleCellCount = 40;      % 有限样品长度 = sampleCellCount*cellSize
pointsPerCell = 20;        % 每个几何长度单元内的 Yee 网格数
backgroundCellCount = 5;   % 样品每侧的静态背景长度
spongeCellCount = 3;       % 每侧海绵吸收层长度；海绵不是 PML
spongeStrength = 0.08;

%% ===================== 用户参数：初值激发 =====================

% 'bulk' 选择体带中的中心波数；'k-gap' 选择动量带隙中的中心波数。
% 材料参数改变后带隙位置也会改变，下面的波数必须相应调整。程序会用
% 方波 TMM 或正弦 PWE 做只读验证，但不会在 FDTD 图中叠加理论结果。
excitationRegion = 'k-gap';
bulkKNormalized = 0.30;     % k/Omega，默认材料下位于体带
kGapKNormalized = 0.70;     % k/Omega，默认材料下位于第一 k-gap

fieldAmplitude = 1;
packetWidthCells = 1.5;     % Gaussian 包络宽度；越窄则激发的 k 范围越宽

%% ===================== 用户参数：时间离散与记录 =====================

stepsPerPeriod = 100;       % dt = temporalPeriod/stepsPerPeriod
simulationTime = 12*temporalPeriod;
recordEvery = 2;            % 每隔多少个 FDTD 时间步记录一次复电场

% 第一模块保存完整复数 E(x,t)，第二模块读取它做 FFT。这里只保存 MAT
% 数组，不保存场分布图片；对 PNG 或 abs(E) 做 FFT 会丢失相位信息。
% MAT 文件已由 .gitignore 排除，不会进入版本库。
dataFile = 'fdtd_field_data.mat';

%% ===================== 参数检查与材料定义 =====================

materialParameters = [epsHigh epsLow muRelative backgroundEpsilon temporalPeriod];
if numel(materialParameters) ~= 5 || ~isreal(materialParameters) || ...
        any(~isfinite(materialParameters)) || epsLow <= 0 || ...
        epsHigh < epsLow || muRelative <= 0 || backgroundEpsilon <= 0 || ...
        temporalPeriod <= 0
    error('必须满足 epsHigh >= epsLow > 0，且 mu、背景 epsilon 和 T 均为正。');
end
integerParameters = [sampleCellCount pointsPerCell backgroundCellCount ...
    spongeCellCount stepsPerPeriod recordEvery];
if any(~isfinite(integerParameters)) || any(integerParameters ~= round(integerParameters)) || ...
        sampleCellCount < 1 || pointsPerCell < 2 || backgroundCellCount < 1 || ...
        spongeCellCount < 1 || stepsPerPeriod < 4 || recordEvery < 1
    error('原胞数、网格数、背景/海绵长度和时间步参数必须是允许范围内的整数。');
end
if ~isnumeric(cellSize) || ~isscalar(cellSize) || ~isreal(cellSize) || ...
        ~isfinite(cellSize) || cellSize <= 0 || ~isnumeric(simulationTime) || ...
        ~isscalar(simulationTime) || ~isreal(simulationTime) || ...
        ~isfinite(simulationTime) || simulationTime <= 0 || ...
        ~isnumeric(packetWidthCells) || ~isscalar(packetWidthCells) || ...
        ~isreal(packetWidthCells) || ~isfinite(packetWidthCells) || ...
        packetWidthCells <= 0 || ~isnumeric(fieldAmplitude) || ...
        ~isscalar(fieldAmplitude) || ~isreal(fieldAmplitude) || ...
        ~isfinite(fieldAmplitude) || fieldAmplitude <= 0 || ...
        ~isnumeric(spongeStrength) || ~isscalar(spongeStrength) || ...
        ~isreal(spongeStrength) || ~isfinite(spongeStrength) || spongeStrength < 0
    error('cellSize、simulationTime、packetWidthCells 和 fieldAmplitude 必须为正。');
end
if mod(stepsPerPeriod,recordEvery) ~= 0
    error('recordEvery 必须整除 stepsPerPeriod，保证记录采样与调制周期可公度。');
end

T = temporalPeriod;
Omega = 2*pi/T;
dt = T/stepsPerPeriod;
nStepsFloat = simulationTime/dt;
if abs(nStepsFloat-round(nStepsFloat)) > 1e-10*max(1,nStepsFloat)
    error('simulationTime 必须是 FDTD 时间步 dt 的整数倍。');
end
nSteps = round(nStepsFloat);
tEnd = nSteps*dt;

switch lower(modulationType)
    case 'square'
        highStepsFloat = dutyCycle*stepsPerPeriod;
        if dutyCycle <= 0 || dutyCycle >= 1 || ...
                abs(highStepsFloat-round(highStepsFloat)) > 1e-10
            error('方波 dutyCycle*stepsPerPeriod 必须为整数，且 dutyCycle 位于 (0,1)。');
        end
        highSteps = round(highStepsFloat);
        lowSteps = stepsPerPeriod-highSteps;
        if highSteps < 2 || lowSteps < 2
            error('方波高、低介电状态均至少需要两个 FDTD 时间步。');
        end

        % 将 t=0 放在高介电时间层内部，而不是放在突变界面上。这样 E(t=0)
        % 和 H(t=-dt/2) 属于同一材料状态，初始阻抗关系不会跨越时间突变。
        phaseSteps = floor(highSteps/2);
        phaseShift = phaseSteps*dt;
        highDuration = highSteps*dt;
        epsilonTime = @(time) epsLow + (epsHigh-epsLow).* ...
            double(mod(time+phaseShift,T) < highDuration);

        % 两组突变时刻分别对应高->低和低->高，且全部严格落在整数时间节点。
        firstHighToLow = highDuration-phaseShift;
        firstLowToHigh = T-phaseShift;
        temporalInterfaces = sort([firstHighToLow:T:tEnd-dt/2, ...
            firstLowToHigh:T:tEnd-dt/2]);

    case 'sinusoidal'
        if ~isnumeric(sinePhase) || ~isscalar(sinePhase) || ...
                ~isreal(sinePhase) || ~isfinite(sinePhase)
            error('sinePhase 必须是有限实标量。');
        end
        epsMean = (epsHigh+epsLow)/2;
        epsAmplitude = (epsHigh-epsLow)/2;
        epsilonTime = @(time) epsMean + ...
            epsAmplitude*cos(Omega*time+sinePhase);
        temporalInterfaces = [];

    otherwise
        error('modulationType 必须是 ''square'' 或 ''sinusoidal''。');
end

switch lower(excitationRegion)
    case 'bulk'
        kCenterNormalized = bulkKNormalized;
    case 'k-gap'
        kCenterNormalized = kGapKNormalized;
    otherwise
        error('excitationRegion 必须是 ''bulk'' 或 ''k-gap''。');
end
if ~isscalar(kCenterNormalized) || ~isreal(kCenterNormalized) || ...
        ~isfinite(kCenterNormalized) || kCenterNormalized == 0
    error('所选中心波数 k/Omega 必须是非零有限实标量。');
end
kCenter = kCenterNormalized*Omega;

% 只读检查中心波数是否真的位于用户声明的区域。方波使用精确两层 TMM；
% 正弦调制使用收敛阶数适中的 PWE。该检查不参与 FFT，也不生成理论曲线。
growthNormalized = validate_excitation_region(modulationType, ...
    excitationRegion,kCenter,epsilonTime,epsHigh,epsLow,muRelative, ...
    dutyCycle,T,Omega);

%% ===================== 构造有限空间区域 =====================

dx = cellSize/pointsPerCell;
nSample = sampleCellCount*pointsPerCell;
nBackground = backgroundCellCount*pointsPerCell;
nSponge = spongeCellCount*pointsPerCell;
Nx = nSample + 2*nBackground + 2*nSponge;
x = (0:Nx-1)*dx;

% 样品界面位于相邻 E 节点之间，使样品恰好包含 nSample 个 E 网格点，
% 几何长度严格等于 sampleCellCount*cellSize。
iSampleFirst = nSponge+nBackground+1;
iSampleLast = iSampleFirst+nSample-1;
sampleEdges = [x(iSampleFirst)-dx/2, x(iSampleLast)+dx/2];
inSample = @(position) position >= sampleEdges(1) & position < sampleEdges(2);
epsilonFunction = @(position,time) backgroundEpsilon + ...
    (epsilonTime(time)-backgroundEpsilon).*double(inSample(position));
muFunction = @(position,time) muRelative*ones(size(position));

%% ===================== 构造满足 Yee 相位关系的初值 =====================

% Gaussian 波包位于样品中心。它不是严格 Floquet 本征态，会投影到两条
% Floquet 分支；在 k-gap 内长时间后增长分支通常占主导。E 和 H 同时设置，
% 避免 H=0 初值无意中产生等强的双向驻波。
xCenter = mean(sampleEdges);
packetWidth = packetWidthCells*cellSize;
E0 = fieldAmplitude*exp(-((x-xCenter)/packetWidth).^2) .* ...
    exp(1i*kCenter*(x-xCenter));

xH = x(1:end-1)+dx/2;
epsInitial = epsilonTime(0);
epsInitialH = epsilonTime(-dt/2);
initialSpeed = 1/sqrt(epsInitial*muRelative);
yeeArgument = initialSpeed*dt/dx*sin(abs(kCenter)*dx/2);
if abs(yeeArgument) > 1+1e-12
    error('初始波数超出当前 Yee 网格的可传播范围；请减小 dt 或 dx。');
end
omegaYee = (2/dt)*asin(max(-1,min(1,yeeArgument)));
propagationSign = sign(kCenter);
Hhalf0 = propagationSign*sqrt(epsInitialH/muRelative)*fieldAmplitude .* ...
    exp(-((xH-xCenter)/packetWidth).^2) .* ...
    exp(1i*kCenter*(xH-xCenter)) .* exp(1i*omegaYee*dt/2);

%% ===================== FDTD 计算 =====================

fdtdCfg = struct();
fdtdCfg.x = x;
fdtdCfg.dt = dt;
fdtdCfg.nSteps = nSteps;
fdtdCfg.epsFun = epsilonFunction;
fdtdCfg.muFun = muFunction;
fdtdCfg.E0 = E0;
fdtdCfg.Hhalf0 = Hhalf0;
fdtdCfg.recordEvery = recordEvery;
fdtdCfg.spongeCells = nSponge;
fdtdCfg.spongeStrength = spongeStrength;
fdtdCfg.temporalInterfaces = temporalInterfaces;
fieldResult = fdtd1d(fdtdCfg);

%% ===================== 保存第二模块所需的复场数据 =====================

% 不保存匿名材料函数，只保存数值数组和解释 FFT 所必需的标量元数据。
fdtdData = struct();
fdtdData.x = fieldResult.x;
fdtdData.t = fieldResult.t;
fdtdData.E = fieldResult.E;
fdtdData.temporalPeriod = T;
fdtdData.Omega = Omega;
fdtdData.sampleEdges = sampleEdges;
fdtdData.cellSize = cellSize;
fdtdData.sampleCellCount = sampleCellCount;
fdtdData.targetRegion = lower(excitationRegion);
fdtdData.kCenterNormalized = kCenterNormalized;
fdtdData.modulationType = lower(modulationType);
fdtdData.material = struct('epsHigh',epsHigh,'epsLow',epsLow, ...
    'muRelative',muRelative,'backgroundEpsilon',backgroundEpsilon, ...
    'dutyCycle',dutyCycle,'sinePhase',sinePhase);
fdtdData.dx = fieldResult.dx;
fdtdData.dt = fieldResult.dt;
fdtdData.courant = fieldResult.courant;
save(dataFile,'fdtdData','-v7.3');

%% ===================== 只绘制场分布 =====================

figure('Color','w','Position',[100 80 980 620]);
imagesc(fieldResult.x,fieldResult.t/T,abs(fieldResult.E));
axis xy;
colormap(parula(256));
fieldColorbar = colorbar;
fieldColorbar.Label.String = '|E(x,t)|';
hold on;
xline(sampleEdges(1),'w--','LineWidth',1.2);
xline(sampleEdges(2),'w--','LineWidth',1.2);
xlabel('空间坐标 x');
ylabel('归一化时间 t/T');
title(sprintf('有限时间晶体场分布：%s，k_0/Ω=%.3f', ...
    excitationRegion,kCenterNormalized));
box on;

% -------------------------------------------------------------------------
function growthNormalized = validate_excitation_region(modulationType, ...
        excitationRegion,kCenter,epsilonTime,epsHigh,epsLow,muRelative, ...
        dutyCycle,T,Omega)
% 用独立理论求解器验证用户指定的 k，而不把理论曲线混入 FDTD/FFT 图。
if strcmpi(modulationType,'square')
    regionCheck = tmm_bands(kCenter,[epsHigh epsLow],muRelative, ...
        T*[dutyCycle 1-dutyCycle]);
else
    % pwe_bands 固定采用 mu_r=1；当实际 mu 为时间不变常数时，频率问题可
    % 等价写成 muRelative*epsilon(t) 与单位磁导率的形式。
    effectiveEpsilonTime = @(time) muRelative*epsilonTime(time);
    checkCfg = struct('T',T,'Mtime',15,'Nt',1024);
    checkFourier = pwe_fourier(effectiveEpsilonTime,checkCfg);
    regionCheck = pwe_bands(checkFourier,kCenter);
end
growthNormalized = max(abs(imag(regionCheck.omega)),[],'all')/Omega;
gapTolerance = 1e-6;
if strcmpi(excitationRegion,'bulk') && growthNormalized > gapTolerance
    error('当前 bulkKNormalized 实际位于 k-gap；请根据材料参数重新设置。');
elseif strcmpi(excitationRegion,'k-gap') && growthNormalized <= gapTolerance
    error('当前 kGapKNormalized 不在可分辨的 k-gap 内；请根据材料参数重新设置。');
end
end
