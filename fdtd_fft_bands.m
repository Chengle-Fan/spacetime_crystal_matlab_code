function result = fdtd_fft_bands(kScan,epsilonTime,muRelative,cfg)
%FDTD_FFT_BANDS  逐 k 推进有限个时间原胞，并由等相位采样重建准频率。
% 每个 k 使用两组独立 D/B 初态；输出只保留第一时间 Floquet 区绘图数据。

if nargin ~= 4 || ~isstruct(cfg) || ~isscalar(cfg)
    error('必须以 kScan、epsilonTime、muRelative 和标量结构体 cfg 调用。');
end
required = {'temporalPeriod','stepsPerPeriod','temporalCellCount','spatialStep'};
missing = required(~isfield(cfg,required));
if ~isempty(missing)
    error('cfg 缺少字段：%s。',strjoin(missing,', '));
end

%% ===================== 波数、材料和离散参数检查 =====================

if ~isnumeric(kScan) || ~isvector(kScan) || isempty(kScan) || ...
        ~isreal(kScan) || any(~isfinite(kScan))
    error('kScan 必须是非空有限实向量。');
end
kScan = kScan(:).';

if ~isa(epsilonTime,'function_handle')
    error('epsilonTime 必须是介电常数随时间变化的函数句柄。');
end
if ~isnumeric(muRelative) || ~isscalar(muRelative) || ~isreal(muRelative) || ...
        ~isfinite(muRelative) || muRelative <= 0
    error('muRelative 必须是有限正实标量。');
end

T = read_positive_scalar(cfg.temporalPeriod,'cfg.temporalPeriod');
dx = read_positive_scalar(cfg.spatialStep,'cfg.spatialStep');
stepsPerPeriod = read_integer(cfg.stepsPerPeriod,'cfg.stepsPerPeriod',4);
temporalCellCount = read_integer(cfg.temporalCellCount, ...
    'cfg.temporalCellCount',4);

if isfield(cfg,'zeroPaddingFactor') && ~isempty(cfg.zeroPaddingFactor)
    zeroPaddingFactor = read_integer(cfg.zeroPaddingFactor, ...
        'cfg.zeroPaddingFactor',1);
else
    zeroPaddingFactor = 1;
end
if isfield(cfg,'dynamicRangeDb') && ~isempty(cfg.dynamicRangeDb)
    dynamicRangeDb = read_positive_scalar(cfg.dynamicRangeDb, ...
        'cfg.dynamicRangeDb');
else
    dynamicRangeDb = 60;
end
if isfield(cfg,'returnTemporalSamples') && ...
        ~isempty(cfg.returnTemporalSamples)
    returnTemporalSamples = cfg.returnTemporalSamples;
    if ~(islogical(returnTemporalSamples) && isscalar(returnTemporalSamples))
        error('cfg.returnTemporalSamples 必须是逻辑标量。');
    end
else
    returnTemporalSamples = false;
end

% 时间突变节点使用一个周期内的零基整数编号：0 表示 t=0，最大允许值为
% stepsPerPeriod-1。平滑调制不需要提供该字段。
if isfield(cfg,'temporalInterfaceSteps') && ...
        ~isempty(cfg.temporalInterfaceSteps)
    interfaceSteps = cfg.temporalInterfaceSteps(:).';
    if ~isnumeric(interfaceSteps) || ~isreal(interfaceSteps) || ...
            any(~isfinite(interfaceSteps)) || ...
            any(interfaceSteps ~= round(interfaceSteps)) || ...
            any(interfaceSteps < 0) || any(interfaceSteps >= stepsPerPeriod) || ...
            any(diff(interfaceSteps) <= 0)
        error(['cfg.temporalInterfaceSteps 必须是 [0,stepsPerPeriod) 内' ...
            '严格递增的整数向量。']);
    end
else
    interfaceSteps = [];
end

dt = T/stepsPerPeriod;
Omega = 2*pi/T;

% 纯时间晶体没有空间晶格，但 Yee 网格仍有空间 Nyquist 上限。等于边界的
% +pi/dx 与 -pi/dx 无法区分，因此这里采用严格不等号。
if any(abs(kScan) >= pi/dx)
    error('所有扫描波数都必须严格满足 |k| < pi/dx。');
end

%% ===================== 一个时间原胞的本构系数 =====================

% 节点 0 和节点 stepsPerPeriod 是相邻原胞的同一调制相位。显式检查其
% 周期一致性，防止非周期材料被误当成 Floquet 系统。
nodeTimes = (0:stepsPerPeriod)*dt;
epsilonNodes = evaluate_epsilon(epsilonTime,nodeTimes);
periodTolerance = 1e-10*max(1,max(abs(epsilonNodes([1 end]))));
if abs(epsilonNodes(end)-epsilonNodes(1)) > periodTolerance
    error('epsilonTime(0) 与 epsilonTime(T) 不一致，材料必须以 T 为周期。');
end

% 通常的 B 更新使用 E^n=D^n/epsilon(n*dt)。若 n*dt 恰是方波时间
% 界面，积分区间在突变前后各占半步，应改用 (E^-+E^+)/2。
inverseEpsilonForB = 1./epsilonNodes(1:end-1);
epsilonMinimum = min(epsilonNodes);
for interfaceStep = interfaceSteps
    interfaceTime = interfaceStep*dt;
    epsilonBefore = evaluate_epsilon(epsilonTime,interfaceTime-dt/4);
    epsilonAfter = evaluate_epsilon(epsilonTime,interfaceTime+dt/4);
    inverseEpsilonForB(interfaceStep+1) = ...
        0.5*(1/epsilonBefore+1/epsilonAfter);
    epsilonMinimum = min([epsilonMinimum epsilonBefore epsilonAfter]);
end

% 这个条件是全 Yee 网格的保守一维 CFL 条件，而不仅是所选 k 点的局部
% 振子条件。因此改变扫描范围不会掩盖空间离散本身的不稳定。
courant = dt/(dx*sqrt(epsilonMinimum*muRelative));
if ~isfinite(courant) || courant >= 1
    error('一维 Courant 数 %.6g 必须严格小于 1；请减小 dt 并重新检查网格收敛。', ...
        courant);
end

%% ===================== 两组独立初态的逐 k D/B-Yee 推进 =====================

% 对 exp(i*k*x) 的 Yee 空间差分并不是连续 k，而是下式的有符号离散
% 波数。保留它可以让本重建与有限空间 FDTD 使用完全相同的数值色散。
kYee = (2/dx)*sin(kScan*dx/2);
nK = numel(kScan);

% 每个 k 的 Maxwell 状态只有 D 和 B 两个自由度。下面采用两组等离散
% 能量的正交基；它们张成完整状态空间，通常能显著降低单一初态造成的
% 漏支风险，也不会任意偏重电初态或磁初态。特殊 E-dark 点仍应通过
% 改变原胞内采样相位检查，不能把任意谱权重当作本征态强度。
D = complex(zeros(2,nK));
B = complex(zeros(2,nK));
D(1,:) = sqrt(epsilonNodes(1));
B(2,:) = sqrt(muRelative);

% 只在 pT（p=0,1,...）记录 E，即有限个时间原胞都取完全相同的调制
% 相位。这对应于有限空间晶体中“逐原胞在同一相对位置采样”的做法。
electricSamples = complex(zeros(temporalCellCount,2,nK));
electricSamples(1,:,:) = reshape(D/epsilonNodes(1),[1 2 nK]);

for cellIndex = 2:temporalCellCount
    for localStep = 1:stepsPerPeriod
        % Faraday：B^(n+1/2)=B^(n-1/2)-i*dt*kYee*E^n。
        electricForB = D*inverseEpsilonForB(localStep);
        B = B-1i*dt*(electricForB.*kYee);

        % Ampere：D^(n+1)=D^n-i*dt*kYee*H^(n+1/2)。mu 在当前
        % 模型中为常数；D/B 在理想时间界面保持连续。
        H = B/muRelative;
        D = D-1i*dt*(H.*kYee);
    end

    electricSamples(cellIndex,:,:) = ...
        reshape(D/epsilonNodes(1),[1 2 nK]);

    % k-gap 中指数增长是 Floquet 物理而非 CFL 发散，但过长记录仍可能
    % 超出浮点范围。此处不逐时归一化，以免人为篡改实际 FFT 线形。
    if any(~isfinite(D(:))) || any(~isfinite(B(:)))
        error(['场在时间原胞推进中超出浮点范围；请减少 temporalCellCount，' ...
            '并检查该 k 区域是否存在强 Floquet 增长。']);
    end
end

%% ===================== 沿时间原胞序列进行 FFT =====================

% 两组完备初态的瞬时功率和不依赖任意相对相位。对记录后半段的对数
% 振幅作最小二乘直线拟合，可得到占主导 Floquet 乘子的每周期增长指数。
% 通带内存在有界拍频，因此其微小非零拟合值只能当数值诊断；k-gap 中
% 经过足够原胞后，该值才近似 max Im(omega*T)，不能恢复衰减支。
sampleIndex = (0:temporalCellCount-1).';
sampleAmplitude = sqrt(reshape(sum(abs(electricSamples).^2,2), ...
    temporalCellCount,nK));
fitRows = floor(temporalCellCount/2)+1:temporalCellCount;
fitCoordinate = sampleIndex(fitRows);
fitCoordinate = fitCoordinate-mean(fitCoordinate);
fitDenominator = sum(fitCoordinate.^2);
logSampleAmplitude = log(max(sampleAmplitude(fitRows,:),realmin));
dominantGrowthPerPeriod = ...
    (fitCoordinate.'*logSampleAmplitude)/fitDenominator;

% 对有限 pT 序列加周期型 Hann 窗以减小截断泄漏。归一化准频率的原生
% bin 间隔为 1/temporalCellCount，Hann 主瓣会进一步展宽；零填充只增加
% 绘图采样点，不能把有限时间记录变成无限时间本征值求解。
window = 0.5-0.5*cos(2*pi*sampleIndex/temporalCellCount);
nFft = zeroPaddingFactor*temporalCellCount;

% 物理约定为 exp(i*k*x-i*omega*t)，所以使用 ifft 的正指数核后，
% e^{-i*omega*pT} 的峰会直接落在正的 omega 轴上。
windowedSamples = electricSamples.*reshape(window,[],1,1);
complexSpectrum = ifft(windowedSamples,nFft,1);
powerSpectrum = abs(complexSpectrum).^2;

% 两组基初态的功率相加，而不是先相干叠加复振幅；这样结果不依赖任意
% 的基矢相对相位，同时保留两条可被任一初态激发的 Floquet 分支。
powerSpectrum = reshape(sum(powerSpectrum,2),nFft,nK);

% 等相位采样率为每 T 一次，故所有 omega+m*Omega 自动混叠到第一时间
% Floquet 区。下式对奇数、偶数 nFft 都给出 [-1/2,1/2) 型有序频轴。
frequencyOrders = (-floor(nFft/2):ceil(nFft/2)-1).';
centeredIndices = mod(frequencyOrders,nFft)+1;
powerSpectrum = powerSpectrum(centeredIndices,:);
omegaNormalized = frequencyOrders/nFft;

% 每个 k 都由两组完整基初态独立激发，因此可以逐列归一化来突出峰位。
% 该操作会删除不同 k 之间的真实增长/响应强度差异，颜色不能解释成
% Im(omega) 或增长率；普通有限序列 FFT 也只能重建 Re(omega)。
columnMaximum = max(powerSpectrum,[],1);
if any(~isfinite(columnMaximum)) || any(columnMaximum <= 0)
    error(['FFT 功率发生溢出或出现空列；请减少 temporalCellCount，' ...
        '并检查材料与采样参数。']);
end
relativePower = powerSpectrum./columnMaximum;
displayFloor = 10^(-dynamicRangeDb/10);
if displayFloor == 0
    error('cfg.dynamicRangeDb 过大，已超出双精度 dB 显示范围。');
end
spectrumDb = 10*log10(max(relativePower,displayFloor));

%% ===================== 最小绘图输出 =====================

result.kNormalized = kScan/Omega;
result.omegaNormalized = omegaNormalized;
result.spectrumDb = spectrumDb;
result.rawColumnPower = columnMaximum;
result.dominantGrowthPerPeriod = dominantGrowthPerPeriod;
result.dominantImagOmegaNormalized = dominantGrowthPerPeriod/(2*pi);
result.dt = dt;
result.dx = dx;
result.courant = courant;
result.temporalPeriod = T;
result.stepsPerPeriod = stepsPerPeriod;
result.temporalCellCount = temporalCellCount;
result.nFft = nFft;
result.nativeOmegaResolutionNormalized = 1/temporalCellCount;
result.windowName = 'periodic-hann';
result.windowCoherentGain = mean(window);
result.windowEnergyGain = mean(window.^2);
if returnTemporalSamples
    result.electricSamples = electricSamples;
    result.sampleCellIndices = sampleIndex;
end
end

% -------------------------------------------------------------------------
function value = evaluate_epsilon(epsilonTime,timeValue)
% 材料函数必须保留输入尺寸，并在所有采样时刻返回有限正实介电常数。
value = epsilonTime(timeValue);
if ~isnumeric(value) || ~isequal(size(value),size(timeValue)) || ...
        ~isreal(value) || any(~isfinite(value)) || any(value <= 0)
    error('epsilonTime 必须返回与时间输入同尺寸的有限正实数组。');
end
end

% -------------------------------------------------------------------------
function value = read_positive_scalar(value,label)
% 统一读取周期、网格步长和动态范围等严格正实标量。
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value <= 0
    error('%s 必须是有限正实标量。',label);
end
end

% -------------------------------------------------------------------------
function value = read_integer(value,label,minimumValue)
% 统一读取时间步数、原胞数和零填充因子等有下界的整数。
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value ~= round(value) || value < minimumValue
    error('%s 必须是不小于 %d 的有限整数。',label,minimumValue);
end
end
