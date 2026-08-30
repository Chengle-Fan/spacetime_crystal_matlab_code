function result = fdtd_fft_bands(electricField,x,t,fftCfg)
%FDTD_FFT_BANDS  用有限样品复电场重建第一时间 Floquet 区的局部能带。
% 输入数组按“时间行、空间列”排列；输出只保留两个坐标轴和 dB 谱。

if nargin ~= 4 || ~isstruct(fftCfg) || ~isscalar(fftCfg)
    error('必须以 electricField、x、t 和标量 fftCfg 调用。');
end
required = {'temporalPeriod','spatialROI','timeROI'};
% temporalPeriod 决定 Floquet 折叠周期；两个 ROI 必须使用原始数组的
% MATLAB 闭区间下标。zeroPaddingFactor 和 dynamicRangeDb 为可选配置。
missing = required(~isfield(fftCfg,required));
if ~isempty(missing)
    error('fftCfg 缺少字段：%s。',strjoin(missing,', '));
end

if ~isnumeric(electricField) || ~ismatrix(electricField) || ...
        any(~isfinite(electricField),'all')
    error('electricField 必须是全部有限的二维数值数组。');
end
[nTimeFull,nSpaceFull] = size(electricField);
if ~isnumeric(x) || ~isvector(x) || numel(x) ~= nSpaceFull || ...
        ~isreal(x) || any(~isfinite(x)) || ~isnumeric(t) || ~isvector(t) || ...
        numel(t) ~= nTimeFull || ~isreal(t) || any(~isfinite(t))
    error('x、t 必须是与 electricField 尺寸匹配的有限实向量。');
end
x = x(:).';
t = t(:);
dx = uniform_step(x,'x');
dtRecord = uniform_step(t,'t');

T = fftCfg.temporalPeriod;
if ~isnumeric(T) || ~isscalar(T) || ~isreal(T) || ~isfinite(T) || T <= 0
    error('fftCfg.temporalPeriod 必须是有限正实标量。');
end
spatialROI = validate_roi(fftCfg.spatialROI,nSpaceFull,'spatialROI');
timeROI = validate_roi(fftCfg.timeROI,nTimeFull,'timeROI');

if isfield(fftCfg,'zeroPaddingFactor') && ~isempty(fftCfg.zeroPaddingFactor)
    zeroPaddingFactor = fftCfg.zeroPaddingFactor;
else
    zeroPaddingFactor = 1;
end
if ~isnumeric(zeroPaddingFactor) || ~isscalar(zeroPaddingFactor) || ...
        ~isreal(zeroPaddingFactor) || ~isfinite(zeroPaddingFactor) || ...
        zeroPaddingFactor < 1 || zeroPaddingFactor ~= round(zeroPaddingFactor)
    error('zeroPaddingFactor 必须是大于等于 1 的整数；小于 1 会截断时域数据。');
end
if isfield(fftCfg,'dynamicRangeDb') && ~isempty(fftCfg.dynamicRangeDb)
    dynamicRangeDb = fftCfg.dynamicRangeDb;
else
    dynamicRangeDb = 60;
end
if ~isnumeric(dynamicRangeDb) || ~isscalar(dynamicRangeDb) || ...
        ~isreal(dynamicRangeDb) || ~isfinite(dynamicRangeDb) || dynamicRangeDb <= 0
    error('dynamicRangeDb 必须是有限正实标量。');
end

% 必须先裁剪用户选择的 ROI，再检查整数周期。这样总仿真时间可以任意，
% 只要被选中的 FFT 时间窗本身包含整数个完整调制周期即可。空间 ROI 只
% 位于样品内部；纯时间晶体没有空间晶格，因此后续保留 signed k 而不折叠。
electricField = electricField(timeROI(1):timeROI(2), ...
    spatialROI(1):spatialROI(2));
[nTime,nSpace] = size(electricField);
if nTime < 4 || nSpace < 4
    error('FFT 的时间和空间 ROI 均至少需要 4 个采样点。');
end

samplesPerPeriodFloat = T/dtRecord;
if abs(samplesPerPeriodFloat-round(samplesPerPeriodFloat)) > 1e-9
    error('记录时间步与调制周期不可公度。');
end
samplesPerPeriod = round(samplesPerPeriodFloat);
if mod(nTime,samplesPerPeriod) ~= 0
    error('timeROI 必须恰好包含整数个周期的无重复端点采样。');
end
periodCount = nTime/samplesPerPeriod;

% 使用两端为零的对称 Hann 窗抑制有限 ROI 边缘泄漏。窗口会牺牲分辨率，
% 但不会像逐 k 归一化那样把弱杂散列人为抬升成假能带。
spaceWindow = symmetric_hann(nSpace).';
timeWindow = symmetric_hann(nTime);
windowedField = electricField.*timeWindow.*spaceWindow;

% 采用 E~exp(i*k*x-i*omega*t)：空间 fft 的峰位对应 +k，时间 ifft 的
% 峰位对应 +omega。若直接用 fft2 而不修正符号，频率轴会被反向解释。
% 时间零填充只在第一维进行，并不添加新的物理信息。
fieldK = fft(windowedField,[],2);
nFftTime = zeroPaddingFactor*nTime;
powerUnfolded = abs(ifft(fieldK,nFftTime,1)).^2;

% 原生频率间隔为 Omega/periodCount，零填充后为
% Omega/(zeroPaddingFactor*periodCount)，所以一个第一 Floquet 区恰含
% nFold=zeroPaddingFactor*periodCount 个频率 bin。相隔 Omega 的 bin 在
% 未移位序列中相隔 nFold；reshape 的第二维枚举 Nyquist 范围内的
% samplesPerPeriod 个 Floquet 谐波副本，对它们的功率求和完成折叠。
nFold = zeroPaddingFactor*periodCount;
powerFolded = reshape(powerUnfolded,nFold,samplesPerPeriod,nSpace);
powerFolded = sum(powerFolded,2);
powerFolded = reshape(powerFolded,nFold,nSpace);

frequencyOrders = -floor(nFold/2):ceil(nFold/2)-1;
frequencyPermutation = mod(frequencyOrders,nFold)+1;
powerFolded = powerFolded(frequencyPermutation,:);
omegaNormalized = frequencyOrders.'/nFold;

spaceOrders = -floor(nSpace/2):ceil(nSpace/2)-1;
spacePermutation = mod(spaceOrders,nSpace)+1;
powerFolded = powerFolded(:,spacePermutation);
Omega = 2*pi/T;
kNormalized = (spaceOrders*2*pi/(nSpace*dx))/Omega;

maximumPower = max(powerFolded,[],'all');
if ~isfinite(maximumPower) || maximumPower <= 0
    error('所选时空 ROI 内没有可用于 FFT 的非零场。');
end
relativePower = powerFolded/maximumPower;
minimumRelativePower = 10^(-dynamicRangeDb/10);
spectrumDb = 10*log10(max(relativePower,minimumRelativePower));

% 输出按全谱最大功率统一归一化，避免逐 k 归一化把弱杂散放大成伪能带。
% k 不折叠；omega 位于 [-Omega/2,Omega/2)。有限 ROI 的真实分辨率约为
% 2*pi/L_ROI 和 2*pi/t_ROI，Hann 窗还会进一步展宽。FFT 不能得到 Im(omega)。
result.kNormalized = kNormalized;
result.omegaNormalized = omegaNormalized;
result.spectrumDb = spectrumDb;
end

% -------------------------------------------------------------------------
function step = uniform_step(grid,label)
% 验证严格递增均匀网格，并返回平均步长。
differences = diff(grid);
if isempty(differences) || any(differences <= 0)
    error('%s 必须严格递增且至少包含两个点。',label);
end
step = mean(differences);
if max(abs(differences-step)) > 1e-10*max(1,abs(step))
    error('%s 必须是均匀网格。',label);
end
end

% -------------------------------------------------------------------------
function roi = validate_roi(value,maximumIndex,label)
% ROI 使用闭区间 MATLAB 下标，并至少包含两个点。
if ~isnumeric(value) || ~isvector(value) || numel(value) ~= 2 || ...
        any(~isfinite(value)) || any(value ~= round(value))
    error('fftCfg.%s 必须是两个有限整数下标。',label);
end
roi = value(:).';
if roi(1) < 1 || roi(2) > maximumIndex || roi(2) <= roi(1)
    error('fftCfg.%s 必须满足 1 <= 起点 < 终点 <= %d。',label,maximumIndex);
end
end

% -------------------------------------------------------------------------
function window = symmetric_hann(pointCount)
% Base MATLAB 实现的对称 Hann 窗，不依赖 Signal Processing Toolbox。
window = 0.5-0.5*cos(2*pi*(0:pointCount-1)'/(pointCount-1));
end
