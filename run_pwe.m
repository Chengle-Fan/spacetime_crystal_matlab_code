%RUN_PWE  用 PWE 计算一维光学时间晶体的复准频率能带。
%
% 通常只需修改下方“用户参数”区域。脚本支持方波和正弦两种时间调制，
% 计算两条复准频率能带，并分别绘制其实部和虚部。采用归一化单位
% c0=eps0=mu0=1，且固定相对磁导率 mu_r=1；因此 k 与频率具有相同量纲。
% 准频率约定为 exp(-i*omega*t)，所以正虚部表示增长，负虚部表示衰减。

%% ===================== 用户参数 =====================

% 时间调制类型：'square' 表示两值方波，'sinusoidal' 表示正弦调制。
modulationType = 'square';

% 两种调制共用的介电常数上下界，必须满足 epsHigh >= epsLow > 0。
epsHigh = 5;
epsLow  = 1;

% 方波参数：每个周期先保持 epsHigh，随后切换到 epsLow。
dutyCycle = 0.5;       % epsHigh 在一个周期内所占的时间比例，范围为 (0,1)
timeShift = 0;         % 整个方波相对于 t=0 的时间平移，单位与 T 相同

% 正弦参数：epsilon(t)=平均值+振幅*cos(Omega*t+sinePhase)。
sinePhase = 0;         % 初始相位，单位为弧度

% 时间周期、Fourier 截断阶数和材料采样数。
% Mtime 越大，保留的时间谐波越多；Nt 用于计算材料 Fourier 系数。
% 方波不连续，通常需要比正弦调制更大的 Mtime 和 Nt，并应通过同时增大
% 两者来检查能带收敛，而不能只依赖某一组截断参数。
T = 1;
Mtime = 19;            % 保留 -Mtime:Mtime 共 2*Mtime+1 个时间谐波
Nt = 4096;             % 一个周期内的无重复端点材料采样数

% 波数扫描使用 k/Omega 归一化，其中 Omega=2*pi/T。
kMinNormalized = 0;
kMaxNormalized =  2;
nK = 301;

% 实部绘图纵轴范围，决定完整 Floquet 本征谱中哪些副本可见。
% Floquet 副本只沿 Re(omega) 相隔整数倍 Omega。
realOmegaYLim = [-0.5 0.5];

%% ===================== 构造时间调制材料 =====================

Omega = 2*pi/T;
switch lower(modulationType)
    case 'square'
        if dutyCycle <= 0 || dutyCycle >= 1
            error('dutyCycle must be strictly between 0 and 1.');
        end
        epsFun = @(t) epsLow + (epsHigh-epsLow) .* ...
            (mod(t-timeShift,T) < dutyCycle*T);
        profileName = '方波调制';

    case 'sinusoidal'
        % 用上下界确定平均值和振幅，从而保证介电常数范围仍为
        % [epsLow,epsHigh]。
        epsMean = (epsHigh+epsLow)/2;
        epsAmplitude = (epsHigh-epsLow)/2;
        epsFun = @(t) epsMean + epsAmplitude*cos(Omega*t+sinePhase);
        profileName = '正弦调制';

    otherwise
        error('modulationType must be ''square'' or ''sinusoidal''.');
end

if epsLow <= 0 || epsHigh < epsLow
    error('Require epsHigh >= epsLow > 0.');
end

%% ===================== PWE 能带计算 =====================

% 第一步把 epsilon(t) 转换为 Fourier 卷积矩阵；第二步对每个守恒波数 k
% 求解广义本征值问题。bands.omega 包含每个 k 点的完整 Floquet 本征谱。
pweCfg = struct('T',T,'Mtime',Mtime,'Nt',Nt);%Cfg struct store data
fourier = pwe_fourier(epsFun,pweCfg);
kNormalized = linspace(kMinNormalized,kMaxNormalized,nK);
kScan = kNormalized*Omega;
bands = pwe_bands(fourier,kScan);

%% ===================== 绘制能带实部和虚部 =====================

% 每个波数对应完整的截断本征谱，因此复制横坐标以便统一绘制散点。
kPlot = repmat(kNormalized,size(bands.omega,1),1);

% 虚部沿 Floquet 副本不具周期性，不能用 ylim 选择。保持原来的选根规则：
% 每个 k 取实部最接近 Omega/4 的两根，并只绘制它们的虚部。
imagOmega = complex(zeros(2,nK));
for ik = 1:nK
    [~,order] = sort(abs(real(bands.omega(:,ik))-Omega/4));
    imagOmega(:,ik) = bands.omega(order(1:2),ik);
end
kImagPlot = repmat(kNormalized,2,1);

figure('Color','w','Position',[100 100 1050 440]);
layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

axReal = nexttile(layout);
scatter(axReal,kPlot(:),real(bands.omega(:))/Omega,12,'filled');
xlabel(axReal,'k/\Omega');
ylabel(axReal,'Re(\omega)/\Omega');
title(axReal,'准频率实部');
ylim(axReal,realOmegaYLim);
grid(axReal,'on');
box(axReal,'on');

axImag = nexttile(layout);
scatter(axImag,kImagPlot(:),imag(imagOmega(:))/Omega,12,'filled');
xlabel(axImag,'k/\Omega');
ylabel(axImag,'Im(\omega)/\Omega');
title(axImag,'准频率虚部');
grid(axImag,'on');
box(axImag,'on');

title(layout,profileName);
