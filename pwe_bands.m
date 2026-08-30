function result = pwe_bands(fourier, kScan)
%PWE_BANDS 用平面波展开法求光学时间晶体的复准频率能带。
%
% 对每个守恒的实波数 k，求解由电场和磁场时间谐波系数组成的广义
% 本征值问题：
%
%   [kI  -W        ] [E] = omega [0         I] [E]
%   [-W*Cepsilon kI] [H]         [Cepsilon 0] [H],
%
% 其中 W=diag(m*Omega)，Cepsilon 是介电常数的 Fourier 卷积矩阵，
% mu_r=1。由于时间周期性，同一个物理解可整体平移任意整数倍 Omega，
% 因而广义本征问题会返回相隔 Omega 的 Floquet 副本。程序从这些原始根
% 中选出两条物理分支的各一个代表，再把其实部折叠到标准第一时间
% Floquet 区间 [-Omega/2,Omega/2)。
%
% 复准频率采用场随时间 exp(-i*omega*t) 演化的约定。因此
% Im(omega)>0 表示指数增长，Im(omega)<0 表示指数衰减。该函数只负责
% 逐个 k 点求解和排序，不计算本征矢、模态权重，也不承诺输出行在跨越
% 简并点时具有连续的分支身份。
%
% 输入：
%   fourier : pwe_fourier 返回的 Fourier 数据结构；
%   kScan   : 非空的有限实波数向量，行向量和列向量均可。
%
% 输出：
%   result.k     : 1 x nK 的波数扫描行向量；
%   result.omega : 2 x nK 的复准频率，两行对应两个物理根；其实部位于
%                  第一 Floquet 区，虚部保持不变。

if nargin ~= 2 || ~isstruct(fourier) || ~isscalar(fourier)
    error('Use pwe_bands(fourier,kScan).');
end
required = {'T','Omega','Mtime','mList','Cepsilon'};
missing = required(~isfield(fourier,required));
if ~isempty(missing)
    error('fourier is missing field(s): %s.',strjoin(missing,', '));
end
if ~isnumeric(kScan) || isempty(kScan) || ~isvector(kScan) || ...
        ~isreal(kScan) || any(~isfinite(kScan))
    error('kScan must be a nonempty vector of finite real wavenumbers.');
end
kScan = kScan(:).';

T = fourier.T;
Omega = fourier.Omega;
Mtime = fourier.Mtime;
if ~isnumeric(T) || ~isscalar(T) || ~isreal(T) || ~isfinite(T) || T <= 0 || ...
        ~isnumeric(Omega) || ~isscalar(Omega) || ~isreal(Omega) || ...
        ~isfinite(Omega) || abs(Omega-2*pi/T) > 64*eps(max(1,Omega)) || ...
        ~isnumeric(Mtime) || ~isscalar(Mtime) || Mtime < 0 || ...
        Mtime ~= round(Mtime)
    error('The Fourier geometry is inconsistent; rerun pwe_fourier.');
end

mList = (-Mtime:Mtime).';
S = numel(mList);
Cepsilon = fourier.Cepsilon;
if ~isequal(fourier.mList,mList) || ~isnumeric(Cepsilon) || ...
        ~isequal(size(Cepsilon),[S S]) || any(~isfinite(Cepsilon),'all')
    error('The Fourier harmonic ordering or convolution matrix is invalid.');
end

I = eye(S);
Z = zeros(S);
% W 的对角元 m*Omega 表示各时间谐波相对于准频率的频移。
W = diag(mList*Omega);
% 广义本征问题 A*v=omega*B*v 中与 k 无关的右端矩阵。
B = [Z,I;Cepsilon,Z];
omega = complex(zeros(2,numel(kScan)));

for ik = 1:numel(kScan)
    % 空间均匀的时间晶体保持波数 k 守恒，因此每个 k 可以独立求解。
    K = kScan(ik)*I;
    A = [K,-W;-W*Cepsilon,K];
    raw = eig(A,B);

    % Floquet 副本之间相差整数倍 Omega。选取实部最接近 Omega/4 的两个根，
    % 等价于优先选取平移代表区间 [-Omega/4,3*Omega/4) 内的两条物理根。
    % 这个区间的边界避开了通常出现动量带隙的 Floquet 区中心和区边界，
    % 能减少物理根恰好落在代表区边界时的跳选。随后再统一折回标准第一
    % Floquet 区。整个选择过程不使用本征矢权重，也不借助 TMM 等外部解。
    [~,order] = sort(abs(real(raw)-Omega/4));
    selected = raw(order(1:2));
    omega(:,ik) = fold_frequency(selected,Omega);
end

result.k = kScan;
result.omega = omega;
end

% -------------------------------------------------------------------------
% 只对准频率实部进行模 Omega 折叠；虚部代表真实增长/衰减率，不能折叠。
function folded = fold_frequency(omega,Omega)
folded = mod(real(omega)+Omega/2,Omega)-Omega/2 + 1i*imag(omega);
end
