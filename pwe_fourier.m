function fourier = pwe_fourier(epsFun, pweCfg)
%PWE_FOURIER 构造光学时间晶体的介电常数 Fourier 卷积矩阵。
%
% 物理模型限定为空间均匀、时间周期、无损且无色散的标量介质，并取
% 相对磁导率 mu_r=1。这里采用的时间 Fourier 展开约定为
%
%   epsilon(t) = sum_m epsilon_m exp(-i*m*Omega*t),  Omega = 2*pi/T.
%
% 因此第 m 阶系数由
%
%   epsilon_m = (1/T) integral_0^T epsilon(t) exp(+i*m*Omega*t) dt
%
% 给出。程序在无重复端点的均匀网格 t=(0:Nt-1)*T/Nt 上用周期求和
% 近似该积分。不能同时采样 t=0 和 t=T，因为它们是同一个周期点，
% 重复端点会破坏离散周期平均并污染 Fourier 系数。
%
% PWE 保留 m=-Mtime:Mtime 共 2*Mtime+1 个场谐波。卷积矩阵元素依赖
% 两个场谐波的差阶 m-m'，所以必须计算 -2*Mtime:2*Mtime 的材料系数。
% 为避免这些差阶在离散 Fourier 网格上发生混叠，强制要求
% Nt >= 4*Mtime+1。方波在时间界面处不连续，其 Fourier 系数衰减较慢，
% 实际计算时通常应令 Nt 远大于这个最低要求，并同时检查 Mtime 与 Nt
% 增大后能带是否收敛。
%
% 输入：
%   epsFun : 介电常数函数句柄 epsFun(t)，返回值尺寸必须与 t 完全相同；
%            每个采样值都必须是有限、实数且严格为正。
%   pweCfg : 标量结构体，包含
%              T      : 严格为正的时间调制周期；
%              Mtime  : 非负整数，时间谐波截断阶数；
%              Nt     : 正整数，材料采样点数，且 Nt >= 4*Mtime+1。
%
% 输出：
%   fourier.T、fourier.Omega、fourier.Mtime : 时间周期、调制角频率和截断阶；
%   fourier.mList    : 列向量形式的谐波阶数 -Mtime:Mtime；
%   fourier.Cepsilon : 介电常数 Fourier 卷积矩阵，行列均按 mList 排列。

if nargin ~= 2 || ~isa(epsFun,'function_handle') || ...
        ~isstruct(pweCfg) || ~isscalar(pweCfg)
    error('Use pwe_fourier(epsFun,pweCfg) with a scalar configuration struct.');
end

T = positive_scalar(pweCfg,'T');
Mtime = nonnegative_integer(pweCfg,'Mtime');
Nt = positive_integer(pweCfg,'Nt');
minNt = 4*Mtime + 1;
if Nt < minNt
    error('pweCfg.Nt must satisfy Nt >= 4*Mtime+1 = %d.',minNt);
end

Omega = 2*pi/T;
% 无重复端点的一个完整时间周期；最后一点是 T-T/Nt，而不是 T。
t = (0:Nt-1)*T/Nt;
epsilon = epsFun(t);
if ~isnumeric(epsilon) || ~isequal(size(epsilon),size(t)) || ...
        ~isreal(epsilon) || any(~isfinite(epsilon)) || any(epsilon <= 0)
    error(['epsFun(t) must return finite, real, strictly positive values ', ...
        'with exactly the same size as t.']);
end

% 按上面的符号约定直接计算材料 Fourier 系数。由于后续卷积需要
% epsilon_(m-m')，这里覆盖场谐波之间可能出现的全部差阶。
mDifference = -2*Mtime:2*Mtime;
epsilonCoeff = complex(zeros(size(mDifference)));
for im = 1:numel(mDifference)
    m = mDifference(im);
    epsilonCoeff(im) = mean(epsilon .* exp(1i*m*Omega*t));
end

mList = (-Mtime:Mtime).';
dm = mList - mList.';
% dm+2*Mtime+1 把物理差阶 [-2M,2M] 映射为 MATLAB 的正整数下标。
Cepsilon = epsilonCoeff(dm + 2*Mtime + 1);

% 对严格为正的实函数 epsilon(t)，任意有限 Fourier 子空间中的乘法算符
% 都应对应 Hermitian 正定矩阵。Hermitian 检查可以尽早发现 Fourier 符号
% 或差阶下标写反的问题；Cholesky 检查则防止采样不足或材料函数异常产生
% 非物理的非正定介电矩阵。先检查，再做对称化以清除舍入级反 Hermitian 量。
scale = max(1,norm(Cepsilon,'fro'));
if norm(Cepsilon-Cepsilon','fro')/scale > 1e-11
    error('The permittivity convolution matrix is not Hermitian.');
end
Cepsilon = (Cepsilon+Cepsilon')/2;
[~,cholFlag] = chol(Cepsilon);
if cholFlag ~= 0
    error(['The permittivity convolution matrix is not positive definite; ', ...
        'increase Nt or inspect epsilon(t).']);
end

fourier.T = T;
fourier.Omega = Omega;
fourier.Mtime = Mtime;
fourier.mList = mList;
fourier.Cepsilon = Cepsilon;
end

% -------------------------------------------------------------------------
% 读取并验证严格为正的有限实标量配置项。
function value = positive_scalar(cfg,name)
if ~isfield(cfg,name)
    error('pweCfg.%s is required.',name);
end
value = cfg.(name);
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value <= 0
    error('pweCfg.%s must be a positive finite real scalar.',name);
end
end

% -------------------------------------------------------------------------
% 读取并验证非负有限整数配置项。
function value = nonnegative_integer(cfg,name)
if ~isfield(cfg,name)
    error('pweCfg.%s is required.',name);
end
value = cfg.(name);
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value < 0 || value ~= round(value)
    error('pweCfg.%s must be a nonnegative finite integer.',name);
end
end

% -------------------------------------------------------------------------
% 正整数是在非负整数检查基础上进一步排除零。
function value = positive_integer(cfg,name)
value = nonnegative_integer(cfg,name);
if value < 1
    error('pweCfg.%s must be a positive integer.',name);
end
end
