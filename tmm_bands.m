function result = tmm_bands(kScan, epsLayers, muLayers, durations)
%TMM_BANDS 用时间传输矩阵法求方波时间晶体的复准频率能带。
%
% 物理模型由两个在时间上依次出现的均匀层组成。每层在空间上均匀，且
% 介电常数和磁导率均为严格正的实数，因此模型是无损、无色散的普通
% 介质。时间界面被视为瞬时切换；由 Maxwell 方程的时间边界条件，
% 电位移 D 和磁感应强度 B 在界面两侧连续，所以选取 [D;B] 作为状态量
% 不需要额外插入界面矩阵。
%
% 对空间因子 exp(i*k*x)，在任意一个常参数时间层内有
%
%   d/dt [D;B] = -i*k*[0 1/mu; 1/eps 0]*[D;B].
%
% 对该常系数方程做解析指数运算可得到精确的 2x2 层传播矩阵，因此 TMM
% 本身没有时间步长离散误差。按时间先后将两个层矩阵左乘，得到一个周期
% 的演化矩阵 U。其本征值 lambda 是 Floquet 乘子，并满足
%
%   lambda = exp(-i*omega*T),    omega = i*log(lambda)/T.
%
% 这里使用复对数主值，使 Re(omega) 位于第一时间 Floquet 区。采用
% exp(-i*omega*t) 约定时，Im(omega)>0 表示场随时间指数增长，
% Im(omega)<0 表示指数衰减；在无动量带隙的通带内两者应为零。
%
% 输入：
%   kScan     : 非空的有限实波数向量，行向量和列向量均可；
%   epsLayers : [eps1 eps2]，两个时间层的相对介电常数；
%   muLayers  : 正标量或 [mu1 mu2]；标量会自动用于两个时间层；
%   durations : [dt1 dt2]，两个时间层严格为正的持续时间，周期
%               T=dt1+dt2。
%
% 输出：
%   result.k     : 1 x nK 的波数扫描行向量；
%   result.omega : 2 x nK 的复准频率。每个 k 点的两根按
%                  [Re(omega),Im(omega)] 独立排序；不承诺跨 k 点连续追踪
%                  同一分支。

if nargin ~= 4
    error('Use tmm_bands(kScan,epsLayers,muLayers,durations).');
end
if ~isnumeric(kScan) || isempty(kScan) || ~isvector(kScan) || ...
        ~isreal(kScan) || any(~isfinite(kScan))
    error('kScan must be a nonempty vector of finite real wavenumbers.');
end
if ~isnumeric(epsLayers) || ~isvector(epsLayers) || numel(epsLayers) ~= 2 || ...
        ~isreal(epsLayers) || any(~isfinite(epsLayers)) || any(epsLayers <= 0)
    error('epsLayers must contain two finite, real, strictly positive values.');
end
if ~isnumeric(muLayers) || ~isvector(muLayers) || ...
        ~(isscalar(muLayers) || numel(muLayers) == 2) || ...
        ~isreal(muLayers) || any(~isfinite(muLayers)) || any(muLayers <= 0)
    error('muLayers must be one or two finite, real, strictly positive values.');
end
if ~isnumeric(durations) || ~isvector(durations) || numel(durations) ~= 2 || ...
        ~isreal(durations) || any(~isfinite(durations)) || any(durations <= 0)
    error('durations must contain two finite, real, strictly positive values.');
end

kScan = kScan(:).';
epsLayers = epsLayers(:).';
durations = durations(:).';
% 允许用户用一个 mu 值表示非磁性或两层磁导率相同的情况。
if isscalar(muLayers)
    muLayers = [muLayers muLayers];
else
    muLayers = muLayers(:).';
end
T = sum(durations);
omega = complex(zeros(2,numel(kScan)));

for ik = 1:numel(kScan)
    % U 初始为单位矩阵。由于第一个时间层先作用，传播矩阵必须写成
    % U=P*U；循环结束后得到 U=P2*P1。
    U = eye(2);
    for layer = 1:2
        epsr = epsLayers(layer);
        mur = muLayers(layer);
        theta = kScan(ik)*durations(layer)/sqrt(epsr*mur);
        c = cos(theta);
        s = sin(theta);

        % theta 是该层内积累的相位。下面是 [D;B] 状态的解析传播矩阵，
        % 直接来自常系数 Maxwell 生成元的矩阵指数。
        P = [c, -1i*s*sqrt(epsr/mur); ...
             -1i*s*sqrt(mur/epsr), c];
        U = P*U;
    end

    % 2x2 特征方程为 lambda^2-tr(U)*lambda+det(U)=0。用半迹书写后，
    % 两个候选根为 halfTrace +/- sqrt(halfTrace^2-determinant)。
    halfTrace = (U(1,1)+U(2,2))/2;
    determinant = U(1,1)*U(2,2)-U(1,2)*U(2,1);
    root = sqrt(halfTrace^2-determinant);
    candidates = [halfTrace+root, halfTrace-root];
    % 强增长带隙中两个 Floquet 乘子的模可能相差很多。直接用较小候选根
    % 会发生严重相消，因此先选模较大的稳定根，再由乘积
    % lambdaLarge*lambdaSmall=det(U) 恢复较小根。
    [~,largeIndex] = max(abs(candidates));
    lambdaLarge = candidates(largeIndex);
    if ~isfinite(lambdaLarge) || lambdaLarge == 0 || ...
            ~isfinite(determinant) || determinant == 0
        error('The one-period evolution is singular or non-finite at k=%.16g.', ...
            kScan(ik));
    end
    lambdaSmall = determinant/lambdaLarge;

    % 主对数自动把准频率实部放入第一 Floquet 区；随后只在当前 k 点
    % 按实部、虚部排序，使输出顺序确定且可重复。
    pair = 1i*log([lambdaLarge;lambdaSmall])/T;
    [~,order] = sortrows([real(pair),imag(pair)],[1 2]);
    omega(:,ik) = pair(order);
end

result.k = kScan;
result.omega = omega;
end
