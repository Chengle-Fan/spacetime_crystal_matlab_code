function out = fdtd1d_db(cfg)
%FDTD1D_DB  一维 D/B-Yee 时域有限差分求解器，适用于时空变化介质。
%

%   ▏  概述

%
%   本函数实现一维 Yee 网格上的 leapfrog（蛙跳）时间推进方案。与传统
%   E/H FDTD 不同，这里直接更新的是电位移场 D 和磁感应场 B，而非电场 E
%   和磁场 H。这是因为：
%
%     1. 在时间界面（介电常数/磁导率发生突变的时刻），D 和 B 自动满足
%        正确的连续性条件 D(t+) = D(t-) 和 B(t+) = B(t-)。
%     2. 在每个时间步通过本构关系 E = D/ε(x,t) 和 H = B/μ(x,t) 恢复
%        电场和磁场，可以自然地处理任意时变材料参数。
%

%   ▏  Yee 网格的空间排布

%
%   在一维情况下，假设电磁波沿 x 方向传播，E 沿 y 方向偏振，H 沿 z 方向
%   偏振。Yee 网格将 E 和 H 在空间上错开半个步长 dx/2：
%
%       E 网格（整数格点）:  x₀    x₁    x₂    x₃    ...   x_{Nx-1}
%                            |     |     |     |     |
%       H 网格（半整数格点）:   xH₀   xH₁   xH₂   ...  xH_{Nx-2}
%
%    其中 xH(i) = x(i) + dx/2，即 H 位于相邻两个 E 格点的正中间。
%    差分格式 ∂E/∂x ≈ [E(i+1)-E(i)]/dx 在 H 格点 xH(i) 处是二阶精确的
%    中心差分。同样，∂H/∂x 在 E 格点处也是中心差分。
%
%   周期边界时，最后一个 E 格点也用于计算环绕的 ∂E/∂x，因此 H 网格长度
%    与 E 相同（Nx 个）。吸收边界时，两端的 ∂E/∂x 可以省略，故 H 网格
%    少一个点（Nx-1 个）。
%

%   ▏  Yee 网格的时间排布（leapfrog 蛙跳格式）

%
%   时间上同样采用交错排布，相隔半个时间步 dt/2：
%
%       时间轴: t=-dt/2    t=0      t=dt/2    t=dt     t=3dt/2   t=2dt
%                  |         |         |         |         |         |
%       已知量:    B⁰        E⁰        B¹        E¹        B²        E²
%                  |         |         |         |         |         |
%       本构关系:  H⁰=B⁰/μ   D⁰=ε·E⁰   H¹=B¹/μ   D¹=ε·E¹   H²=B²/μ
%
%   因此初始条件必须分别在两个时刻给出：
%     - cfg.E0:     电场在 t = 0 时刻，定义在 E 网格上
%     - cfg.Hhalf0: 磁场在 t = -dt/2 时刻，定义在 H 网格上
%
%   如果初始场是行波脉冲，Hhalf0 需要根据传播速度和 dt/2 的时移来构造。
%
%   对于与 E 整数时间节点对齐的介电常数突变，可通过
%   cfg.temporalInterfaces 显式给出界面时刻。在以界面为中心、从
%   B^{n-1/2} 推进到 B^{n+1/2} 的更新中，求解器使用
%
%       E_interface = 0.5·[D/ε(t^-) + D/ε(t^+)]
%
%   计算空间旋度，避免只使用界面一侧 E 所造成的一阶时间误差。
%

%   ▏  Leapfrog 更新方程

%
%   法拉第定律 ∂B/∂t = -∂E/∂x 的离散（步骤 1）:
%
%       B^{n+1/2}(i) = B^{n-1/2}(i) - (Δt/Δx)·[E^n(i+1) - E^n(i)]
%
%   本构关系恢复 H:
%
%       H^{n+1/2}(i) = B^{n+1/2}(i) / μ(xH(i), t^{n+1/2})
%
%   安培定律 ∂D/∂t = -∂H/∂x 的离散（步骤 2）:
%
%       D^{n+1}(i) = D^n(i) - (Δt/Δx)·[H^{n+1/2}(i) - H^{n+1/2}(i-1)]
%
%   本构关系恢复 E:
%
%       E^{n+1}(i) = D^{n+1}(i) / ε(x(i), t^{n+1})
%
%   注意 Δt/Δx 是无量纲数，与 Courant 条件直接相关。

%   ▏  CFL（Courant-Friedrichs-Lewy）稳定性条件
%
%   一维 Yee FDTD 的稳定性要求 Courant 数严格小于 1：
%
%       S = (v_max · Δt) / Δx < 1
%
%   其中 v_max 是介质中的最大相速度 v_max = max[1/√(ε_r·μ_r)]。
%   对于真空（ε_r = μ_r = 1），即要求 Δt < Δx/c₀。
%
%   本代码通过"采样 CFL 审计"来检查稳定性：
%     - 在仿真时间范围内均匀采样若干个时刻
%     - 在每个采样时刻，遍历整个空间网格计算 1/√(ε_r·μ_r)
%     - 取所有采样中的最大值作为 v_max
%     - 若 Courant 数 ≥ 1，直接报错退出
%
%   对于色散介质、负折射率介质或有源介质，标量相速度的 CFL 条件不足以
%   保证稳定性，此时需要通过 cfg.maxWaveSpeed 手动指定最大波速，或
%   使用更复杂的稳定性分析。
%

%   ▏  Sponge 吸收边界条件

%
%   计算域两端各留 spongeCells 个格点作为海绵吸收层。每个时间步将
%   靠近边界的场乘以衰减因子:
%
%       damp = exp(-strength · s³)
%
%   其中 s = (spongeCells - 到吸收层内边界的距离) / spongeCells。
%   三次方剖面从内部（s=0, damp=1, 无衰减）平滑过渡到边界
%   （s=1, damp=exp(-strength), 最大衰减），减小界面反射。
%   注意：这不是严格的 CPML（卷积完美匹配层），而是简化的乘法吸收。

%   ▏  输入参数

%
%   必填字段 (cfg):
%     x              - E/D 场的均匀空间网格（行向量或列向量），至少 5 点
%                      建议每波长至少 32 个网格点
%     dt             - 时间步长，必须满足 CFL 条件 Δt < Δx/v_max
%                      建议值: 0.8·dx/c₀（均匀介质）或更保守
%     nSteps         - 总时间步数，总仿真时间 = nSteps·dt
%     epsFun(x,t)    - 相对介电常数函数句柄
%                      输入: x（空间坐标数组）和 t（标量时间）
%                      输出: ε_r 数组，长度与 x 相同
%
%   可选字段 (cfg):
%     muFun(x,t)       - 相对磁导率函数句柄，默认 @(x,t) ones(size(x))
%     E0               - 初始电场（t=0, E 网格），默认全零
%     Hhalf0           - 初始磁场（t=-dt/2, H 网格），默认全零
%     boundary         - 边界条件: 'sponge'（默认）或 'periodic'
%     spongeCells      - 吸收层宽度（格点数），默认 min(80, floor(Nx/8))
%     spongeStrength   - 吸收强度，默认 0.12
%     recordEvery      - 场记录间隔（每 N 步记录一次），默认 1（每步记录）
%     storeFields      - 是否存储完整 E/H 历史，默认 true
%                        设为 false 可大幅降低内存占用
%     storeD           - 是否独立存储 E 网格上的 D 历史，默认 false
%                        与 storeFields 相互独立；即使 storeFields=false
%                        也可以令 storeD=true 保存 D
%     temporalInterfaces - 介电常数突变时刻的严格递增实数向量，默认空
%                        每个时刻必须在 temporalInterfaceTolerance 内与
%                        E 的整数时间节点 n·dt 对齐，并且小于 nSteps·dt
%     temporalInterfaceTolerance - temporalInterfaces 的绝对对齐容差；
%                        默认由 dt 和总仿真时长的浮点精度确定，且小于 dt/4
%     spectralFilterMask - 可选的周期网格 FFT 掩码（长度 Nx，默认空）
%                        若给出，则在每个 temporalInterfaces 更新后同时
%                        过滤 D 和 B。它用于解析信号/窄带脉冲，避免有源
%                        时变介质把舍入噪声中的非物理高 k 分量指数放大
%     probeIndices     - 探针点索引（E 网格），即使 storeFields=false 也会
%                        记录这些点的场值
%     stabilityTimes   - CFL 审计的采样时刻数组
%     maxWaveSpeed     - 手动指定的最大波速，跳过采样审计
%     sourceD(x,t,n)   - D 场增量源函数句柄，在每个时间步更新 D 后调用
%     precision       - 数值精度: 'double'（默认）或 'single'
%                       'single' 可节省约 50% 内存，适合大规模/长时仿真
%     progressBar     - 是否显示进度条: true 或 false（默认）
%                       长仿真时建议开启，方便监控推进进度

%   ▏  输出结构体
%
%   out.x               - E 空间网格
%   out.xH              - H 空间网格（相对于 E 偏移 dx/2）
%   out.t               - 记录时刻数组
%   out.E               - 电场历史矩阵 (nRecords × Nx)
%   out.H               - 磁场历史矩阵 (nRecords × Nx)，已插值到 E 网格
%   out.D               - 电位移历史矩阵 (nRecords × Nx)；storeD=false 时为空
%   out.probeIndices    - 探针索引
%   out.probeX          - 探针位置
%   out.probeE          - 探针处电场记录
%   out.probeH          - 探针处磁场记录
%   out.energy          - 瞬时电磁能量时间序列
%   out.finalD/B/E/H    - 最终时刻的场
%   out.dx              - 空间步长
%   out.dt              - 时间步长
%   out.boundary        - 使用的边界条件类型
%   out.sampledMaxWaveSpeed - CFL 审计得到的最大波速
%   out.sampledCourant  - CFL 审计的 Courant 数
%   out.processedTemporalInterfaceCount - 主循环实际处理的时间界面数

%==========================================================================
% 第 1 部分: 空间网格验证
%==========================================================================

% 将 x 转为行向量，方便后续处理
x = cfg.x(:).';
Nx = numel(x);                  % E 网格总点数
if Nx < 5
    error('The E grid must contain at least five points.');
end

% 验证网格的均匀性 —— Yee 算法要求均匀网格
% 计算相邻格点间距的向量 diff(x)，检查所有间距是否一致
dxVector = diff(x);
dx = mean(dxVector);            % 空间步长 Δx
if max(abs(dxVector-dx)) > 1e-10*max(1,abs(dx))
    error('cfg.x must be uniformly spaced.');
end

% 提取时间步长和总步数
dt = cfg.dt;                    % 时间步长 Δt
nSteps = cfg.nSteps;            % 总推进步数 N_t
% 总仿真时间 = nSteps · dt
if ~isscalar(dt) || ~isreal(dt) || ~isfinite(dt) || dt <= 0
    error('cfg.dt must be a finite positive real scalar.');
end
if ~isscalar(nSteps) || ~isreal(nSteps) || ~isfinite(nSteps) || ...
        nSteps < 0 || nSteps ~= round(nSteps)
    error('cfg.nSteps must be a nonnegative integer scalar.');
end

%==========================================================================
% 第 2 部分: 可选参数的默认值处理
%==========================================================================

% --- 磁导率函数 ---
% 若未指定 μ(x,t)，默认真空磁导率 μ_r = 1
if ~isfield(cfg,'muFun') || isempty(cfg.muFun)
    muFun = @(xq,tq) ones(size(xq));
else
    muFun = cfg.muFun;
end

% --- 边界条件 ---
% 默认使用海绵吸收层（sponge），也可选择周期边界（periodic）
if ~isfield(cfg,'boundary') || isempty(cfg.boundary)
    boundary = 'sponge';
else
    boundary = lower(cfg.boundary);
end

% --- 记录间隔 ---
% recordEvery=N 表示每 N 步记录一次场和能量，减少输出数据量
if ~isfield(cfg,'recordEvery') || isempty(cfg.recordEvery)
    recordEvery = 1;            % 默认每步都记录
else
    recordEvery = cfg.recordEvery;
end

% --- 场存储开关 ---
% storeFields=false 时只记录探针点，不存储全域场，大幅降低内存占用
if ~isfield(cfg,'storeFields') || isempty(cfg.storeFields)
    storeFields = true;
else
    storeFields = logical(cfg.storeFields);
end

% --- D 场历史存储开关 ---
% 与 storeFields 独立，以便只保存复现图真正需要的 D 场。
if ~isfield(cfg,'storeD') || isempty(cfg.storeD)
    storeD = false;
else
    storeDValue = cfg.storeD;
    if ~(islogical(storeDValue) || isnumeric(storeDValue)) || ...
            ~isscalar(storeDValue) || ~isreal(storeDValue) || ...
            ~isfinite(double(storeDValue))
        error('cfg.storeD must be a finite real logical or numeric scalar.');
    end
    storeD = logical(storeDValue);
end

% --- 探针点 ---
% 即使在低内存模式（storeFields=false）下，探针点的场也会被记录
% 探针索引必须是 1..Nx 范围内的整数
if ~isfield(cfg,'probeIndices') || isempty(cfg.probeIndices)
    probeIndices = zeros(1,0);
else
    probeIndices = unique(cfg.probeIndices(:).');
    if any(probeIndices < 1) || any(probeIndices > Nx) || ...
            any(probeIndices ~= round(probeIndices))
        error('cfg.probeIndices must contain valid integer E-grid indices.');
    end
end

% --- 数值精度 ---
% 默认双精度；可选 'single' 以节省内存并可能加速
if ~isfield(cfg,'precision') || isempty(cfg.precision)
    precision = 'double';
else
    precision = lower(cfg.precision);
    if ~ismember(precision, {'double','single'})
        error('cfg.precision must be ''double'' or ''single''.');
    end
end
useSingle = strcmp(precision,'single');

% --- 进度条 ---
% 长仿真时显示进度条，方便监控推进进度
if ~isfield(cfg,'progressBar') || isempty(cfg.progressBar)
    showProgress = false;
else
    showProgress = logical(cfg.progressBar);
end

% --- 显式介电时间界面 ---
% 界面位于 E/D 的整数时间节点。循环第 step 次更新 B 时，其中心时刻为
% (step-1)*dt，因此节点 n 对应循环索引 n+1。
timeScale = max(1,abs(nSteps*dt));
defaultInterfaceTolerance = max(128*eps(timeScale),1e-12*dt);
defaultInterfaceTolerance = min(defaultInterfaceTolerance,dt/100);
if isfield(cfg,'temporalInterfaceTolerance') && ...
        ~isempty(cfg.temporalInterfaceTolerance)
    temporalInterfaceTolerance = cfg.temporalInterfaceTolerance;
    if ~isscalar(temporalInterfaceTolerance) || ...
            ~isreal(temporalInterfaceTolerance) || ...
            ~isfinite(temporalInterfaceTolerance) || ...
            temporalInterfaceTolerance <= 0 || ...
            temporalInterfaceTolerance >= dt/4
        error(['cfg.temporalInterfaceTolerance must be a finite positive ', ...
            'real scalar smaller than cfg.dt/4.']);
    end
else
    temporalInterfaceTolerance = defaultInterfaceTolerance;
end

if ~isfield(cfg,'temporalInterfaces') || isempty(cfg.temporalInterfaces)
    temporalInterfaces = zeros(1,0);
    temporalInterfaceNodes = zeros(1,0);
else
    temporalInterfacesInput = cfg.temporalInterfaces;
    if ~isnumeric(temporalInterfacesInput) || ...
            ~isvector(temporalInterfacesInput) || ...
            ~isreal(temporalInterfacesInput) || ...
            any(~isfinite(temporalInterfacesInput(:)))
        error('cfg.temporalInterfaces must be a finite real numeric vector.');
    end
    temporalInterfacesInput = double(temporalInterfacesInput(:).');
    if any(diff(temporalInterfacesInput) <= 0)
        error('cfg.temporalInterfaces must be strictly increasing.');
    end

    temporalInterfaceNodes = round(temporalInterfacesInput/dt);
    temporalInterfaces = temporalInterfaceNodes*dt;
    alignmentError = abs(temporalInterfacesInput-temporalInterfaces);
    if any(alignmentError > temporalInterfaceTolerance)
        error(['Each cfg.temporalInterfaces value must align with an E-grid ', ...
            'time n*cfg.dt within cfg.temporalInterfaceTolerance. ', ...
            'Maximum alignment error is %.6g.'],max(alignmentError));
    end
    if any(temporalInterfaceNodes < 0) || ...
            any(temporalInterfaceNodes >= nSteps)
        error(['cfg.temporalInterfaces must lie in [0, cfg.nSteps*cfg.dt). ', ...
            'An interface at the final time has no following B half-step.']);
    end
    if any(diff(temporalInterfaceNodes) <= 0)
        error(['Distinct cfg.temporalInterfaces values must map to distinct ', ...
            'E-grid time nodes.']);
    end
end

temporalInterfaceIdByStep = zeros(1,nSteps,'uint32');
if ~isempty(temporalInterfaceNodes)
    temporalInterfaceIdByStep(temporalInterfaceNodes+1) = ...
        uint32(1:numel(temporalInterfaceNodes));
end
temporalInterfaceSideOffset = dt/4;

% --- 时间界面处的可选空间谱投影 ---
% 该选项只对周期网格定义；海绵边界破坏空间平移对称性，不能使用 FFT
% 掩码而仍保持明确的物理动量含义。
if ~isfield(cfg,'spectralFilterMask') || isempty(cfg.spectralFilterMask)
    spectralFilterMask = zeros(1,0);
else
    spectralFilterMask = cfg.spectralFilterMask(:).';
    if ~isnumeric(spectralFilterMask) && ~islogical(spectralFilterMask)
        error('cfg.spectralFilterMask must be a numeric or logical vector.');
    end
    if numel(spectralFilterMask) ~= Nx || ...
            any(~isfinite(double(spectralFilterMask))) || ...
            any(real(spectralFilterMask) < 0) || ...
            any(real(spectralFilterMask) > 1) || ...
            any(imag(spectralFilterMask) ~= 0)
        error(['cfg.spectralFilterMask must be a finite real vector of ' ...
            'length Nx with values in [0,1].']);
    end
    if ~strcmp(boundary,'periodic')
        error('cfg.spectralFilterMask requires boundary=''periodic''.');
    end
    if isempty(temporalInterfaces)
        error(['cfg.spectralFilterMask requires at least one explicit ' ...
            'cfg.temporalInterfaces value.']);
    end
    if useSingle
        spectralFilterMask = single(spectralFilterMask);
    else
        spectralFilterMask = double(spectralFilterMask);
    end
end

%==========================================================================
% 第 3 部分: 构造 Yee 网格（空间交错排布）
%==========================================================================

% H 网格位于相邻 E 格点的中点: xH(i) = x(i) + dx/2
switch boundary
    case 'periodic'
        % 周期边界: 最后一个 E 格点也参与 curl 计算（通过 circshift 环绕）
        % 因此 H 网格与 E 网格等长
        xH = x + dx/2;
    case 'sponge'
        % 吸收边界: 两端的 curl 不需要，H 网格比 E 网格少一个点
        % xH(i) 位于 x(i) 和 x(i+1) 的中点
        xH = x(1:end-1) + dx/2;
    otherwise
        error('cfg.boundary must be ''periodic'' or ''sponge''.');
end

%==========================================================================
% 第 4 部分: 初始条件设置
%==========================================================================

% --- E 场初始条件（t=0, E 网格）---
if isfield(cfg,'E0') && ~isempty(cfg.E0)
    E = cfg.E0(:).';            % 转为行向量
else
    E = zeros(1,Nx);            % 默认零场
end
if numel(E) ~= Nx
    error('cfg.E0 must have the same length as cfg.x.');
end

% --- H 场初始条件（t=-dt/2, H 网格）---
% 注意：因为 leapfrog 格式中 H 定义在半整数时间步，初始条件对应
% t = -dt/2 时刻。如果初始场是行波，需要通过波速和 dt/2 时移来构造。
if isfield(cfg,'Hhalf0') && ~isempty(cfg.Hhalf0)
    H = cfg.Hhalf0(:).';        % 转为行向量
else
    H = zeros(size(xH));        % 默认零场
end
if numel(H) ~= numel(xH)
    error('cfg.Hhalf0 has the wrong length for the selected boundary.');
end

%==========================================================================
% 第 5 部分: 通过本构关系初始化 D 和 B 场
%==========================================================================

% D(x, t=0) = ε(x, t=0) · E(x, t=0)
% 在 E 网格上计算 t=0 时刻的 ε_r
epsNow = cfg.epsFun(x,0);

% B(xH, t=-dt/2) = μ(xH, t=-dt/2) · H(xH, t=-dt/2)
% 在 H 网格上计算 t=-dt/2 时刻的 μ_r
muHalf = muFun(xH,-dt/2);

% 验证材料函数的输出尺寸正确
if numel(epsNow) ~= Nx || numel(muHalf) ~= numel(xH)
    error('epsFun or muFun returned an array with the wrong grid size.');
end

% 通过本构关系从 E/H 得到 D/B —— 之后只更新 D 和 B
D = epsNow.*E;                  % D = ε_r · E
B = muHalf.*H;                  % B = μ_r · H

% --- 若启用单精度，将工作数组转换为 single ---
if useSingle
    D = single(D);
    B = single(B);
    E = single(E);
    H = single(H);
end

%==========================================================================
% 第 6 部分: 采样 CFL 稳定性审计
%==========================================================================
%
% 对普通正折射率、无色散介质，稳定性条件为:
%   Courant 数 S = v_max · Δt/Δx < 1
%   其中 v_max = max[1/√(ε_r(x,t)·μ_r(x,t))]
%
% 此审计通过在仿真时间范围内采样若干时刻来估算 v_max。
% 如果用户已通过 cfg.maxWaveSpeed 指定了最大波速，则跳过采样。
% 对于色散/负折射率/有源介质（ε_r 或 μ_r 有虚部或实部为负），
% 审计返回 NaN，不进行 CFL 检查。

if isfield(cfg,'maxWaveSpeed') && ~isempty(cfg.maxWaveSpeed)
    % 用户手动指定了最大波速（例如已知材料参数或通过更严格的分析）
    sampledMaxWaveSpeed = cfg.maxWaveSpeed;
else
    % --- 确定采样时刻 ---
    if isfield(cfg,'stabilityTimes') && ~isempty(cfg.stabilityTimes)
        stabilityTimes = cfg.stabilityTimes(:).';
    else
        % 默认在 [0, nSteps·dt] 范围内均匀采样 9 个时刻
        stabilityTimes = linspace(0,nSteps*dt,9);
    end

    % --- 在所有采样时刻上计算整个空间网格的最大相速度 ---
    sampledMaxWaveSpeed = 0;
    for tq = stabilityTimes
        % 在当前采样时刻 tq，计算整个空间网格上的 ε_r 和 μ_r
        epsSample = cfg.epsFun(x,tq);
        muSample = muFun(x,tq);
        if numel(epsSample) ~= Nx || numel(muSample) ~= Nx
            error(['epsFun and muFun must return the same size as their ', ...
                'input grid during the CFL audit.']);
        end

        % 检查是否为普通正折射率介质:
        %   1. ε_r 和 μ_r 的虚部可忽略（无显著色散/损耗/增益）
        %   2. ε_r 和 μ_r 的实部为正（正折射率，无非正介质）
        ordinaryPositive = all(abs(imag(epsSample)) < 1e-12) ...
            && all(abs(imag(muSample)) < 1e-12) ...
            && all(real(epsSample) > 0) && all(real(muSample) > 0);

        if ~ordinaryPositive
            % 非正折射率或有色散/增益 —— 标量 CFL 不适用
            % 设为 NaN 以跳过后续的 Courant 数检查
            sampledMaxWaveSpeed = nan;
            break;
        end

        % 当前采样时刻的最大相速度 v = 1/√(ε_r·μ_r)
        % 对所有空间点取 max，再对所有采样时刻取 max
        sampledMaxWaveSpeed = max(sampledMaxWaveSpeed, ...
            max(1./sqrt(real(epsSample).*real(muSample))));
    end
end

% --- 计算 Courant 数并检查 ---
% Courant 数 S = v_max · Δt / Δx
% 对于真空（ε_r = μ_r = 1），退化为 S = c₀·Δt/Δx
sampledCourant = sampledMaxWaveSpeed*dt/dx;

if isfinite(sampledCourant) && sampledCourant >= 1
    % Courant 数 ≥ 1 意味着数值不稳定——场将指数增长发散
    % 用户必须减小 Δt 或增大 Δx（即增大 dx）
    % 对于色散/有源介质，可能需要提供自定义的本构更新方案
    error(['Sampled one-dimensional CFL number is %.6g >= 1. ', ...
        'Reduce cfg.dt or provide a validated constitutive update.'], ...
        sampledCourant);
end

%==========================================================================
% 第 7 部分: 构造 Sponge 吸收边界衰减剖面
%==========================================================================
%
% 海绵吸收层位于计算域的两端。每个时间步将靠近边界的场乘以衰减因子
% 来抑制反射波。衰减因子随距离呈三次方平滑变化：
%
%   damp(dist) = exp( -strength · s³ )
%   s = (spongeCells - dist) / spongeCells
%
% 其中 dist 是该点到吸收层内边界的距离（以网格点数为单位）。
% 三次方剖面确保衰减从零开始平滑过渡，使反射最小化。

dampE = ones(size(x));          % E 网格上的衰减因子，初始全 1（无衰减）
dampH = ones(size(xH));         % H 网格上的衰减因子，初始全 1（无衰减）

if strcmp(boundary,'sponge')
    % --- 吸收层参数 ---
    if ~isfield(cfg,'spongeCells') || isempty(cfg.spongeCells)
        % 默认吸收层宽度: 80 个格点或总格点数的 1/8，取较小值
        spongeCells = min(80,floor(Nx/8));
    else
        spongeCells = cfg.spongeCells;
    end
    if ~isfield(cfg,'spongeStrength') || isempty(cfg.spongeStrength)
        spongeStrength = 0.12;  % 默认吸收强度
    else
        spongeStrength = cfg.spongeStrength;
    end

    % --- 计算每个格点到最近边界的距离（以 dx 为单位）---
    % edgeDistanceE(i) = min( i-1, Nx-1-(i-1) )  即到左端或右端的较小者
    edgeDistanceE = min((x-x(1))/dx, (x(end)-x)/dx);
    edgeDistanceH = min((xH-x(1))/dx, (x(end)-xH)/dx);

    % --- 标记吸收层内的格点 ---
    % 只有距离小于 spongeCells 的格点才被衰减
    maskE = edgeDistanceE < spongeCells;
    maskH = edgeDistanceH < spongeCells;

    % --- 计算归一化穿透深度 s ∈ [0,1] ---
    % s=0 在吸收层内边界（无衰减），s=1 在计算域边界（最大衰减）
    sE = (spongeCells-edgeDistanceE(maskE))/spongeCells;
    sH = (spongeCells-edgeDistanceH(maskH))/spongeCells;

    % --- 三次方衰减剖面 ---
    % exp(-strength · s³): 比线性或二次方剖面反射更小
    dampE(maskE) = exp(-spongeStrength*sE.^3);
    dampH(maskH) = exp(-spongeStrength*sH.^3);
end

% --- 若启用单精度，将衰减因子也转为 single ---
if useSingle
    dampE = single(dampE);
    dampH = single(dampH);
end

%==========================================================================
% 第 8 部分: 分配历史记录数组
%==========================================================================

% 记录总次数 = floor(N_t / recordEvery) + 1（含第 0 步的初始状态）
nRecords = floor(nSteps/recordEvery) + 1;

% 根据精度选择基础类型，用于历史数组分配
if useSingle
    baseZero = single(0);
else
    baseZero = double(0);
end

if storeFields
    % 存储完整 E/H 历史: 每行是一个时刻，每列是一个空间格点
    EHistory = complex(zeros(nRecords,Nx,'like',baseZero));
    HHistory = complex(zeros(nRecords,Nx,'like',baseZero));
else
    % 低内存模式: 不分配全域历史空间
    EHistory = complex(zeros(0,Nx,'like',baseZero));
    HHistory = complex(zeros(0,Nx,'like',baseZero));
end

if storeD
    DHistory = complex(zeros(nRecords,Nx,'like',baseZero));
else
    DHistory = complex(zeros(0,Nx,'like',baseZero));
end

% 探针记录: 无论 storeFields 如何，探针总是被记录
probeE = complex(zeros(nRecords,numel(probeIndices),'like',baseZero));
probeH = complex(zeros(nRecords,numel(probeIndices),'like',baseZero));

% 记录时刻和能量（始终为 double，因为它们是时间标签和标量汇总）
tHistory = zeros(nRecords,1);
energy = zeros(nRecords,1);

%==========================================================================
% 第 9 部分: 记录初始状态（第 0 步）
%==========================================================================

recordId = 1;  % 当前记录索引

% 将 H/B 从 H 网格插值到 E 网格，以便在同一空间位置上记录 E 和 H
[HOnE, BOnE] = yee_h_to_e_grid(H,B,boundary,Nx);

if storeFields
    EHistory(recordId,:) = E;
    HHistory(recordId,:) = HOnE;  % 注意：存储的是插值到 E 网格的 H
end
if storeD
    DHistory(recordId,:) = D;
end

% 记录探针点的场
probeE(recordId,:) = E(probeIndices);
probeH(recordId,:) = HOnE(probeIndices);

% 记录初始时刻的瞬时电磁能量
% 能量 = (Δx/2) · Σ Re[ E*·D + H*·B ]
% 这是连续形式 ∫ (1/2)(E·D + H·B) dx 的离散近似
energy(recordId) = 0.5*dx*sum(real(conj(E).*D + conj(HOnE).*BOnE));

%==========================================================================
% 第 10 部分: Leapfrog 时间推进主循环
%==========================================================================
%
%  每个完整时间步（step）包含两个半时间步的子步骤：
%
%   ┌─────────────────────────────────────────────────────────┐
%   │ 步骤 A: 法拉第定律 —— 由 E^n 更新 B^{n+1/2}           │
%   │         ∂B/∂t = -∂E/∂x                                 │
%   │         ↓ 离散化                                        │
%   │         B(i)^{n+1/2} = B(i)^{n-1/2} - (Δt/Δx)·curlE(i) │
%   │         然后: H(i)^{n+1/2} = B(i)^{n+1/2} / μ(xH,tHalf)│
%   ├─────────────────────────────────────────────────────────┤
%   │ 步骤 B: 安培定律 —— 由 H^{n+1/2} 更新 D^{n+1}         │
%   │         ∂D/∂t = -∂H/∂x                                 │
%   │         ↓ 离散化                                        │
%   │         D(i)^{n+1} = D(i)^n - (Δt/Δx)·curlH(i)         │
%   │         然后: E(i)^{n+1} = D(i)^{n+1} / ε(x,tNow)      │
%   └─────────────────────────────────────────────────────────┘
%

% --- 初始化进度条 ---
if showProgress
    % 使用文本进度条（fprintf），兼容无图形界面的环境
    fprintf('FDTD 1D D/B-Yee: 共 %d 步, 开始推进...\n', nSteps);
    progressInterval = max(1, floor(nSteps/50));  % 每 2% 更新一次
    progressNextTick = progressInterval;
    ticProg = tic;
end

processedTemporalInterfaceCount = 0;
for step = 1:nSteps
    % =====================================================================
    % 步骤 A: 法拉第定律 —— 更新 B（从 n-1/2 到 n+1/2）
    % =====================================================================
    % 用 t = n·dt 时刻的 E 计算 B 的增量。若该整数节点是显式介电
    % 时间界面，则使用 D/ε(t^-) 与 D/ε(t^+) 的平均场。

    temporalInterfaceId = temporalInterfaceIdByStep(step);
    if temporalInterfaceId ~= 0
        processedTemporalInterfaceCount = ...
            processedTemporalInterfaceCount+1;
        tInterface = temporalInterfaces(double(temporalInterfaceId));
        epsBefore = cfg.epsFun(x,tInterface-temporalInterfaceSideOffset);
        epsAfter = cfg.epsFun(x,tInterface+temporalInterfaceSideOffset);
        if numel(epsBefore) ~= Nx || numel(epsAfter) ~= Nx
            error(['epsFun returned an array with the wrong grid size ', ...
                'around temporal interface %.16g.'],tInterface);
        end
        epsBefore = epsBefore(:).';
        epsAfter = epsAfter(:).';
        if any(~isfinite(real(epsBefore))) || ...
                any(~isfinite(imag(epsBefore))) || ...
                any(~isfinite(real(epsAfter))) || ...
                any(~isfinite(imag(epsAfter))) || ...
                any(abs(epsBefore) == 0) || any(abs(epsAfter) == 0)
            error(['epsFun must return finite nonzero permittivity on both ', ...
                'sides of temporal interface %.16g.'],tInterface);
        end
        if useSingle
            epsBefore = single(epsBefore);
            epsAfter = single(epsAfter);
        end
        EForBUpdate = 0.5*(D./epsBefore + D./epsAfter);
    else
        EForBUpdate = E;
    end

    if strcmp(boundary,'periodic')
        % 周期边界: curlE(i) = E(i+1) - E(i)，最后一点用 circshift 环绕
        % circshift(E,-1) 将 E 向左循环移位 1 位，即 E(2:end) 接 E(1)
        curlE = circshift(EForBUpdate,-1) - EForBUpdate;
    else
        % 吸收边界: curlE = diff(E)，即 ∂E/∂x ≈ [E(i+1)-E(i)]/dx
        % 仅在内部点计算，边界由 sponge 处理
        curlE = diff(EForBUpdate);
    end

    % B^{n+1/2} = B^{n-1/2} - (Δt/Δx)·curlE
    B = B - (dt/dx)*curlE;

    % 对 H 网格施加 sponge 衰减
    B = B.*dampH;

    % --- 通过本构关系恢复 H ---
    % 时间: tHalf = (step - 0.5)·dt，即 n+1/2 步
    tHalf = (step-0.5)*dt;
    muHalf = muFun(xH,tHalf);       % 在 H 网格、当前时刻采样 μ_r
    if useSingle; muHalf = single(muHalf); end
    H = B./muHalf;                  % H = B / (μ₀·μ_r)，本构关系

    % =====================================================================
    % 步骤 B: 安培定律 —— 更新 D（从 n 到 n+1）
    % =====================================================================
    % 用 t = (n+1/2)·dt 时刻的 H 计算 D 的增量

    if strcmp(boundary,'periodic')
        % 周期边界: curlH(i) = H(i) - H(i-1)
        % circshift(H,1) 将 H 向右循环移位 1 位，即 H(end) 接 H(1:end-1)
        curlH = H - circshift(H,1);
        D = D - (dt/dx)*curlH;
    else
        % 吸收边界: 只更新内部点（排除两端）
        % ∂H/∂x 在 E 格点 x(i) 处用中心差分 [H(i) - H(i-1)]/dx
        % 注意: D(1) 和 D(end) 不参与 curl 更新，仅靠 sponge 衰减
        D(2:end-1) = D(2:end-1) ...
            - (dt/dx)*(H(2:end)-H(1:end-1));
    end

    % --- 当前时刻 tNow = step·dt = n+1 ---
    tNow = step*dt;

    % --- 注入源项（可选）---
    % 通过向 D 添加增量实现"软源"注入，避免了硬源对反射波的阻挡
    if isfield(cfg,'sourceD') && ~isempty(cfg.sourceD)
        srcVal = cfg.sourceD(x,tNow,step);
        if useSingle; srcVal = single(srcVal); end
        D = D + srcVal;
    end

    % --- 施加 sponge 衰减 ---
    D = D.*dampE;

    % --- 通过本构关系恢复 E ---
    epsNow = cfg.epsFun(x,tNow);    % 在 E 网格、当前时刻采样 ε_r
    if useSingle; epsNow = single(epsNow); end
    E = D./epsNow;                  % E = D / (ε₀·ε_r)，本构关系

    % 有源时变介质会把任何落在高阶动量带隙中的舍入噪声指数放大。
    % 对已知为窄带解析信号的周期问题，可在真实材料时间界面处投影回
    % 用户指定的物理 k 支撑；D 和 B 必须一起投影以保持 Maxwell 状态。
    if temporalInterfaceId ~= 0 && ~isempty(spectralFilterMask)
        D = ifft(fft(D).*spectralFilterMask);
        B = ifft(fft(B).*spectralFilterMask);
        E = D./epsNow;
        H = B./muHalf;
    end

    % =====================================================================
    % 步骤 C: 按记录间隔保存场和能量
    % =====================================================================
    if mod(step,recordEvery) == 0
        recordId = recordId + 1;

        % 将 H/B 插值到 E 网格，以便统一记录
        [HOnE, BOnE] = yee_h_to_e_grid(H,B,boundary,Nx);

        if storeFields
            EHistory(recordId,:) = E;
            HHistory(recordId,:) = HOnE;
        end
        if storeD
            DHistory(recordId,:) = D;
        end

        % 探针记录
        probeE(recordId,:) = E(probeIndices);
        probeH(recordId,:) = HOnE(probeIndices);

        % 记录时刻
        tHistory(recordId) = tNow;

        % 计算瞬时电磁能量
        % 对非色散介质，能量密度 = (1/2)(E·D + H·B)
        % 注意：对时变介质，调制会与场交换能量，总能一般不守恒
        energy(recordId) = 0.5*dx*sum(real( ...
            conj(E).*D + conj(HOnE).*BOnE));
    end

    % --- 更新进度条 ---
    if showProgress && step >= progressNextTick
        elapsed = toc(ticProg);
        pct = step/nSteps*100;
        eta = elapsed/step*(nSteps-step);  % 预计剩余时间
        fprintf('  进度: %5.1f%% | 当前步: %d/%d | 已用: %.1fs | 预计剩余: %.1fs\n', ...
            pct, step, nSteps, elapsed, eta);
        progressNextTick = step + progressInterval;
    end
end

if processedTemporalInterfaceCount ~= numel(temporalInterfaces)
    error(['Processed %d temporal interfaces, but %d aligned interfaces ' ...
        'were configured.'],processedTemporalInterfaceCount, ...
        numel(temporalInterfaces));
end

% --- 关闭进度条 ---
if showProgress
    totalTime = toc(ticProg);
    fprintf('FDTD 1D D/B-Yee: 完成 %d 步, 总用时 %.1f s.\n', nSteps, totalTime);
end

%==========================================================================
% 第 11 部分: 组装输出结构体
%==========================================================================

out.x = x;                          % E 空间网格
out.xH = xH;                        % H 空间网格（偏移 dx/2）
out.t = tHistory;                   % 记录时刻数组
out.E = EHistory;                   % 电场历史（每行一个时刻）
out.H = HHistory;                   % 磁场历史（已插值到 E 网格，每行一个时刻）
out.D = DHistory;                   % 电位移历史（每行一个时刻）
out.probeIndices = probeIndices;    % 探针索引
out.probeX = x(probeIndices);       % 探针空间位置
out.probeE = probeE;                % 探针电场记录
out.probeH = probeH;                % 探针磁场记录
out.energy = energy;                % 能量时间序列
out.finalD = D;                     % 最终时刻 D 场
out.finalB = B;                     % 最终时刻 B 场
out.finalE = E;                     % 最终时刻 E 场
out.finalH = H;                     % 最终时刻 H 场
out.dx = dx;                        % 空间步长
out.dt = dt;                        % 时间步长
out.boundary = boundary;            % 边界条件类型
out.storeFields = storeFields;      % 场存储模式
out.storeD = storeD;                % D 历史存储模式
out.temporalInterfaces = temporalInterfaces;
out.processedTemporalInterfaceCount = processedTemporalInterfaceCount;
out.temporalInterfaceTolerance = temporalInterfaceTolerance;
out.spectralFilterMask = spectralFilterMask;
out.sampledMaxWaveSpeed = sampledMaxWaveSpeed;  % CFL 审计波速
out.sampledCourant = sampledCourant;            % CFL 审计 Courant 数
out.precision = precision;                      % 数值精度

end

%==========================================================================
% 局部辅助函数: 将 H/B 从 H 网格（半整数）插值到 E 网格（整数）
%==========================================================================
%
% 由于 Yee 网格中 E 和 H 在空间上错开 dx/2，我们通常希望在同一空间
% 位置上（E 网格）输出两个场以便绘图和分析。此函数做线性插值：
%
%   E 网格的 H(i) ≈ [H(i-1) + H(i)]/2    （内部点）
%   E 网格的 H(1) ≈ H(1)                 （左边界，零阶外推）
%   E 网格的 H(end) ≈ H(end)             （右边界，零阶外推）
%
% 对于周期边界，所有点都可以用环绕的平均值。
%
function [HOnE, BOnE] = yee_h_to_e_grid(H,B,boundary,Nx)
if strcmp(boundary,'periodic')
    % 周期边界: 所有点都可以通过相邻 H 值的算术平均来插值
    % circshift(H,1) 得到 H(i-1)，与 H(i) 平均得到 E 网格 x(i) 处的值
    HOnE = 0.5*(H + circshift(H,1));
    BOnE = 0.5*(B + circshift(B,1));
else
    % 吸收边界: 内部点用算术平均，边界点用单侧值
    % 预分配 E 网格大小的数组
    HOnE = zeros(1,Nx,'like',H);
    BOnE = zeros(1,Nx,'like',B);

    % 内部点 (i = 2 ... Nx-1):
    %   HOnE(i) ≈ [H(i-1) + H(i)] / 2
    %   H(i-1) 在 x(i)-dx/2 处，H(i) 在 x(i)+dx/2 处
    %   平均给出 x(i) 处的二阶插值
    HOnE(2:end-1) = 0.5*(H(1:end-1)+H(2:end));
    BOnE(2:end-1) = 0.5*(B(1:end-1)+B(2:end));

    % 端点: 无法做中心插值，直接使用最近的 H 值
    % 这只影响 sponge 吸收层内的（不重要）边界点
    HOnE(1)   = H(1);
    HOnE(end) = H(end);
    BOnE(1)   = B(1);
    BOnE(end) = B(end);
end
end
