%% test_fdtd_optimizations.m
% 测试 B1 (single precision), B4 (progress bar), A1 (animation)
% 使用归一化单位 (c0 = 1)，与项目 demo 保持一致

fprintf('=== 测试 1: 双精度 + 进度条 ===\n');

% --- 归一化单位仿真参数 ---
% 参照 demo05/demo06 的归一化约定：c0 = 1, lambda0 = 1
c0 = 1;
lambda0 = 1;
Lf = 32;           % 每波长采样点数
Lx = 30;           % 波长数
k0 = 2*pi/lambda0;

dx = lambda0 / Lf;
x = (0:dx:Lx*lambda0).';   % 列向量
dt = 0.80 * dx / c0;       % CFL = 0.80
tEnd = 15 * lambda0 / c0;  % 15 个周期
nSteps = ceil(tEnd / dt);

% --- 高斯波包初始条件 ---
% 参照 demo05 的约定：正向传播波 H = n * E (阻抗匹配)
x0 = 5 * lambda0;
sigma = 1.5 * lambda0;
profile = @(xq) exp(-((xq-x0)/sigma).^2) .* exp(1i*k0*(xq-x0));
E0 = profile(x).';

% H 初始条件 (t = -dt/2, H 网格), 正向波: H(x, -dt/2) = E(x + c*dt/2, 0)
xH = x(1:end-1) + dx/2;
Hhalf0 = profile(xH.' + c0*dt/2).';  % 真空 n=1, 所以 H = 1 * E

% --- 配置 ---
cfg.x = x;
cfg.dt = dt;
cfg.nSteps = nSteps;
cfg.epsFun = @(xq, tq) ones(size(xq));
cfg.muFun = @(xq, tq) ones(size(xq));
cfg.E0 = E0;
cfg.Hhalf0 = Hhalf0;
cfg.boundary = 'sponge';
cfg.spongeCells = 80;
cfg.spongeStrength = 0.10;
cfg.recordEvery = 2;
cfg.progressBar = true;      % B4: 测试进度条
cfg.precision = 'double';    % B1: 双精度基准
cfg.storeFields = true;

out_double = fdtd1d_db(cfg);

% --- 诊断: 检查波包是否确实传播 ---
nRec = size(out_double.E, 1);
[~, peak0] = max(abs(out_double.E(1, :)));
[~, peakMid] = max(abs(out_double.E(round(nRec/2), :)));
[~, peakEnd] = max(abs(out_double.E(end, :)));
distTraveled = out_double.x(peakEnd) - out_double.x(peak0);
distExpected = c0 * out_double.t(end);
fprintf('\n诊断 —— 波包传播:\n');
fprintf('  峰值位置:  t=0: idx=%d (x=%.2f),  t=T/2: idx=%d (x=%.2f),  t=T: idx=%d (x=%.2f)\n', ...
    peak0, out_double.x(peak0), peakMid, out_double.x(peakMid), peakEnd, out_double.x(peakEnd));
fprintf('  实际传播距离: %.2f  (期望: %.2f)\n', distTraveled, distExpected);
fprintf('  误差: %.2f%%\n', abs(distTraveled - distExpected)/distExpected * 100);

if distTraveled < 0.01 * distExpected
    fprintf('  ✗ 错误: 波包几乎不动，仿真未正确传播！\n');
else
    fprintf('  ✓ 波包正常传播\n');
end

fprintf('\n=== 测试 2: 单精度 + 进度条 ===\n');

cfg.precision = 'single';    % B1: 单精度
out_single = fdtd1d_db(cfg);

% --- 比较双精度与单精度的差异 ---
diffE = abs(out_double.E(:) - double(out_single.E(:)));
maxDiff = max(diffE);
relDiff = maxDiff / max(abs(out_double.E(:)));
fprintf('\n双精度 vs 单精度:\n');
fprintf('  最大绝对误差: %.3e\n', maxDiff);
fprintf('  最大相对误差: %.3e\n', relDiff);
if relDiff < 1e-3
    fprintf('  ✓ 单精度误差在可接受范围内 (< 0.1%%)\n');
else
    fprintf('  ⚠ 单精度误差较大，建议检查\n');
end

% --- 内存对比 ---
mem_double = whos('out_double');
mem_single = whos('out_single');
fprintf('\n内存占用:\n');
fprintf('  double: %.2f MB\n', mem_double.bytes / 1e6);
fprintf('  single: %.2f MB\n', mem_single.bytes / 1e6);
if mem_double.bytes > 0
    fprintf('  节省:   %.1f%%\n', (1 - mem_single.bytes/mem_double.bytes)*100);
end

fprintf('\n=== 测试 3: 动画生成 ===\n');

% Line 模式动画 → GIF
fprintf('  生成 line 模式动画...\n');
animate_fdtd1d(out_double, struct(...
    'saveAs', 'test_wave_line.gif', ...
    'viewMode', 'line', ...
    'fields', 'E', ...
    'fps', 15, ...
    'skipFrames', 4, ...
    'showTime', true, ...
    'speed', 10));

% Waterfall 模式动画 → GIF
fprintf('  生成 waterfall 模式动画...\n');
animate_fdtd1d(out_double, struct(...
    'saveAs', 'test_wave_waterfall.gif', ...
    'viewMode', 'waterfall', ...
    'fps', 15, ...
    'skipFrames', 4, ...
    'speed', 10));

% --- 诊断图: 绘制 E(x) 在三个时刻的快照 ---
fprintf('\n=== 生成诊断快照 ===\n');
figDiag = figure('Color', 'w', 'Position', [100 100 900 400], 'Visible', 'off');
hold on;
plot(out_double.x, real(out_double.E(1, :)), 'b-', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('t = %.2f', out_double.t(1)));
plot(out_double.x, real(out_double.E(round(nRec/2), :)), 'r-', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('t = %.2f', out_double.t(round(nRec/2))));
plot(out_double.x, real(out_double.E(end, :)), 'k-', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('t = %.2f', out_double.t(end)));
xlabel('x / \lambda_0');
ylabel('Re(E_y)');
title('波包传播诊断: t=0, t=T/2, t=T');
legend('Location', 'best');
grid on;
saveas(figDiag, 'test_diagnostic_snapshot.png');
close(figDiag);
fprintf('  快照已保存至: test_diagnostic_snapshot.png\n');

fprintf('\n=== 全部测试完成 ===\n');
fprintf('生成文件:\n');
fprintf('  test_wave_line.gif       — line 模式动画\n');
fprintf('  test_wave_waterfall.gif  — waterfall 模式动画\n');
fprintf('  test_diagnostic_snapshot.png — 诊断快照\n');
