%RUN_TMM Entry script: Floquet band structure of a binary photonic time
%  crystal from the exact 2x2 D/B monodromy (transfer-matrix method).
%
% A plane wave exp(i*k*x) propagates through layers that are uniform in
% space but switch in time. Within each layer the D/B state evolves by
%     d/dt [D;B] = -i*k*[0 1/mu; 1/eps 0]*[D;B],
% and D,B are continuous at every time interface, so the one-period
% evolution U is a product of layer exponentials. Its eigenvalues lambda
% give the Floquet quasifrequencies omega = i*log(lambda)/T (principal
% log, Re(omega) in the first temporal Brillouin zone). A momentum gap is
% where |trace(U)/2| > 1, i.e. omega picks up an imaginary part.
%
% Plan: 定义材料 -> 周期/截断参数 -> 均匀介质解析检验 -> 计算 -> 作图
%
% A uniform-medium analytic check (omega = k/sqrt(eps*mu), no gap) is
% asserted before the main computation (acceptance repair.md 7.308).
%
% Results remain in the workspace; by default nothing is written to disk.
% Set doSave = true (end of file) to export the figure to a local dir.

clear; clc; close all;

%% ======================= 定义材料 (material) =============================
% Binary photonic time crystal: two lossless layers switching in time.
epsLayers = [1 4];          % permittivity of each temporal layer
muLayers  = 1;              % permeability (scalar -> broadcast to all layers)
durations = [0.5 0.5];      % duration of each layer (arbitrary units)

% A uniform-medium analytic check (omega = k/sqrt(eps*mu), no momentum gap)
% runs below as a separate self-contained block; it does not alter the demo
% material configured above.

%% ================= 周期 / 截断参数 (period & scan) ========================
T     = sum(durations);     % total modulation period
Omega = 2*pi/T;             % temporal modulation frequency

kNorm   = linspace(0, 1.5, 301);   % normalized wavenumber k/Omega (c0 = 1)
kValues = kNorm*(2*pi/T);          % physical wavenumbers k = kNorm*Omega

%% ================== 均匀介质解析检验 (uniform-medium check) =================
% Acceptance repair.md 7.308: for eps = const, mu = 1, the TMM monodromy
% must recover the folded dispersion omega = k/sqrt(eps*mu) to machine
% precision (no momentum gap). Two identical half-period layers of a
% uniform medium are exactly equivalent to one uniform period.
epsU   = 2.25;  nU = sqrt(epsU);
kNormU = linspace(0, 1.2, 101);                 % normalized scan, same as PWE check
kU     = kNormU*Omega;
wExact = mod(kU/nU + Omega/2, Omega) - Omega/2; % analytic Re(omega), folded into [-Omega/2, Omega/2]
bandsU = tmm_bands(kU, [epsU epsU], 1, [T/2 T/2]);
dU     = mod(real(bandsU.omega) - wExact + Omega/2, Omega) - Omega/2;   % wrap onto period-Omega circle
minErrU = min(abs(dU), [], 1)/Omega;            % nearest branch, like the PWE-vs-TMM comparison
maxErrU = max(minErrU);
uniformTol = 1e-6;
assert(maxErrU <= uniformTol, ...
    ['TMM uniform-medium recovery max |Delta Re(omega)|/Omega = %.3e exceeds %.1e ', ...
    '(acceptance 7.308); monodromy does not return omega = k/sqrt(eps*mu).'], ...
    maxErrU, uniformTol);
fprintf('TMM uniform-medium analytic check: max |Delta Re(omega)|/Omega = %.3e (<= %.1e) PASSED\n', ...
    maxErrU, uniformTol);

%% ============================ 计算 (compute) =============================
result = tmm_bands(kValues, epsLayers, muLayers, durations);

% Convenience handles (band structure is 2 x nK, one row per branch).
omega  = result.omega;               % 2 x nK complex quasifrequencies
halfTr = result.halfTrace;           % 1 x nK complex trace(U)/2
gap    = result.gapMask;             % 1 x nK logical momentum-gap mask

% Locate contiguous momentum-gap segments (runs of true in gapMask).
d        = diff([false, gap(:).', false]);
gapStart = find(d ==  1);            % first index of each run
gapStop  = find(d == -1) - 1;        % last  index of each run
nGaps    = numel(gapStart);

%% ========================== 作图 (plot) ==================================
% White figure, three tiles: Re(omega), Im(omega), and the half-trace
% criterion |trace(U)/2| vs 1.
fig = figure('Color','w','Position',[80 80 1320 400]);
tl  = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');

c1 = [0.00 0.45 0.75];   % branch 1 (blue)
c2 = [0.85 0.33 0.10];   % branch 2 (orange)
gapShade = [0.97 0.90 0.90];  % light tint for momentum-gap regions

% ---- panel 1: real part of the Floquet quasifrequency -------------------
ax1 = nexttile(tl);
hold(ax1,'on');
yl1 = [-0.55 0.55];                      % folded Re(omega)/Omega band range
for s = 1:nGaps                          % shade each momentum-gap segment
    xg = kNorm([gapStart(s) gapStop(s) gapStop(s) gapStart(s)]);
    patch(ax1, xg, [yl1(1) yl1(1) yl1(2) yl1(2)], gapShade, ...
        'EdgeColor','none', 'HandleVisibility','off');
end
plot(ax1, kNorm, real(omega(1,:))/Omega, '-', 'Color',c1, 'LineWidth',1.6);
plot(ax1, kNorm, real(omega(2,:))/Omega, '-', 'Color',c2, 'LineWidth',1.6);
xlabel(ax1, 'k/\Omega');
ylabel(ax1, 'Re(\omega)/\Omega');
title(ax1, 'Floquet bands (momentum gaps shaded)');
xlim(ax1,[kNorm(1) kNorm(end)]); ylim(ax1, yl1);
grid(ax1,'on'); box(ax1,'on');

% ---- panel 2: imaginary part -> growth/decay inside the gap -------------
ax2 = nexttile(tl);
hold(ax2,'on');
plot(ax2, kNorm, imag(omega(1,:))/Omega, '-', 'Color',c1, 'LineWidth',1.6);
plot(ax2, kNorm, imag(omega(2,:))/Omega, '-', 'Color',c2, 'LineWidth',1.6);
yline(ax2, 0, 'k:');
xlabel(ax2, 'k/\Omega');
ylabel(ax2, 'Im(\omega)/\Omega');
title(ax2, 'Growth / decay within momentum gap');
xlim(ax2,[kNorm(1) kNorm(end)]);
grid(ax2,'on'); box(ax2,'on');

% ---- panel 3: gap criterion |trace(U)/2| vs 1 ---------------------------
ax3 = nexttile(tl);
hold(ax3,'on');
plot(ax3, kNorm, real(halfTr), '-', 'Color',[0.30 0.30 0.30], ...
    'LineWidth',1.4);
yline(ax3,  1, 'k--', 'HandleVisibility','off');
yline(ax3, -1, 'k--', 'HandleVisibility','off');
xlabel(ax3, 'k/\Omega');
ylabel(ax3, 'Re(tr U / 2)');
title(ax3, 'Momentum-gap criterion');
xlim(ax3,[kNorm(1) kNorm(end)]);
yl3 = [min([real(halfTr), -1]) - 0.2, max([real(halfTr), 1]) + 0.2];
ylim(ax3, yl3);
grid(ax3,'on'); box(ax3,'on');

title(tl, sprintf('Binary photonic time crystal  \\epsilon = [%s],  T = %.2f', ...
    num2str(epsLayers,'%.3g '), T));

drawnow;

%% ========================= 控制台输出 (console) ===========================
fprintf('============================================================\n');
fprintf('T     = %.6f\n', T);
fprintf('Omega = %.6f  (2*pi/T)\n', Omega);
fprintf('epsLayers = [%s]\n', num2str(epsLayers,'%.6g '));
fprintf('durations = [%s]\n', num2str(durations,'%.6g '));
fprintf('Momentum gaps (|Re(tr U / 2)| > 1): %d found\n', nGaps);
if nGaps > 0
    for s = 1:nGaps
        fprintf('  gap %d: k/Omega in [%.4f, %.4f]   (k in [%.4f, %.4f])\n', ...
            s, kNorm(gapStart(s)), kNorm(gapStop(s)), ...
               kValues(gapStart(s)), kValues(gapStop(s)));
    end
    maxImOverOmega = max(abs(imag(omega(:, gap))), [], 'all') / Omega;
    fprintf('max |Im(omega)|/Omega inside gaps = %.6e\n', maxImOverOmega);
else
    maxImOverOmega = 0;
    fprintf('  (no momentum gap in the scanned k-range)\n');
end
fprintf('Base MATLAB only (no Toolbox); declared minimum R2020a (untested on this machine), verified R2026a.\n');
fprintf('============================================================\n');

%% ====================== 可选保存 (optional save) =========================
% By default nothing is written and no directory is created. Flip this
% guard to true to export the figure to a local 'output' dir.
doSave = false;
if doSave
    outDir = fullfile(pwd, 'output');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    outFile = fullfile(outDir, 'run_tmm_bands.png');
    try
        exportgraphics(fig, outFile, 'Resolution', 220);
    catch
        print(fig, outFile, '-dpng', '-r220');
    end
    fprintf('Saved figure to %s\n', outFile);
end
