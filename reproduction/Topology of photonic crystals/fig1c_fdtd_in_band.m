function fig1c_fdtd_in_band()
%FIG1C_FDTD_IN_BAND  Reproduce Fig. 1(c) of Lustig et al. (arXiv:1803.08731).
%
%   FDTD simulation of a wave packet whose central momentum falls within
%   a Floquet band of the binary PTC. Upon entering the finite-duration
%   PTC, the pulse splits into forward- and backward-propagating Floquet
%   modes. Each mode splits again upon exit, yielding four emergent pulses.
%
%   Reference:
%     E. Lustig, Y. Sharabi, and M. Segev,
%     "Topology of photonic time-crystals," arXiv:1803.08731v1 (2018).

% --- Add parent toolbox to path ---
rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(rootDir, 'startup_stm.m'));

% =========================================================================
% PTC parameters
% =========================================================================
eps1 = 3;   eps2 = 1;
mu1  = 1;   mu2  = 1;
T    = 2*pi;
t1   = 0.5*T;
t2   = 0.5*T;
epsBg = (eps1 + eps2)/2;
muBg  = 1;
nBg   = sqrt(epsBg);

% =========================================================================
% Find an in-band k value
% =========================================================================
fprintf('=== Identifying in-band momentum ===\n');
kScan = linspace(0.01, 2.5, 2500);
bands = temporal_crystal_bands(kScan, [eps1 eps2], [mu1 mu2], [t1 t2]);
halfTr = bands.halfTrace;
inBand = abs(halfTr) <= 1;

% Find band segments (no toolbox needed)
d = diff([false, inBand, false]);
bStarts = find(d == 1);
bEnds   = find(d == -1) - 1;
nBands = length(bStarts);
fprintf('Found %d bands\n', nBands);

% Use the 3rd band (or the one nearest to a target k ~ 1.0)
targetK = 1.0;
bestBand = 1;
bestDist = inf;
for b = 1:nBands
    kMid = (kScan(bStarts(b)) + kScan(bEnds(b))) / 2;
    dist = abs(kMid - targetK);
    if dist < bestDist && (bEnds(b) - bStarts(b)) > 20
        bestDist = dist;
        bestBand = b;
    end
end
bandIdx = bestBand;

% Pick k at band center
kMidIdx = round((bStarts(bandIdx) + bEnds(bandIdx)) / 2);
k0 = kScan(kMidIdx);
lambda0 = 2*pi / k0;
bandKRange = [kScan(bStarts(bandIdx)), kScan(bEnds(bandIdx))];
fprintf('Band %d: k ∈ [%.4f, %.4f], chosen k0 = %.4f\n', ...
    bandIdx, bandKRange(1), bandKRange(2), k0);
fprintf('  |halfTr(k0)| = %.4f, wavelength = %.4f\n', ...
    abs(halfTr(kMidIdx)), lambda0);

% =========================================================================
% FDTD grid
% =========================================================================
dx = lambda0 / 28;
xDomain = 140;        % spatial extent
x = 0:dx:xDomain;
Nx = numel(x);
fprintf('FDTD grid: Nx = %d, dx = %.4f\n', Nx, dx);

dt = 0.8 * dx;
fprintf('dt = %.4f, Courant = %.4f\n', dt, dt/dx);

% =========================================================================
% Initial pulse
% =========================================================================
x0 = xDomain * 0.22;
sigma = 2.5;  % spatial half-width (FWHM ≈ 5.9)
pulseFn = @(xq) exp(-((xq-x0)/sigma).^2) .* exp(1i*k0*(xq-x0));

E0 = pulseFn(x);
vBg = 1/nBg;
xH = x(1:end-1) + dx/2;
Hhalf0 = nBg * pulseFn(xH + vBg*dt/2);

% =========================================================================
% PTC temporal window
% =========================================================================
tStart = 220;
tEnd   = 340;
nPeriods = round((tEnd - tStart) / T);
tEnd = tStart + nPeriods * T;
fprintf('PTC: t ∈ [%.1f, %.1f], N = %d periods\n', tStart, tEnd, nPeriods);

tSim = tEnd + 120;
nSteps = ceil(tSim / dt);
fprintf('Simulation: %d steps, t_max = %.1f\n', nSteps, nSteps*dt);

% =========================================================================
% Run FDTD
% =========================================================================
cfg.x = x;
cfg.dt = dt;
cfg.nSteps = nSteps;
cfg.epsFun = @(xq, tq) ptcEps(xq, tq, epsBg, eps1, eps2, t1, t2, tStart, tEnd);
cfg.muFun = @(xq, tq) ones(size(xq));
cfg.E0 = E0;
cfg.Hhalf0 = Hhalf0;
cfg.boundary = 'sponge';
cfg.spongeCells = 140;
cfg.spongeStrength = 0.10;
cfg.recordEvery = 5;
cfg.storeFields = true;
cfg.progressBar = true;

fprintf('\n=== FDTD: in-band pulse ===\n');
out = fdtd1d_db(cfg);

% =========================================================================
% Create figure
% =========================================================================
fig = figure('Color', 'w', 'Position', [60 60 1150 780]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact');

% --- Main panel: |E(x,t)| ---
ax1 = nexttile(tl, 1, [1 2]);
Eabs = abs(out.E);
Eabs(Eabs < 1e-10) = 1e-10;
imagesc(ax1, out.x, out.t, log10(Eabs));
set(ax1, 'YDir', 'normal');
hold(ax1, 'on');

% PTC window highlighting
yline(ax1, tStart, 'w--', 'LineWidth', 1.5);
yline(ax1, tEnd, 'w--', 'LineWidth', 1.5);
xL = ax1.XLim;
fill(ax1, [xL(1) xL(2) xL(2) xL(1)], [tStart tStart tEnd tEnd], ...
    'white', 'FaceAlpha', 0.10, 'EdgeColor', 'none');

cbar = colorbar(ax1);
cbar.Label.String = 'log_{10}|E|';
cbar.Label.FontSize = 11;
xlabel(ax1, 'Space  x', 'FontSize', 12);
ylabel(ax1, 'Time  t', 'FontSize', 12);
title(ax1, 'Fig. 1(c): |E(x,t)| — Pulse in a Floquet band', ...
    'FontSize', 14, 'FontWeight', 'bold');
colormap(ax1, jet(256));
set(ax1, 'FontSize', 11);

% --- Probe time trace ---
ax2 = nexttile(tl);
[~, probeIdx] = min(abs(out.x - x0));
semilogy(ax2, out.t, abs(out.E(:, probeIdx)), 'b-', 'LineWidth', 1.2);
hold(ax2, 'on');
xline(ax2, tStart, 'k--', 'LineWidth', 1);
xline(ax2, tEnd, 'k--', 'LineWidth', 1);
xlabel(ax2, 'Time  t', 'FontSize', 11);
ylabel(ax2, '|E| at pulse center', 'FontSize', 11);
title(ax2, 'On-axis field amplitude', 'FontSize', 12);
grid(ax2, 'on'); box(ax2, 'on'); set(ax2, 'FontSize', 10);

% --- Spatial snapshots ---
ax3 = nexttile(tl);
[~, beforeId] = min(abs(out.t - (tStart - 8)));
[~, afterId]  = min(abs(out.t - (tEnd + 30)));
plot(ax3, out.x, abs(out.E(beforeId, :)), 'k-', 'LineWidth', 1.1);
hold(ax3, 'on');
plot(ax3, out.x, abs(out.E(afterId, :)), 'r-', 'LineWidth', 1.4);
xlabel(ax3, 'Space  x', 'FontSize', 11);
ylabel(ax3, '|E|', 'FontSize', 11);
legend(ax3, sprintf('Before (t=%.0f)', out.t(beforeId)), ...
    sprintf('After (t=%.0f)', out.t(afterId)), ...
    'Location', 'best', 'FontSize', 9);
title(ax3, 'Spatial profiles: before vs after PTC', 'FontSize', 11);
grid(ax3, 'on'); box(ax3, 'on'); set(ax3, 'FontSize', 10);

% --- Energy ---
ax4 = nexttile(tl);
semilogy(ax4, out.t, out.energy/out.energy(1), 'b-', 'LineWidth', 1.3);
hold(ax4, 'on');
xline(ax4, tStart, 'k--', 'LineWidth', 1);
xline(ax4, tEnd, 'k--', 'LineWidth', 1);
xlabel(ax4, 'Time  t', 'FontSize', 11);
ylabel(ax4, 'E / E_0', 'FontSize', 11);
title(ax4, 'Normalized energy (band → stable)', 'FontSize', 12);
grid(ax4, 'on'); box(ax4, 'on'); set(ax4, 'FontSize', 10);

% --- Overall title ---
title(tl, sprintf(['In-band: k_0=%.4f, \\lambda_0=%.3f, ' ...
    'band %d, |Tr/2|=%.3f'], k0, lambda0, bandIdx, abs(halfTr(kMidIdx))), ...
    'FontSize', 10, 'FontWeight', 'normal');

% =========================================================================
% Save
% =========================================================================
outputDir = fullfile(fileparts(mfilename('fullpath')), 'output');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

outputFile = fullfile(outputDir, 'fig1c_fdtd_in_band.png');
try
    exportgraphics(fig, outputFile, 'Resolution', 150);
catch
    print(fig, outputFile, '-dpng', '-r150');
end
fprintf('\nSaved: %s\n', outputFile);

save(fullfile(outputDir, 'fig1c_data.mat'), ...
    'out', 'k0', 'lambda0', 'bandIdx', 'tStart', 'tEnd', ...
    'eps1', 'eps2', 'T', 't1', 't2', 'epsBg');
fprintf('Saved data to fig1c_data.mat\n');
end

% =========================================================================
% Permittivity function
% =========================================================================
function epsVal = ptcEps(x, t, epsBg, eps1, eps2, t1, t2, tStart, tEnd)
if t < tStart || t > tEnd
    epsVal = epsBg * ones(size(x));
else
    tLocal = t - tStart;
    Tper = t1 + t2;
    tPhase = mod(tLocal, Tper);
    if tPhase < t1
        epsVal = eps1 * ones(size(x));
    else
        epsVal = eps2 * ones(size(x));
    end
end
end
