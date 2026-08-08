function fig2_fdtd_simulations()
%FIG2_FDTD_SIMULATIONS  Reproduce Fig. 2 of Lustig et al., Optica 5, 1390 (2018).
%
%   Two-panel figure from FDTD simulations:
%     (a) Pulse in a momentum band — enters PTC, splits into two Floquet
%         modes, exits, splits again → four output pulses.
%     (b) Pulse in a momentum bandgap — exponential growth during PTC,
%         then two output pulses.
%
%   Both pulses have FWHM ~ 45 fs. The in-band pulse has a center
%   wavelength of 1.4 um; the in-gap pulse has center wavelength 0.93 um.
%   The PTC starts at t=220 fs and ends at t=340 fs (60 periods).
%
%   NOTE: FDTD simulations are time-consuming (~3-10 min each).
%   Set doFDTD = false to load previously saved data.
%
%   Reference:
%     E. Lustig, Y. Sharabi, and M. Segev,
%     "Topological aspects of photonic time crystals,"
%     Optica 5, 1390-1395 (2018).  DOI: 10.1364/OPTICA.5.001390

% --- Add parent toolbox to path ---
rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(rootDir, 'startup_stm.m'));

% =========================================================================
% Configuration
% =========================================================================
doFDTD = true;   % set to false to load saved data
outputDir = fullfile(fileparts(mfilename('fullpath')), 'output');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

% =========================================================================
% PTC parameters (matching the paper)
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
% Part A: In-band pulse
% =========================================================================
dataFileInBand = fullfile(outputDir, 'fig2_inband_data.mat');

if doFDTD || ~exist(dataFileInBand, 'file')
    fprintf('=== Part A: In-band FDTD simulation ===\n');

    % Find in-band k value (~1.4 um wavelength → k ≈ 2π/1.4 ≈ 4.5...
    % but with dimensionless units, we find k from band structure)
    kScan = linspace(0.01, 2.5, 2500);
    bands = temporal_crystal_bands(kScan, [eps1 eps2], [mu1 mu2], [t1 t2]);
    halfTr = bands.halfTrace;
    inBand = abs(halfTr) <= 1;

    d = diff([false, inBand, false]);
    bStarts = find(d == 1);
    bEnds   = find(d == -1) - 1;
    nBands = length(bStarts);

    % Target the 3rd band for good visibility
    targetK = 1.0;
    bestBand = 1;  bestDist = inf;
    for b = 1:nBands
        kMid = (kScan(bStarts(b)) + kScan(bEnds(b))) / 2;
        dist = abs(kMid - targetK);
        if dist < bestDist && (bEnds(b) - bStarts(b)) > 20
            bestDist = dist;  bestBand = b;
        end
    end
    bandIdx = bestBand;
    kMidIdx = round((bStarts(bandIdx) + bEnds(bandIdx)) / 2);
    k0_inband = kScan(kMidIdx);
    lambda0_inband = 2*pi / k0_inband;
    fprintf('  Band %d: k0=%.4f, lambda=%.4f, |Tr/2|=%.3f\n', ...
        bandIdx, k0_inband, lambda0_inband, abs(halfTr(kMidIdx)));

    % FDTD grid
    dx = lambda0_inband / 28;
    xDomain = 140;
    x = 0:dx:xDomain;
    Nx = numel(x);
    dt = 0.8 * dx;
    fprintf('  Grid: Nx=%d, dx=%.4f, dt=%.4f\n', Nx, dx, dt);

    % Initial pulse
    x0 = xDomain * 0.22;
    sigma = 2.5;
    pulseFn = @(xq) exp(-((xq-x0)/sigma).^2) .* exp(1i*k0_inband*(xq-x0));
    E0_inband = pulseFn(x);
    vBg = 1/nBg;
    xH = x(1:end-1) + dx/2;
    H0_inband = nBg * pulseFn(xH + vBg*dt/2);

    % PTC window
    tStart = 220;  tEnd = 340;
    nPeriods = round((tEnd - tStart) / T);
    tEnd = tStart + nPeriods * T;
    tSim = tEnd + 120;
    nSteps = ceil(tSim / dt);
    fprintf('  PTC: [%.1f, %.1f], %d periods, %d steps\n', ...
        tStart, tEnd, nPeriods, nSteps);

    cfg = struct('x', x, 'dt', dt, 'nSteps', nSteps, ...
        'epsFun', @(xq,tq) ptcEps(xq,tq,epsBg,eps1,eps2,t1,t2,tStart,tEnd), ...
        'muFun', @(xq,tq) ones(size(xq)), ...
        'E0', E0_inband, 'Hhalf0', H0_inband, ...
        'boundary', 'sponge', 'spongeCells', 140, 'spongeStrength', 0.10, ...
        'recordEvery', 5, 'storeFields', true, 'progressBar', true);
    out_inband = fdtd1d_db(cfg);
    halfTr_inband = abs(halfTr(kMidIdx));
    save(dataFileInBand, 'out_inband', 'k0_inband', 'lambda0_inband', ...
        'bandIdx', 'tStart', 'tEnd', 'nPeriods', 'x0', 'halfTr_inband');
    fprintf('  Saved in-band data.\n');
else
    fprintf('=== Loading saved in-band data ===\n');
    load(dataFileInBand, 'out_inband', 'k0_inband', 'lambda0_inband', ...
        'bandIdx', 'tStart', 'tEnd', 'nPeriods', 'x0', 'halfTr_inband');
end

% =========================================================================
% Part B: In-gap pulse
% =========================================================================
dataFileInGap = fullfile(outputDir, 'fig2_ingap_data.mat');

if doFDTD || ~exist(dataFileInGap, 'file')
    fprintf('\n=== Part B: In-gap FDTD simulation ===\n');

    kScan = linspace(0.01, 2.5, 2500);
    bands = temporal_crystal_bands(kScan, [eps1 eps2], [mu1 mu2], [t1 t2]);
    halfTr = bands.halfTrace;
    inGap = abs(halfTr) > 1;

    dG = diff([false, inGap, false]);
    gStarts = find(dG == 1);
    gEnds   = find(dG == -1) - 1;
    nGaps = length(gStarts);

    % Use 2nd gap (cleaner exponential growth)
    if nGaps >= 2, gapIdx = 2; else, gapIdx = 1; end
    kMidIdxGap = round((gStarts(gapIdx) + gEnds(gapIdx)) / 2);
    k0_ingap = kScan(kMidIdxGap);
    lambda0_ingap = 2*pi / k0_ingap;
    fprintf('  Gap %d: k0=%.4f, lambda=%.4f, |Tr/2|=%.3f\n', ...
        gapIdx, k0_ingap, lambda0_ingap, abs(halfTr(kMidIdxGap)));

    % Expected growth
    omegaF0 = bands.omegaF(:, kMidIdxGap);
    [~, growIdx] = max(abs(imag(omegaF0)));
    gammaFDTD = abs(imag(omegaF0(growIdx)));
    nPeriodsFDTD = round((340-220)/T);
    expectedGain = exp(gammaFDTD * nPeriodsFDTD * T);
    fprintf('  Expected gain over %d periods: e^{%.2f} ≈ %.1e\n', ...
        nPeriodsFDTD, gammaFDTD*nPeriodsFDTD*T, expectedGain);

    % FDTD grid (same as in-band)
    dx = lambda0_ingap / 28;
    xDomain = 140;
    x = 0:dx:xDomain;
    dt = 0.8 * dx;
    fprintf('  Grid: Nx=%d, dx=%.4f, dt=%.4f\n', numel(x), dx, dt);

    % Pulse
    x0 = xDomain * 0.22;
    sigma = 2.5;
    pulseFn = @(xq) exp(-((xq-x0)/sigma).^2) .* exp(1i*k0_ingap*(xq-x0));
    E0_ingap = pulseFn(x);
    xH = x(1:end-1) + dx/2;
    H0_ingap = nBg * pulseFn(xH + vBg*dt/2);

    tStart = 220;  tEnd = 340;
    nPeriods = round((tEnd - tStart) / T);
    tEnd = tStart + nPeriods * T;
    tSim = tEnd + 120;
    nSteps = ceil(tSim / dt);

    cfg = struct('x', x, 'dt', dt, 'nSteps', nSteps, ...
        'epsFun', @(xq,tq) ptcEps(xq,tq,epsBg,eps1,eps2,t1,t2,tStart,tEnd), ...
        'muFun', @(xq,tq) ones(size(xq)), ...
        'E0', E0_ingap, 'Hhalf0', H0_ingap, ...
        'boundary', 'sponge', 'spongeCells', 140, 'spongeStrength', 0.10, ...
        'recordEvery', 5, 'storeFields', true, 'progressBar', true);
    out_ingap = fdtd1d_db(cfg);
    halfTr_ingap = abs(halfTr(kMidIdxGap));
    save(dataFileInGap, 'out_ingap', 'k0_ingap', 'lambda0_ingap', ...
        'gapIdx', 'tStart', 'tEnd', 'nPeriods', 'x0', 'halfTr_ingap', 'expectedGain');
    fprintf('  Saved in-gap data.\n');
else
    fprintf('=== Loading saved in-gap data ===\n');
    load(dataFileInGap, 'out_ingap', 'k0_ingap', 'lambda0_ingap', ...
        'gapIdx', 'tStart', 'tEnd', 'nPeriods', 'x0', 'halfTr_ingap', 'expectedGain');
end

% =========================================================================
% Create Figure 2
% =========================================================================
fig = figure('Color', 'w', 'Position', [30 30 1400 600]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel (a): In-band ---
axA = nexttile(tl);
Eabs = abs(out_inband.E);
Eabs(Eabs < 1e-10) = 1e-10;
imagesc(axA, out_inband.x, out_inband.t, log10(Eabs));
set(axA, 'YDir', 'normal');
hold(axA, 'on');

yline(axA, tStart, 'w--', 'LineWidth', 1.5);
yline(axA, tEnd, 'w--', 'LineWidth', 1.5);
xL = axA.XLim;
fill(axA, [xL(1) xL(2) xL(2) xL(1)], [tStart tStart tEnd tEnd], ...
    'white', 'FaceAlpha', 0.10, 'EdgeColor', 'none');

cbar = colorbar(axA);
cbar.Label.String = 'log_{10}|D|';
cbar.Label.FontSize = 10;
xlabel(axA, 'Space  x', 'FontSize', 12);
ylabel(axA, 'Time  t', 'FontSize', 12);
title(axA, sprintf('(a)  Pulse in a band  (k_0=%.3f, band %d)', ...
    k0_inband, bandIdx), 'FontSize', 13, 'FontWeight', 'bold');
colormap(axA, jet(256));
set(axA, 'FontSize', 10);

% --- Panel (b): In-gap ---
axB = nexttile(tl);
Eabs = abs(out_ingap.E);
Eabs(Eabs < 1e-10) = 1e-10;
imagesc(axB, out_ingap.x, out_ingap.t, log10(Eabs));
set(axB, 'YDir', 'normal');
hold(axB, 'on');

yline(axB, tStart, 'w--', 'LineWidth', 1.5);
yline(axB, tEnd, 'w--', 'LineWidth', 1.5);
xL = axB.XLim;
fill(axB, [xL(1) xL(2) xL(2) xL(1)], [tStart tStart tEnd tEnd], ...
    'white', 'FaceAlpha', 0.10, 'EdgeColor', 'none');

cbar = colorbar(axB);
cbar.Label.String = 'log_{10}|D|';
cbar.Label.FontSize = 10;
xlabel(axB, 'Space  x', 'FontSize', 12);
ylabel(axB, 'Time  t', 'FontSize', 12);
title(axB, sprintf('(b)  Pulse in a bandgap  (k_0=%.3f, gap %d)', ...
    k0_ingap, gapIdx), 'FontSize', 13, 'FontWeight', 'bold');
colormap(axB, jet(256));
set(axB, 'FontSize', 10);

% Overall title
title(tl, sprintf(['Fig. 2: FDTD — amplitude of |D(x,t)|  ' ...
    '(\\epsilon_1=%d, \\epsilon_2=%d, T=%.4g, %d periods)'], ...
    eps1, eps2, T, nPeriods), 'FontSize', 13, 'FontWeight', 'bold');

% =========================================================================
% Save
% =========================================================================
outputFile = fullfile(outputDir, 'fig2_fdtd_simulations.png');
try
    exportgraphics(fig, outputFile, 'Resolution', 150);
catch
    print(fig, outputFile, '-dpng', '-r150');
end
fprintf('\nSaved: %s\n', outputFile);
end

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
