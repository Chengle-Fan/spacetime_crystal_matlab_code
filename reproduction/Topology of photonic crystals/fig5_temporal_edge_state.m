function fig5_temporal_edge_state()
%FIG5_TEMPORAL_EDGE_STATE  Reproduce Fig. 5 of Lustig et al., Optica 5, 1390 (2018).
%
%   Three-panel figure showing temporal topological edge states:
%     (a) Schematic of two cascaded PTCs with different Zak phases
%         (PTC1: eps=[3,1], PTC2: eps=[1,3]), interface at t = 8T.
%     (b) Temporal domain wall mode — field amplitude as a function
%         of time, showing exponential growth toward the interface
%         in PTC1, exponential decay away from it in PTC2, creating
%         a localized temporal peak at the topological boundary.
%     (c) Schematic of smooth (non-step-like) time modulation
%         preserving the topological features.
%
%   Panel (b) can use either:
%     - Analytical domain-wall matching (fast, default)
%     - Full FDTD simulation (set doFDTD=true, slow)
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
doFDTD = false;   % set true for FDTD simulation of panel (b) (slow)
outputDir = fullfile(fileparts(mfilename('fullpath')), 'output');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

% =========================================================================
% PTC parameters
% =========================================================================
epsA = 3;   epsB = 1;
T    = 2*pi;
durations = [0.5, 0.5] * T;  % equal segments
nPeriodsLeft  = 8;   % PTC1 periods before interface
nPeriodsRight = 8;   % PTC2 periods after interface

% PTC1: eps = [epsA, epsB] = [3, 1]
% PTC2: eps = [epsB, epsA] = [1, 3]  (swapped → different topology)
epsLeft  = [epsA, epsB];
epsRight = [epsB, epsA];
mu       = [1, 1];

% =========================================================================
% Part 1: Find temporal domain wall mode (analytical)
% =========================================================================
fprintf('=== Temporal domain wall mode analysis ===\n');
fprintf('PTC1 (left):  eps = [%d, %d]  for %d periods\n', epsLeft, nPeriodsLeft);
fprintf('PTC2 (right): eps = [%d, %d]  for %d periods\n', epsRight, nPeriodsRight);
fprintf('Interface at t = %dT = %.4g\n', nPeriodsLeft, nPeriodsLeft*T);

% Scan k to find the common momentum gap and domain wall mode
kNormRange = [0.30, 0.80];  % k/(2*pi/T) range for gap search
kNorms = linspace(kNormRange(1), kNormRange(2), 2000);
kValues = kNorms * 2*pi/T;

mode = temporal_domain_wall_mode(kValues, ...
    epsLeft, mu, durations, ...
    epsRight, mu, durations, nPeriodsLeft, nPeriodsRight);

kBest = mode.k;
kBestNorm = kBest / (2*pi/T);
fprintf('\nDomain wall mode found:\n');
fprintf('  k = %.6f (k_norm = %.4f)\n', kBest, kBestNorm);
fprintf('  Left  multiplier  |lambda_grow| = %.4f\n', abs(mode.leftMultiplier));
fprintf('  Right multiplier  |lambda_decay| = %.4f\n', abs(mode.rightMultiplier));
fprintf('  Eigenspace mismatch = %.3e\n', mode.bestMismatch);

% =========================================================================
% Part 2: FDTD simulation (optional, for panel b)
% =========================================================================
fdtDataFile = fullfile(outputDir, 'fig5_fdtd_data.mat');

if doFDTD
    fprintf('\n=== FDTD simulation of temporal edge state ===\n');

    % FDTD grid parameters
    lambda0 = 2*pi / kBest;
    dx = lambda0 / 30;
    xDomain = 120;
    x = 0:dx:xDomain;
    dt = 0.8 * dx;
    fprintf('Grid: Nx=%d, dx=%.4f, dt=%.4f\n', numel(x), dx, dt);

    % Plane-wave-like initial condition (narrowband around kBest)
    x0 = xDomain * 0.25;
    sigma = 12.0;  % wide spatial envelope → narrow k-space
    pulseFn = @(xq) exp(-((xq-x0)/sigma).^2) .* exp(1i*kBest*(xq-x0));
    E0 = pulseFn(x);

    nBg = sqrt((epsA + epsB)/2);
    xH = x(1:end-1) + dx/2;
    vBg = 1/nBg;
    H0 = nBg * pulseFn(xH + vBg*dt/2);

    % Temporal schedule: free space → PTC1 → PTC2 → free space
    tFree1 = 120;                          % free propagation before PTC1
    tPTC1_start = tFree1;
    tPTC1_end   = tPTC1_start + nPeriodsLeft * T;   % interface at t = 8T
    tPTC2_end   = tPTC1_end + nPeriodsRight * T;
    tSim        = tPTC2_end + 80;
    nSteps      = ceil(tSim / dt);

    fprintf('Schedule: free [0,%.1f] → PTC1 [%.1f,%.1f] → PTC2 [%.1f,%.1f] → free\n', ...
        tFree1, tPTC1_start, tPTC1_end, tPTC1_end, tPTC2_end);
    fprintf('Total steps: %d, t_max: %.1f\n', nSteps, nSteps*dt);

    cfg = struct('x', x, 'dt', dt, 'nSteps', nSteps, ...
        'epsFun', @(xq,tq) ptcDomainWallEps(xq, tq, epsA, epsB, ...
            tPTC1_start, tPTC1_end, tPTC2_end, durations), ...
        'muFun', @(xq,tq) ones(size(xq)), ...
        'E0', E0, 'Hhalf0', H0, ...
        'boundary', 'sponge', 'spongeCells', 120, 'spongeStrength', 0.08, ...
        'recordEvery', 3, 'storeFields', true, 'progressBar', true);

    fprintf('Running FDTD...\n');
    out_fdtd = fdtd1d_db(cfg);

    save(fdtDataFile, 'out_fdtd', 'kBest', 'kBestNorm', ...
        'tPTC1_start', 'tPTC1_end', 'tPTC2_end', ...
        'nPeriodsLeft', 'nPeriodsRight', 'epsA', 'epsB', 'T', 'durations');
    fprintf('Saved FDTD data.\n');
else
    fprintf('\n=== Skipping FDTD (set doFDTD=true to run) ===\n');
    fprintf('Using analytical domain wall mode for panel (b).\n');
end

% =========================================================================
% Create Figure 5
% =========================================================================
fig = figure('Color', 'w', 'Position', [30 30 1400 500]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel (a): Two cascaded PTCs schematic -------------------------------
axA = nexttile(tl);

% Build epsilon profile
nShowLeft = nPeriodsLeft;
nShowRight = min(nPeriodsRight, 8);
tLeft  = linspace(0, nShowLeft*T, nShowLeft*200);
tRight = linspace(nShowLeft*T, (nShowLeft + nShowRight)*T, nShowRight*200);
tAll = [tLeft, tRight];
epsAll = zeros(size(tAll));

% PTC1: eps=[epsA, epsB] for t < nPeriodsLeft*T
for i = 1:length(tLeft)
    tMod = mod(tLeft(i), T);
    if tMod < durations(1)
        epsAll(i) = epsA;
    else
        epsAll(i) = epsB;
    end
end

% PTC2: eps=[epsB, epsA] for t >= nPeriodsLeft*T
offset = length(tLeft);
for i = 1:length(tRight)
    tMod = mod(tRight(i) - nShowLeft*T, T);
    if tMod < durations(1)
        epsAll(offset + i) = epsB;  % swapped!
    else
        epsAll(offset + i) = epsA;
    end
end

% Fill the PTC regions
fill(axA, [tLeft, nShowLeft*T, 0], [epsAll(1:length(tLeft)), 0, epsAll(1)], ...
    [0.7 0.85 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.35);
hold(axA, 'on');
fill(axA, [nShowLeft*T, tRight, (nShowLeft+nShowRight)*T, nShowLeft*T], ...
    [epsAll(offset+1), epsAll(offset+1:end), 0, 0], ...
    [1.0 0.8 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.35);

stairs(axA, tAll, epsAll, 'b-', 'LineWidth', 1.6);

% Interface marker
xline(axA, nShowLeft*T, 'r-', 'LineWidth', 2.2);
text(axA, nShowLeft*T + 0.5, epsA + 0.3, ...
    sprintf('Interface\n(t = %dT)', nPeriodsLeft), ...
    'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'top');

% Labels
text(axA, nShowLeft*T/2, epsA + 0.6, ...
    sprintf('PTC 1\n\\epsilon=[%d,%d]', epsA, epsB), ...
    'FontSize', 10, 'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold', 'Color', [0.1 0.3 0.6]);
text(axA, nShowLeft*T + nShowRight*T/2, epsA + 0.6, ...
    sprintf('PTC 2\n\\epsilon=[%d,%d]', epsB, epsA), ...
    'FontSize', 10, 'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold', 'Color', [0.7 0.2 0.1]);

ylim(axA, [0, epsA + 1.2]);
xlabel(axA, 'Time  t', 'FontSize', 12);
ylabel(axA, 'Permittivity  \epsilon(t)', 'FontSize', 12);
title(axA, '(a)  Two cascaded PTCs', 'FontSize', 13, 'FontWeight', 'bold');
set(axA, 'FontSize', 10);
box(axA, 'on');

% --- Panel (b): Temporal edge state --------------------------------------
axB = nexttile(tl);
hold(axB, 'on');

if doFDTD && exist(fdtDataFile, 'file')
    % FDTD-based: field amplitude at a fixed x
    load(fdtDataFile, 'out_fdtd', 'tPTC1_start', 'tPTC1_end', 'tPTC2_end');

    [~, probeIdx] = min(abs(out_fdtd.x - out_fdtd.x(round(end*0.3))));
    Damp = abs(out_fdtd.E(:, probeIdx));

    semilogy(axB, out_fdtd.t, Damp, 'b-', 'LineWidth', 1.3);
    xline(axB, tPTC1_start, 'k--', 'LineWidth', 1);
    xline(axB, tPTC1_end, 'r-', 'LineWidth', 2);
    xline(axB, tPTC2_end, 'k--', 'LineWidth', 1);

    legend(axB, '|D(t)|', 'PTC start', 'Interface', 'PTC end', ...
        'Location', 'best', 'FontSize', 8);
    title(axB, '(b)  FDTD: temporal edge state', ...
        'FontSize', 13, 'FontWeight', 'bold');
else
    % Analytical: plot the domain wall mode envelope
    cellIndex = mode.cellIndex;
    stateNorm = mode.stateNorm;

    % Plot the state norm (envelope) across temporal cells
    % Cell index 0 = interface
    semilogy(axB, cellIndex, stateNorm, 'o-', ...
        'Color', [0.12 0.42 0.78], 'LineWidth', 2.2, ...
        'MarkerFaceColor', [0.12 0.42 0.78], 'MarkerSize', 8);

    % Interface marker
    xline(axB, 0, 'r-', 'LineWidth', 2.2);

    % Exponential fit guides for visual reference
    leftIdx = cellIndex <= 0;
    rightIdx = cellIndex >= 0;
    if any(leftIdx) && any(rightIdx)
        % Fit exponential growth on the left
        leftVals = log(stateNorm(leftIdx));
        leftVals = leftVals(isfinite(leftVals));
        % Fit exponential decay on the right (first few cells before re-growth)
        rightFitIdx = rightIdx & cellIndex <= 4;
        if sum(rightFitIdx) >= 2
            rightVals = log(stateNorm(rightFitIdx));
            rightVals = rightVals(isfinite(rightVals));
        end
    end

    % Annotate
    text(axB, -nPeriodsLeft/2, stateNorm(1)*1.5, ...
        {'Exponential'; 'growth (PTC1)'}, ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'Color', [0.1 0.3 0.6]);
    text(axB, nPeriodsRight/2, stateNorm(end)*0.7, ...
        {'Re-growth'; '(PTC2)'}, ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'Color', [0.7 0.2 0.1]);
    text(axB, 0.5, stateNorm(mode.interfaceId)*1.3, ...
        'Peak at\ninterface', ...
        'FontSize', 9, 'Color', 'r', 'FontWeight', 'bold');

    legend(axB, '||[D,B]|| / interface', 'Interface (t=8T)', ...
        'Location', 'best', 'FontSize', 9);
    title(axB, sprintf('(b)  Domain-wall mode  (k=%.4f)', kBestNorm), ...
        'FontSize', 13, 'FontWeight', 'bold');
end

xlabel(axB, 'Temporal cell index  (0 = interface)', 'FontSize', 11);
ylabel(axB, 'Normalized field amplitude  (log scale)', 'FontSize', 11);
set(axB, 'FontSize', 10);
grid(axB, 'on'); box(axB, 'on');

% --- Panel (c): Smooth modulation schematic -------------------------------
axC = nexttile(tl);

% Build a smoothed version using raised-cosine transitions
nShowSmooth = 6;
tSmooth = linspace(0, nShowSmooth*T, nShowSmooth*300);
epsSmooth = zeros(size(tSmooth));

% Smooth transition parameter (fraction of T/2)
smoothFrac = 0.15;

for i = 1:length(tSmooth)
    tMod = mod(tSmooth(i), T);
    halfT = T/2;

    if tMod < halfT
        % Segment 1 (eps=epsA) with smooth transitions at boundaries
        if tMod < smoothFrac*halfT
            % Rising edge from epsB to epsA
            frac = tMod / (smoothFrac*halfT);
            epsSmooth(i) = epsB + (epsA - epsB) * (sin(frac*pi - pi/2) + 1) / 2;
        elseif tMod > halfT - smoothFrac*halfT
            % Falling edge from epsA to epsB
            frac = (tMod - (halfT - smoothFrac*halfT)) / (smoothFrac*halfT);
            epsSmooth(i) = epsA + (epsB - epsA) * (sin(frac*pi - pi/2) + 1) / 2;
        else
            epsSmooth(i) = epsA;
        end
    else
        % Segment 2 (eps=epsB) with smooth transitions
        tMod2 = tMod - halfT;
        if tMod2 < smoothFrac*halfT
            frac = tMod2 / (smoothFrac*halfT);
            epsSmooth(i) = epsA + (epsB - epsA) * (sin(frac*pi - pi/2) + 1) / 2;
        elseif tMod2 > halfT - smoothFrac*halfT
            frac = (tMod2 - (halfT - smoothFrac*halfT)) / (smoothFrac*halfT);
            epsSmooth(i) = epsB + (epsA - epsB) * (sin(frac*pi - pi/2) + 1) / 2;
        else
            epsSmooth(i) = epsB;
        end
    end
end

% Step-like reference (faint)
tStep = linspace(0, nShowSmooth*T, nShowSmooth*50);
epsStep = zeros(size(tStep));
for i = 1:length(tStep)
    tMod = mod(tStep(i), T);
    if tMod < durations(1)
        epsStep(i) = epsA;
    else
        epsStep(i) = epsB;
    end
end
plot(axC, tStep, epsStep, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
hold(axC, 'on');

% Smooth profile
fill(axC, [tSmooth, tSmooth(end), tSmooth(1)], ...
    [epsSmooth, epsA, epsA], ...
    [0.7 0.85 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.35);
plot(axC, tSmooth, epsSmooth, 'b-', 'LineWidth', 2.2);

ylim(axC, [epsB-0.2, epsA+0.5]);
xlabel(axC, 'Time  t', 'FontSize', 12);
ylabel(axC, 'Permittivity  \epsilon(t)', 'FontSize', 12);
title(axC, '(c)  Smooth time modulation', 'FontSize', 13, 'FontWeight', 'bold');

% Annotate
text(axC, T/4, epsA+0.2, sprintf('\\epsilon_1=%d', epsA), ...
    'FontSize', 10, 'HorizontalAlignment', 'center');
text(axC, 3*T/4, epsB-0.1, sprintf('\\epsilon_2=%d', epsB), ...
    'FontSize', 10, 'HorizontalAlignment', 'center');
legend(axC, 'Step-like (ref.)', 'Smooth \epsilon(t)', ...
    'Location', 'best', 'FontSize', 9);
set(axC, 'FontSize', 10);
grid(axC, 'on'); box(axC, 'on');

% Overall title
title(tl, sprintf(['Fig. 5: Temporal topological edge states ' ...
    '(\\epsilon_1=%d, \\epsilon_2=%d, interface at t=%dT)'], ...
    epsA, epsB, nPeriodsLeft), 'FontSize', 13, 'FontWeight', 'bold');

% =========================================================================
% Save
% =========================================================================
outputFile = fullfile(outputDir, 'fig5_temporal_edge_state.png');
try
    exportgraphics(fig, outputFile, 'Resolution', 200);
catch
    print(fig, outputFile, '-dpng', '-r200');
end
fprintf('\nSaved: %s\n', outputFile);

save(fullfile(outputDir, 'fig5_data.mat'), ...
    'mode', 'epsA', 'epsB', 'T', 'durations', ...
    'nPeriodsLeft', 'nPeriodsRight', 'kBest', 'kBestNorm', ...
    'epsLeft', 'epsRight', 'mu');
fprintf('Saved data to fig5_data.mat\n');
end

% =========================================================================
% Permittivity for domain-wall FDTD
% =========================================================================
function epsVal = ptcDomainWallEps(x, t, epsA, epsB, tPTC1s, tPTC1e, tPTC2e, durations)
T = sum(durations);
if t < tPTC1s || t > tPTC2e
    % Free space (background)
    epsVal = (epsA + epsB)/2 * ones(size(x));
elseif t < tPTC1e
    % PTC1: eps = [epsA, epsB]
    tLocal = t - tPTC1s;
    tPhase = mod(tLocal, T);
    if tPhase < durations(1)
        epsVal = epsA * ones(size(x));
    else
        epsVal = epsB * ones(size(x));
    end
else
    % PTC2: eps = [epsB, epsA] (swapped)
    tLocal = t - tPTC1e;
    tPhase = mod(tLocal, T);
    if tPhase < durations(1)
        epsVal = epsB * ones(size(x));  % swapped
    else
        epsVal = epsA * ones(size(x));  % swapped
    end
end
end
