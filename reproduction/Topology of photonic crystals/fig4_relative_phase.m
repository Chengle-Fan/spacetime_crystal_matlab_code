function fig4_relative_phase()
%FIG4_RELATIVE_PHASE  Reproduce Fig. 4 of Lustig et al., Optica 5, 1390 (2018).
%
%   Six panels (a)-(f) showing the relative phase phi = arg(E^- / E^+)
%   between forward-propagating (time-refracted) and backward-propagating
%   (time-reflected) Floquet modes for each of the first six momentum
%   band gaps of the binary PTC.
%
%   The sign of phi in each gap is dictated by the Zak phases of all
%   lower-lying bands via Eq. (6), proving that the relative phase
%   is a topological observable.
%
%   Method: exact transfer-matrix calculation using
%   temporal_finite_crystal_response.
%
%   Reference:
%     E. Lustig, Y. Sharabi, and M. Segev,
%     "Topological aspects of photonic time crystals,"
%     Optica 5, 1390-1395 (2018).  DOI: 10.1364/OPTICA.5.001390

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
nPeriods = 19;   % finite PTC duration

% =========================================================================
% Part 1: Band structure and gap identification
% =========================================================================
fprintf('=== Computing band structure and TMM relative phases ===\n');
kMin = 0.01;  kMax = 3.5;  Nk = 5000;
kScan = linspace(kMin, kMax, Nk);
bands = temporal_crystal_bands(kScan, [eps1 eps2], [mu1 mu2], [t1 t2]);
halfTr = bands.halfTrace;
inGap  = abs(halfTr) > 1;

% Find gap segments
dGap = diff([false, inGap, false]);
gapStarts = find(dGap == 1);
gapEnds   = find(dGap == -1) - 1;
nGaps = length(gapStarts);
nGapsToShow = min(6, nGaps);

fprintf('Found %d momentum gaps, displaying first %d:\n', nGaps, nGapsToShow);
for g = 1:nGapsToShow
    fprintf('  Gap %d: k ∈ [%.4f, %.4f]\n', ...
        g, kScan(gapStarts(g)), kScan(gapEnds(g)));
end

% =========================================================================
% Part 2: TMM relative phase for each gap
% =========================================================================
fprintf('\n=== TMM relative phase (phi = arg(E^- / E^+)) ===\n');

phiData  = cell(1, nGapsToShow);
kGapData = cell(1, nGapsToShow);

for g = 1:nGapsToShow
    ks = gapStarts(g);  ke = gapEnds(g);
    kGapVals = kScan(ks:ke);

    % Dense sampling within each gap for smooth curves
    kGapSample = linspace(kGapVals(1), kGapVals(end), max(80, length(kGapVals)));

    response = temporal_finite_crystal_response(kGapSample, ...
        [eps1 eps2], [mu1 mu2], [t1 t2], nPeriods, ...
        epsBg, muBg, [1; 0]);

    kGapData{g} = kGapSample;
    phiData{g}   = response.relativePhase;

    kMid = (kGapVals(1) + kGapVals(end)) / 2;
    [~, idxMid] = min(abs(kGapSample - kMid));
    phiMid = response.relativePhase(idxMid);
    fprintf('  Gap %d: phi(k_mid=%.4f) = %+.4f rad, sgn = %+d\n', ...
        g, kMid, phiMid, sign(phiMid));
end

% =========================================================================
% Part 3: Topological verification via Eq. (6)
% =========================================================================
fprintf('\n=== Topological verification (Eq. 6) ===\n');
fprintf('Eq. (6): sgn(phi_s) = delta * (-1)^s * (-1)^l * exp(i * sum Zak)\n');
fprintf('With delta = sgn(1-eps1/eps2) = sgn(%d) = %+d\n', 1-eps1/eps2, sign(1-eps1/eps2));

% Paper's Zak phases: alternating 0, pi, 0, pi, ...
% Derived from the gap sign pattern: (+1,+1,-1,-1,-1,+1) for first 6 gaps
delta = sign(1 - eps1/eps2);
paperSigns = [1, 1, -1, -1, -1, 1];  % from paper

fprintf('\nGap | TMM phi/pi | TMM sgn | Paper sgn | Match?\n');
fprintf('----|------------|---------|----------|-------\n');
for g = 1:nGapsToShow
    kMid = (kScan(gapStarts(g)) + kScan(gapEnds(g))) / 2;
    [~, idxMid] = min(abs(kGapData{g} - kMid));
    phiMid = phiData{g}(idxMid);
    tmmSign = sign(phiMid);
    paperSign = paperSigns(min(g, length(paperSigns)));
    matchStr = iif(tmmSign == paperSign, 'YES ✓', 'NO ✗');

    fprintf('  %d  | %+10.4f |    %+2d    |    %+2d     | %s\n', ...
        g, phiMid/pi, tmmSign, paperSign, matchStr);
end

% =========================================================================
% Create Figure 4
% =========================================================================
fig = figure('Color', 'w', 'Position', [30 30 1350 800]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Colors for the gaps (consistent palette)
gapColors = lines(nGapsToShow);

for g = 1:nGapsToShow
    ax = nexttile(tl);
    hold(ax, 'on');

    ks = gapStarts(g);  ke = gapEnds(g);
    kLeft  = kScan(ks);
    kRight = kScan(ke);

    % Shade gap region
    fill(ax, [kLeft kRight kRight kLeft], [-pi -pi pi pi], ...
        [0.9 0.9 0.92], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    % Plot TMM phase
    plot(ax, kGapData{g}, phiData{g}/pi, '-', ...
        'Color', gapColors(g, :), 'LineWidth', 2.0);

    % Gap edges
    xline(ax, kLeft, 'r:', 'LineWidth', 0.8);
    xline(ax, kRight, 'r:', 'LineWidth', 0.8);

    % Zero line
    yline(ax, 0, 'k-', 'LineWidth', 0.5);

    % Sign label at gap center
    kMid = (kLeft + kRight) / 2;
    [~, idxMid] = min(abs(kGapData{g} - kMid));
    phiMid = phiData{g}(idxMid);
    yPos = phiMid/pi + 0.15 * sign(phiMid);
    if abs(yPos) > 0.9, yPos = phiMid/pi - 0.15 * sign(phiMid); end
    text(ax, kMid, yPos, sprintf('sgn = %+d', sign(phiMid)), ...
        'FontSize', 12, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'w');

    xlabel(ax, 'Momentum  k', 'FontSize', 11);
    ylabel(ax, '\phi / \pi', 'FontSize', 11);
    title(ax, sprintf('(%c)  Gap %d', char('a'+g-1), g), ...
        'FontSize', 12, 'FontWeight', 'bold');

    ylim(ax, [-1.05, 1.05]);
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontSize', 10);
end

title(tl, sprintf(['Fig. 4: Relative phase \\phi = arg(E^-/E^+) ' ...
    'for the first %d momentum gaps\n' ...
    '\\epsilon_1=%d, \\epsilon_2=%d, N=%d periods'], ...
    nGapsToShow, eps1, eps2, nPeriods), ...
    'FontSize', 13, 'FontWeight', 'bold');

% =========================================================================
% Save
% =========================================================================
outputDir = fullfile(fileparts(mfilename('fullpath')), 'output');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

outputFile = fullfile(outputDir, 'fig4_relative_phase.png');
try
    exportgraphics(fig, outputFile, 'Resolution', 200);
catch
    print(fig, outputFile, '-dpng', '-r200');
end
fprintf('\nSaved: %s\n', outputFile);

save(fullfile(outputDir, 'fig4_data.mat'), ...
    'kScan', 'halfTr', 'inGap', ...
    'gapStarts', 'gapEnds', ...
    'kGapData', 'phiData', 'nGapsToShow', ...
    'eps1', 'eps2', 'T', 't1', 't2', 'nPeriods', 'epsBg');
fprintf('Saved data to fig4_data.mat\n');
end

% =========================================================================
function result = iif(condition, trueVal, falseVal)
if condition, result = trueVal; else, result = falseVal; end
end
