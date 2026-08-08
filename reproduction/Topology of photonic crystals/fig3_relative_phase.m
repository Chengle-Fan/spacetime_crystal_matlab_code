function fig3_relative_phase()
%FIG3_RELATIVE_PHASE  Reproduce Fig. 3 of Lustig et al. (arXiv:1803.08731).
%
%   The phase difference phi = arg(E^- / E^+) between forward-
%   propagating (time-refracted) and backward-propagating (time-reflected)
%   Floquet modes, plotted for each of the first six momentum band gaps
%   of the binary photonic time-crystal.
%
%   The sign of phi in each gap is dictated by the Zak phases of all
%   lower-lying bands via Eq. (6), proving that the relative phase
%   is a topological observable.
%
%   Method: exact transfer-matrix calculation using
%   temporal_finite_crystal_response (no FDTD required).
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
nPeriods = 19;   % finite-duration PTC periods

% =========================================================================
% Part 1: Band structure and gap identification
% =========================================================================
fprintf('=== Computing band structure ===\n');
kMin = 0.01;
kMax = 3.5;
Nk = 5000;
kScan = linspace(kMin, kMax, Nk);
bands = temporal_crystal_bands(kScan, [eps1 eps2], [mu1 mu2], [t1 t2]);
halfTr = bands.halfTrace;
inBand = abs(halfTr) <= 1;
inGap  = ~inBand;

% Find gap segments (no toolbox required)
dGap = diff([false, inGap, false]);
gapStarts = find(dGap == 1);
gapEnds   = find(dGap == -1) - 1;
nGaps = length(gapStarts);
nGapsToShow = min(6, nGaps);

fprintf('Found %d momentum gaps, will show first %d\n', nGaps, nGapsToShow);
for g = 1:nGapsToShow
    ks = gapStarts(g); ke = gapEnds(g);
    fprintf('  Gap %d: k ∈ [%.4f, %.4f], |halfTr|_max = %.3f, Nk = %d\n', ...
        g, kScan(ks), kScan(ke), ...
        max(abs(halfTr(ks:ke))), ke-ks+1);
end

% =========================================================================
% Part 2: Zak phases for bands
% =========================================================================
fprintf('\n=== Computing Zak phases ===\n');

% Find band segments
dBand = diff([false, inBand, false]);
bandStarts = find(dBand == 1);
bandEnds   = find(dBand == -1) - 1;
nBands = length(bandStarts);

zakPhases = zeros(1, nBands);

for b = 1:nBands
    ks = bandStarts(b); ke = bandEnds(b);
    NkBand = ke - ks + 1;
    if NkBand < 10
        zakPhases(b) = 0;
        fprintf('  Band %d: too narrow (%d pts), Zak ≈ 0\n', b, NkBand);
        continue;
    end

    % Subsample
    sampleStep = max(1, floor(NkBand / 120));
    idx = ks:sampleStep:ke;
    if idx(end) ~= ke, idx = [idx, ke]; end
    kSample = kScan(idx);
    Ns = length(kSample);

    rightStates = complex(zeros(2, Ns));
    leftStates  = complex(zeros(2, Ns));

    for ik = 1:Ns
        kVal = kSample(ik);
        U = temporal_crystal_monodromy(kVal, [eps1 eps2], [mu1 mu2], [t1 t2]);
        [V, D] = eig(U);
        lam = diag(D);
        omegaVals = 1i*log(lam)/T;
        [~, sortIdx] = sort(real(omegaVals));
        V = V(:, sortIdx);
        lam = lam(sortIdx);

        if ik == 1
            branchIdx = 1;
        else
            ov1 = abs(rightStates(:, ik-1)' * V(:, 1));
            ov2 = abs(rightStates(:, ik-1)' * V(:, 2));
            if ov1 >= ov2, branchIdx = 1; else, branchIdx = 2; end
        end
        rightStates(:, ik) = V(:, branchIdx);

        % Left eigenvector
        [Vt, Dt] = eig(U.');
        [~, idxL] = min(abs(diag(Dt) - lam(branchIdx)));
        leftStates(:, ik) = conj(Vt(:, idxL));
        leftStates(:, ik) = leftStates(:, ik) / norm(leftStates(:, ik));
    end

    % Wilson loop
    links = complex(zeros(1, Ns-1));
    for ik = 1:Ns-1
        ovlp = leftStates(:, ik)' * rightStates(:, ik+1);
        if abs(ovlp) > 1e-12
            links(ik) = ovlp / abs(ovlp);
        else
            links(ik) = 1;
        end
    end
    zak = -angle(prod(links));
    zak = mod(zak + pi, 2*pi) - pi;
    zakPhases(b) = zak;

    fprintf('  Band %d: Zak = %+6.3f rad (%+5.0f°) → %s  [%d pts]\n', ...
        b, zak, zak*180/pi, ...
        iif(abs(zak) < 0.15, '0 (trivial)', '\pi (topological)'), Ns);
end

% =========================================================================
% Part 3: TMM phase for each gap
% =========================================================================
fprintf('\n=== TMM relative phase for each gap ===\n');

phiData  = cell(1, nGapsToShow);
kGapData = cell(1, nGapsToShow);

for g = 1:nGapsToShow
    ks = gapStarts(g); ke = gapEnds(g);
    kGapVals = kScan(ks:ke);
    NkGap = length(kGapVals);

    % Sample gap for TMM (dense enough for smooth curves)
    kGapSample = linspace(kGapVals(1), kGapVals(end), max(80, NkGap));

    response = temporal_finite_crystal_response(kGapSample, ...
        [eps1 eps2], [mu1 mu2], [t1 t2], nPeriods, ...
        epsBg, muBg, [1; 0]);

    kGapData{g} = kGapSample;
    phiData{g} = response.relativePhase;

    % Sign at gap center
    kMid = (kGapVals(1) + kGapVals(end)) / 2;
    [~, idxMid] = min(abs(kGapSample - kMid));
    phiMid = response.relativePhase(idxMid);
    fprintf('  Gap %d: phi(k_mid=%.4f) = %+.4f rad → sgn = %+d\n', ...
        g, kMid, phiMid, sign(phiMid));
end

% =========================================================================
% Part 4: Analytical prediction via Eq. (6)
% =========================================================================
fprintf('\n=== Topological prediction (Eq. 6) ===\n');
% sgn(phi_s) = (-1)^s * (-1)^(s-1) * exp(i * sum_{m=1}^{s-1} theta_m^Zak)
%            = -1 * exp(i * sum_{m=1}^{s-1} theta_m^Zak)

fprintf('Gap | TMM sgn | Pred sgn | Zak sum | Match?\n');
fprintf('----|---------|----------|---------|-------\n');
for g = 1:nGapsToShow
    % s = gap index, l = number of band crossings below gap = g-1
    l = g - 1;
    if g-1 > 0
        zakSum = sum(zakPhases(1:g-1));
    else
        zakSum = 0;
    end
    predPhase = (-1)^g * (-1)^l * exp(1i * zakSum);
    predSign = sign(angle(predPhase));
    if abs(angle(predPhase)) < 0.05
        predSign = 0;  % near-zero → sign ambiguous
    end

    kMid = (kScan(gapStarts(g)) + kScan(gapEnds(g))) / 2;
    [~, idxMid] = min(abs(kGapData{g} - kMid));
    tmmSign = sign(phiData{g}(idxMid));

    matchStr = '---';
    if predSign ~= 0
        matchStr = iif(tmmSign == predSign, 'YES', 'NO');
    end

    fprintf('  %d  |    %+2d     |    %+2d     |  %+5.2f  | %s\n', ...
        g, tmmSign, predSign, zakSum, matchStr);
end

% =========================================================================
% Create Figure 3
% =========================================================================
fig = figure('Color', 'w', 'Position', [30 30 1350 800]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for g = 1:nGapsToShow
    ax = nexttile(tl);
    hold(ax, 'on');

    % Gap k-range
    ks = gapStarts(g); ke = gapEnds(g);
    kLeft  = kScan(ks);
    kRight = kScan(ke);

    % Shade gap region
    fill(ax, [kLeft kRight kRight kLeft], [-pi -pi pi pi], ...
        [0.9 0.9 0.92], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    % Plot TMM phase
    plot(ax, kGapData{g}, phiData{g}/pi, 'b-', 'LineWidth', 1.8);

    % Gap edge markers
    xline(ax, kLeft, 'r:', 'LineWidth', 0.8);
    xline(ax, kRight, 'r:', 'LineWidth', 0.8);

    % Zero line
    yline(ax, 0, 'k-', 'LineWidth', 0.5);

    % Mark sign at gap center
    kMid = (kLeft + kRight) / 2;
    [~, idxMid] = min(abs(kGapData{g} - kMid));
    phiMid = phiData{g}(idxMid);
    yPos = phiMid/pi + 0.15 * sign(phiMid);
    if abs(yPos) > 0.9, yPos = phiMid/pi - 0.15 * sign(phiMid); end
    text(ax, kMid, yPos, sprintf('sgn= %+d', sign(phiMid)), ...
        'FontSize', 11, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'w');

    xlabel(ax, 'k', 'FontSize', 11);
    ylabel(ax, '\phi / \pi', 'FontSize', 11);
    title(ax, sprintf('Gap %d  [k: %.3f–%.3f]', g, kLeft, kRight), ...
        'FontSize', 12, 'FontWeight', 'bold');

    ylim(ax, [-1.05, 1.05]);
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontSize', 10);
end

title(tl, sprintf(['Fig. 3: Relative phase \\phi = arg(E^-/E^+) ' ...
    'for first %d momentum gaps\n' ...
    '\\epsilon_1=%d, \\epsilon_2=%d, T=%.3g, ' ...
    'N=%d periods, \\epsilon_{bg}=%g'], ...
    nGapsToShow, eps1, eps2, T, nPeriods, epsBg), ...
    'FontSize', 13, 'FontWeight', 'bold');

% =========================================================================
% Save
% =========================================================================
outputDir = fullfile(fileparts(mfilename('fullpath')), 'output');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

outputFile = fullfile(outputDir, 'fig3_relative_phase.png');
try
    exportgraphics(fig, outputFile, 'Resolution', 200);
catch
    print(fig, outputFile, '-dpng', '-r200');
end
fprintf('\nSaved: %s\n', outputFile);

save(fullfile(outputDir, 'fig3_data.mat'), ...
    'kScan', 'halfTr', 'inBand', 'inGap', ...
    'gapStarts', 'gapEnds', 'bandStarts', 'bandEnds', ...
    'zakPhases', 'kGapData', 'phiData', 'nGapsToShow', ...
    'eps1', 'eps2', 'T', 't1', 't2', 'nPeriods', 'epsBg');
fprintf('Saved data to fig3_data.mat\n');
end

% =========================================================================
% Helpers
% =========================================================================
function result = iif(condition, trueVal, falseVal)
if condition, result = trueVal; else, result = falseVal; end
end
