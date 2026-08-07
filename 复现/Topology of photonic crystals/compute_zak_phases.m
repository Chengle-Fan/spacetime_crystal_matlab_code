function compute_zak_phases()
%COMPUTE_ZAK_PHASES  Compute and visualize the Zak phases of all Floquet bands.
%
%   This script performs a rigorous, self-contained computation of the
%   biorthogonal Zak phase for each band of the binary photonic time-crystal
%   described in Lustig et al., "Topology of photonic time-crystals,"
%   arXiv:1803.08731v1 (2018).
%
%   Key steps:
%     1. Floquet band structure via temporal_crystal_bands.
%     2. Contiguous band-interval detection in k-space.
%     3. Per-band eigenstate tracking using continuity of overlap.
%     4. Wilson-loop Zak phase from biorthogonal left/right eigenvectors
%        of the 2×2 monodromy matrix U(k).
%     5. Verification: sign of relative phase phi in each gap via Eq. (6).
%
%   Output:
%     - Console summary table of all Zak phases.
%     - Figure with band structure + Zak phase labels.
%     - Figure with Wilson-loop convergence diagnostics.
%     - output/compute_zak_phases.png, output/fig_zak_wilson_diagnostics.png
%     - output/compute_zak_phases.mat (full data).

% --- Path setup ---------------------------------------------------------
rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(rootDir, 'startup_stm.m'));

% =========================================================================
% 1.  Parameters (identical to the paper)
% =========================================================================
eps1 = 3;   eps2 = 1;          % permittivities
mu1  = 1;   mu2  = 1;          % permeabilities
T    = 2*pi;                   % modulation period  →  Omega = 1
t1   = 0.5 * T;                % segment durations (equal duty cycle)
t2   = 0.5 * T;
epsBg = (eps1 + eps2) / 2;     % background  (impedance-matched average)
muBg  = 1;

fprintf('==========  Zak Phase Computation  ==========\n');
fprintf('PTC:  eps = [%d, %d],  T = %.4g,  t1 = t2 = T/2\n\n', eps1, eps2, T);

% =========================================================================
% 2.  Band structure on a dense k-mesh
% =========================================================================
kMax  = 3.5;
Nk    = 5000;
kGrid = linspace(1e-6, kMax, Nk);   % avoid k = 0 exactly

bands  = temporal_crystal_bands(kGrid, [eps1 eps2], [mu1 mu2], [t1 t2]);
omegaF = bands.omegaF;              % 2 × Nk   complex Floquet frequencies
halfTr = bands.halfTrace;           % 1 × Nk   Tr(U)/2

Omega = 2*pi / T;   % = 1

% =========================================================================
% 3.  Contiguous band segments (|Tr(U)/2| <= 1)
% =========================================================================
inBand = abs(halfTr) <= 1;
inGap  = ~inBand;

dBand = diff([false, inBand, false]);
bandStarts = find(dBand == 1);       % indices into kGrid
bandEnds   = find(dBand == -1) - 1;
nBands = length(bandStarts);

dGap = diff([false, inGap, false]);
gapStarts = find(dGap == 1);
gapEnds   = find(dGap == -1) - 1;
nGaps = length(gapStarts);

fprintf('Detected %d bands and %d gaps in k ∈ [0, %.3g]\n\n', nBands, nGaps, kMax);

% =========================================================================
% 4.  Per-band biorthogonal Zak phase
% =========================================================================
fprintf('--- Band properties ---\n');
fprintf('%4s  %10s  %10s  %8s  %10s  %8s  %5s\n', ...
    'Band', 'k_min', 'k_max', 'Nk_raw', 'Zak(rad)', 'Zak(°)', 'Type');

zakPhases     = zeros(1, nBands);
bandMidK      = zeros(1, nBands);
bandWidths    = zeros(1, nBands);
wilsonLoops   = complex(zeros(1, nBands));
allLinks      = cell(1, nBands);     % store Wilson links for diagnostics

for b = 1:nBands
    ks = bandStarts(b);
    ke = bandEnds(b);
    kBand = kGrid(ks:ke);
    NkRaw = ke - ks + 1;
    bandMidK(b)   = (kGrid(ks) + kGrid(ke)) / 2;
    bandWidths(b) = kGrid(ke) - kGrid(ks);

    if NkRaw < 15
        % Too few points for reliable computation
        zakPhases(b) = 0;
        fprintf('%4d  %10.4f  %10.4f  %8d  %10s  %8s  %5s\n', ...
            b, kGrid(ks), kGrid(ke), NkRaw, '---', '---', 'narrow');
        continue;
    end

    % ---- 4a.  Subsample to ~120 points for efficiency ------------------
    sampleStep = max(1, floor(NkRaw / 120));
    idxSample  = ks:sampleStep:ke;
    if idxSample(end) ~= ke
        idxSample = [idxSample, ke];
    end
    kSample = kGrid(idxSample);
    Ns = length(kSample);

    % ---- 4b.  Build right & left eigenstates along the band ------------
    rightStates = complex(zeros(2, Ns));
    leftStates  = complex(zeros(2, Ns));
    trackedBranch = 1;   % which eigenvalue branch of U(k) we follow

    for ik = 1:Ns
        kVal = kSample(ik);
        U = temporal_crystal_monodromy(kVal, [eps1 eps2], [mu1 mu2], [t1 t2]);

        % Right eigenvectors
        [V, D] = eig(U);
        lamVals = diag(D);
        omegaVals = 1i * log(lamVals) / T;
        [~, order] = sort(real(omegaVals));
        V = V(:, order);
        lamVals = lamVals(order);

        % Track the correct branch via maximal overlap with previous state
        if ik == 1
            trackedBranch = 1;
        else
            ov1 = abs(rightStates(:, ik-1)' * V(:, 1));
            ov2 = abs(rightStates(:, ik-1)' * V(:, 2));
            if ov1 >= ov2
                trackedBranch = 1;
            else
                trackedBranch = 2;
            end
        end
        rightStates(:, ik) = V(:, trackedBranch);

        % Left eigenvector:  w' * U = lambda * w'
        % i.e. w is the conjugate of the right eigenvector of U.'
        [Vt, Dt] = eig(U.');
        [~, idxL] = min(abs(diag(Dt) - lamVals(trackedBranch)));
        w = conj(Vt(:, idxL));
        leftStates(:, ik) = w / (norm(w) + eps);
    end

    % ---- 4c.  Wilson loop (product of biorthogonal unit-modulus links) -
    links = complex(zeros(1, Ns - 1));
    minAbsLink = inf;
    for ik = 1:(Ns - 1)
        ovlp = leftStates(:, ik)' * rightStates(:, ik + 1);
        absOvlp = abs(ovlp);
        if absOvlp < 1e-14
            links(ik) = 1;        % singular case (should not happen in band)
        else
            links(ik) = ovlp / absOvlp;
        end
        if absOvlp < minAbsLink
            minAbsLink = absOvlp;
        end
    end
    wilsonLoop = prod(links);
    allLinks{b} = links;

    % ---- 4d.  Zak phase = -arg(Wilson loop), wrapped to [-pi, pi] ------
    zak = -angle(wilsonLoop);
    zak = mod(zak + pi, 2*pi) - pi;
    zakPhases(b)   = zak;
    wilsonLoops(b) = wilsonLoop;

    % Determine "type":  0 ↔ trivial,  π ↔ topological
    if abs(zak) < 0.15
        ztype = '0 (trivial)';
    elseif abs(abs(zak) - pi) < 0.25
        ztype = 'pi (topo)';
    else
        ztype = sprintf('%.2f', zak);
    end

    fprintf('%4d  %10.4f  %10.4f  %8d  %+10.4f  %+8.1f  %5s  [max|link|=%.4f]\n', ...
        b, kGrid(ks), kGrid(ke), NkRaw, zak, zak*180/pi, ztype, minAbsLink);
end

% =========================================================================
% 5.  Verify Eq. (6) — topological prediction of gap-phase sign
% =========================================================================
fprintf('\n--- Topological verification via Eq. (6) ---\n');
fprintf('sgn(phi_s) = (-1)^s * (-1)^{s-1} * exp(i * sum_{m=1}^{s-1} theta_m)\n\n');

nPeriods = 19;   % finite-crystal periods
fprintf('%4s  %12s  %12s  %10s  %6s\n', ...
    'Gap', 'phi TMM', 'phi pred', 'Zak sum', 'Match');

for g = 1:min(6, nGaps)
    % TMM reference phase at gap centre
    ks = gapStarts(g);  ke = gapEnds(g);
    kGapVals = kGrid(ks:ke);
    kMid = (kGapVals(1) + kGapVals(end)) / 2;

    % Compute TMM response at a few points around gap centre, interpolate
    kFine = linspace(kGapVals(1), kGapVals(end), 120);
    resp = temporal_finite_crystal_response(kFine, ...
        [eps1 eps2], [mu1 mu2], [t1 t2], nPeriods, epsBg, muBg, [1;0]);
    [~, imid] = min(abs(kFine - kMid));
    phiTMM = resp.relativePhase(imid);
    sgnTMM = sign(phiTMM);

    % Prediction from Eq. (6)
    l = g - 1;                                    % band crossings below gap
    if g > 1
        zakSum = sum(zakPhases(1:(g-1)));
    else
        zakSum = 0;
    end
    predComplex = (-1)^g * (-1)^l * exp(1i * zakSum);
    sgnPred = sign(angle(predComplex));
    if abs(angle(predComplex)) < 0.05, sgnPred = 0; end

    matchStr = '---';
    if sgnPred ~= 0
        matchStr = iif(sgnTMM == sgnPred, 'YES ✓', 'NO ✗');
    end

    fprintf('%4d  %+12.4f  %+12.4f  %+10.2f  %6s\n', ...
        g, phiTMM, angle(predComplex), zakSum, matchStr);
end

% =========================================================================
% 6.  Figure 1 — Band structure with Zak phase labels
% =========================================================================
fig1 = figure('Color','w','Position',[60 60 1100 600]);
hold on;

% Shade gaps
for g = 1:nGaps
    ks = kGrid(gapStarts(g));
    ke = kGrid(gapEnds(g));
    fill([ks ke ke ks], [-0.58 -0.58 0.58 0.58], ...
        [0.88 0.88 0.88], 'EdgeColor','none', 'FaceAlpha', 0.6);
end

% Bands (real part of omega_F, only in-band points)
re1 = real(omegaF(1,:));
re2 = real(omegaF(2,:));
plot(kGrid(inBand), re1(inBand), 'b-', 'LineWidth', 1.6);
plot(kGrid(inBand), re2(inBand), 'b-', 'LineWidth', 1.6);

% ---- Annotate each band with its Zak phase ----
for b = 1:nBands
    km = bandMidK(b);

    % Which omega_F branch is active at km?
    [~, ikm] = min(abs(kGrid - km));
    omegaMid = re1(ikm);
    % Distinguish the two branches
    if abs(re2(ikm)) < abs(re1(ikm))
        omegaMid = re2(ikm);
    end

    if abs(zakPhases(b)) < 0.15
        zakLab = '0';
    elseif abs(abs(zakPhases(b)) - pi) < 0.25
        zakLab = '\pi';
    else
        zakLab = sprintf('%.2f', zakPhases(b));
    end

    yOff = 0.07;
    yPos = omegaMid + yOff;
    if yPos > 0.47, yPos = omegaMid - 0.12; end

    text(km, yPos, sprintf('Zak=%s', zakLab), ...
        'FontSize', 10, 'FontWeight', 'bold', ...
        'Color', [0.85 0.15 0.15], ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', [1 1 0.85]);
end

yline( 0.5, 'k:', 'LineWidth', 0.6);
yline(-0.5, 'k:', 'LineWidth', 0.6);
xlabel('Momentum  k  (a.u.)', 'FontSize', 13);
ylabel('Floquet frequency  \omega_F / \Omega', 'FontSize', 13);
title(sprintf('Floquet bands with biorthogonal Zak phases  ' + ...
    '(\\epsilon_1=%d, \\epsilon_2=%d, T=%.4g)', eps1, eps2, T), ...
    'FontSize', 14, 'FontWeight', 'bold');
ylim([-0.58, 0.58]);
xlim([0, kMax]);
set(gca, 'FontSize', 12);
grid on; box on;

% =========================================================================
% 7.  Figure 2 — Wilson-loop convergence diagnostics
% =========================================================================
fig2 = figure('Color','w','Position',[100 100 1300 750]);
tl = tiledlayout(fig2, 2, 3, 'TileSpacing', 'compact');
nDiag = min(6, nBands);

for b = 1:nDiag
    ax = nexttile(tl);
    if isempty(allLinks{b})
        title(ax, sprintf('Band %d — too narrow', b));
        continue;
    end

    % Phase of each Wilson link (should vary smoothly)
    linkPhases = angle(allLinks{b});

    % Cumulative Wilson-loop phase
    cumPhase = -cumsum(linkPhases);   % -sum(arg(link_j)) → Zak phase as k→end

    yyaxis left;
    plot(ax, 1:length(linkPhases), linkPhases/pi, 'b.-', ...
        'MarkerSize', 6, 'LineWidth', 1.0);
    ylabel(ax, 'arg(link_j) / \pi', 'FontSize', 10);
    ylim(ax, [-1.2, 1.2]);
    yline(ax, 0, 'k-', 'LineWidth', 0.5);

    yyaxis right;
    plot(ax, 1:length(cumPhase), cumPhase/pi, 'r-', 'LineWidth', 1.8);
    ylabel(ax, '-\Sigma arg(link) / \pi', 'FontSize', 10);
    ylim(ax, [-1.2, 1.2]);

    zakLab = '?';
    if abs(zakPhases(b)) < 0.15
        zakLab = '0';
    elseif abs(abs(zakPhases(b)) - pi) < 0.25
        zakLab = '\pi';
    end
    title(ax, sprintf('Band %d — Zak = %s  (%.2f°)', ...
        b, zakLab, zakPhases(b)*180/pi), 'FontSize', 11, 'FontWeight', 'bold');
    xlabel(ax, 'k-index  j', 'FontSize', 9);
    grid(ax, 'on'); box(ax, 'on');
end

title(tl, 'Wilson-loop convergence diagnostics (per band)', ...
    'FontSize', 13, 'FontWeight', 'bold');

% =========================================================================
% 8.  Save outputs
% =========================================================================
outputDir = fullfile(fileparts(mfilename('fullpath')), 'output');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

% Figure 1
f1 = fullfile(outputDir, 'compute_zak_phases.png');
try
    exportgraphics(fig1, f1, 'Resolution', 220);
catch
    print(fig1, f1, '-dpng', '-r220');
end
fprintf('\nSaved: %s\n', f1);

% Figure 2
f2 = fullfile(outputDir, 'fig_zak_wilson_diagnostics.png');
try
    exportgraphics(fig2, f2, 'Resolution', 180);
catch
    print(fig2, f2, '-dpng', '-r180');
end
fprintf('Saved: %s\n', f2);

% Full data
dataFile = fullfile(outputDir, 'compute_zak_phases.mat');
save(dataFile, ...
    'kGrid', 'omegaF', 'halfTr', 'inBand', 'inGap', ...
    'bandStarts', 'bandEnds', 'gapStarts', 'gapEnds', ...
    'zakPhases', 'bandMidK', 'bandWidths', 'wilsonLoops', 'allLinks', ...
    'eps1', 'eps2', 'T', 't1', 't2', 'epsBg', 'muBg');
fprintf('Saved: %s\n', dataFile);

% =========================================================================
% 9.  Final summary table
% =========================================================================
fprintf('\n==========  Summary: Zak Phases of PTC Bands  ==========\n');
fprintf('Band |  k_range         |  Zak (rad)  |  Zak (°)  |  Type\n');
fprintf('-----|------------------|-------------|-----------|--------\n');
for b = 1:nBands
    if abs(zakPhases(b)) < 0.15
        ztype = '0  (trivial)';
    elseif abs(abs(zakPhases(b)) - pi) < 0.25
        ztype = 'pi (topological)';
    else
        ztype = sprintf('%.2f', zakPhases(b));
    end
    fprintf('  %2d  | [%.4f, %.4f]  |  %+9.4f  |  %+7.1f  |  %s\n', ...
        b, kGrid(bandStarts(b)), kGrid(bandEnds(b)), ...
        zakPhases(b), zakPhases(b)*180/pi, ztype);
end
fprintf('=========================================================\n');
end

% =========================================================================
% Helper
% =========================================================================
function result = iif(condition, trueVal, falseVal)
if condition, result = trueVal; else, result = falseVal; end
end
