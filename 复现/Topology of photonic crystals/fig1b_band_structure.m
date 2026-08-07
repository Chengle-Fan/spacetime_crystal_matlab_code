function fig1b_band_structure()
%FIG1B_BAND_STRUCTURE  Reproduce Fig. 1(b) of Lustig et al. (arXiv:1803.08731).
%
%   Floquet band structure of the binary photonic time-crystal.
%   The Floquet frequency omega_F (real in bands, complex in gaps)
%   is plotted as a function of the conserved momentum k.
%   Each band is labeled with its quantized Zak phase (0 or pi).
%
%   Parameters: eps1=3, eps2=1, t1=t2=0.5*T, T=2*pi (Omega=1).
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
T    = 2*pi;       % modulation period (Omega = 1)
t1   = 0.5*T;      % equal segments
t2   = 0.5*T;

% =========================================================================
% Compute band structure
% =========================================================================
kMax  = 2.5;
Nk    = 3000;
kGrid = linspace(1e-6, kMax, Nk);  % avoid k=0 exactly

bands = temporal_crystal_bands(kGrid, [eps1 eps2], [mu1 mu2], [t1 t2]);
omegaF = bands.omegaF;    % 2 x Nk complex
halfTr = bands.halfTrace; % 1 x Nk

Omega = 2*pi/T;  % = 1

% =========================================================================
% Identify band intervals (no toolbox required)
% =========================================================================
inBand = abs(halfTr) <= 1;
inGap  = ~inBand;

% Find contiguous band segments
d = diff([false, inBand, false]);
bandStarts = find(d == 1);   % indices into kGrid
bandEnds   = find(d == -1) - 1;
nBands = length(bandStarts);

fprintf('Detected %d bands in k ∈ [%.6f, %.3f]\n', nBands, kGrid(1), kMax);
for b = 1:min(nBands, 12)
    fprintf('  Band %d: k ∈ [%.4f, %.4f]  (Nk = %d)\n', ...
        b, kGrid(bandStarts(b)), kGrid(bandEnds(b)), ...
        bandEnds(b) - bandStarts(b) + 1);
end

% Find contiguous gap segments
dGap = diff([false, inGap, false]);
gapStarts = find(dGap == 1);
gapEnds   = find(dGap == -1) - 1;
nGaps = length(gapStarts);
fprintf('  (%d gaps detected)\n', nGaps);

% =========================================================================
% Compute Zak phases for each band
% =========================================================================
fprintf('\n--- Zak phase computation ---\n');
zakPhases = zeros(1, nBands);

for b = 1:nBands
    ks = bandStarts(b);
    ke = bandEnds(b);
    kBand = kGrid(ks:ke);
    NkBand = length(kBand);

    if NkBand < 10
        fprintf('  Band %d: too few points (%d), skipping\n', b, NkBand);
        zakPhases(b) = 0;
        continue;
    end

    % Subsample to ~150 points per band for efficiency
    sampleStep = max(1, floor(NkBand / 150));
    sampleIdx = 1:sampleStep:NkBand;
    if sampleIdx(end) ~= NkBand
        sampleIdx = [sampleIdx, NkBand];
    end
    kSample = kBand(sampleIdx);
    Ns = length(kSample);

    rightStates = complex(zeros(2, Ns));
    leftStates  = complex(zeros(2, Ns));
    trackedBandIdx = 1;  % which eigenvalue branch we track

    for ik = 1:Ns
        kVal = kSample(ik);
        U = temporal_crystal_monodromy(kVal, [eps1 eps2], [mu1 mu2], [t1 t2]);

        % Right eigenvectors
        [V, D] = eig(U);
        lam = diag(D);
        omegaVals = 1i*log(lam)/T;

        % Sort by real(omega_F) for consistent ordering
        [~, sortIdx] = sort(real(omegaVals));
        V = V(:, sortIdx);
        lam = lam(sortIdx);

        % Track the correct band by continuity
        if ik == 1
            % For the first k-point, determine which eigenvalue branch
            % belongs to this band by checking real(omega_F) proximity
            % to the mid-gap edges
            trackedBandIdx = 1;  % start with branch 1
            rightStates(:, ik) = V(:, trackedBandIdx);
        else
            % Track by maximum overlap with previous eigenstate
            ovlp1 = abs(rightStates(:, ik-1)' * V(:, 1));
            ovlp2 = abs(rightStates(:, ik-1)' * V(:, 2));
            if ovlp1 >= ovlp2
                trackedBandIdx = 1;
            else
                trackedBandIdx = 2;
            end
            rightStates(:, ik) = V(:, trackedBandIdx);
        end

        % Left eigenvector for the tracked eigenvalue
        leftStates(:, ik) = computeLeftEigenvector(U, lam(trackedBandIdx));
    end

    % Compute biorthogonal Zak phase = Wilson loop
    links = complex(zeros(1, Ns-1));
    for ik = 1:Ns-1
        ovlp = leftStates(:, ik)' * rightStates(:, ik+1);
        if abs(ovlp) < 1e-12
            links(ik) = 1;  % singular, but shouldn't happen in a band
        else
            links(ik) = ovlp / abs(ovlp);
        end
    end
    wilsonLoop = prod(links);
    zak = -angle(wilsonLoop);
    zak = mod(zak + pi, 2*pi) - pi;

    zakPhases(b) = zak;
    fprintf('  Band %d: Zak = %+.4f rad (%+.1f deg) → %s  [%d k-points]\n', ...
        b, zak, zak*180/pi, ...
        iif(abs(zak) < 0.15, '0 (trivial)', '\pi (topological)'), Ns);
end

% =========================================================================
% Create figure
% =========================================================================
fig = figure('Color', 'w', 'Position', [100 100 1050 600]);

hold on;

% Shade momentum gaps in gray
for g = 1:nGaps
    ks = kGrid(gapStarts(g));
    ke = kGrid(gapEnds(g));
    fill([ks ke ke ks], [-0.55 -0.55 0.55 0.55], ...
        [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.55);
end

% Plot Floquet bands (only where inBand is true)
re1 = real(omegaF(1, :));
re2 = real(omegaF(2, :));

% Band 1 (lower branch)
plot(kGrid(inBand), re1(inBand), 'b-', 'LineWidth', 1.8);
% Band 2 (upper branch)
plot(kGrid(inBand), re2(inBand), 'b-', 'LineWidth', 1.8);

% Add Zak phase labels to bands
for b = 1:nBands
    ks = bandStarts(b);
    ke = bandEnds(b);
    kMid = (kGrid(ks) + kGrid(ke)) / 2;

    % Find omega_F at band center (interpolate to get the right branch)
    [~, kmIdx] = min(abs(kGrid - kMid));
    % Determine which omega_F branch is inside the band
    if inBand(kmIdx)
        omegaMid = re1(kmIdx);
        % Check if band 2 also has this k in-band
        if abs(re2(kmIdx)) < abs(omegaMid) && inBand(kmIdx)
            omegaMid = re2(kmIdx);
        end
    else
        omegaMid = 0;
    end

    % Choose label based on Zak phase
    if abs(zakPhases(b)) < 0.15
        zakLabel = '0';
    elseif abs(abs(zakPhases(b)) - pi) < 0.2
        zakLabel = '\pi';
    else
        zakLabel = sprintf('%.2f', zakPhases(b));
    end

    % Place label slightly above the band
    yPos = omegaMid + 0.07;
    if yPos > 0.48, yPos = omegaMid - 0.10; end

    text(kMid, yPos, sprintf('Zak=%s', zakLabel), ...
        'FontSize', 10, 'FontWeight', 'bold', ...
        'Color', [0.8 0.15 0.15], ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', [1 1 0.88]);
end

% Formatting
xlabel('Momentum  k  (a.u.)', 'FontSize', 13);
ylabel('Floquet frequency  \omega_F / \Omega', 'FontSize', 13);
title(sprintf(['Fig. 1(b): Floquet bands of binary PTC  ' ...
    '(\\epsilon_1=%d, \\epsilon_2=%d, T=%.4g, \\Omega=1)'], ...
    eps1, eps2, T), 'FontSize', 14, 'FontWeight', 'bold');

ylim([-0.55, 0.55]);
xlim([0, kMax]);
set(gca, 'FontSize', 12);
grid on;
box on;

% Zone boundary lines
yline(0.5, 'k:', 'LineWidth', 0.7);
yline(-0.5, 'k:', 'LineWidth', 0.7);

% =========================================================================
% Save
% =========================================================================
outputDir = fullfile(fileparts(mfilename('fullpath')), 'output');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

outputFile = fullfile(outputDir, 'fig1b_band_structure.png');
try
    exportgraphics(fig, outputFile, 'Resolution', 200);
catch
    print(fig, outputFile, '-dpng', '-r200');
end
fprintf('\nSaved: %s\n', outputFile);

save(fullfile(outputDir, 'fig1b_data.mat'), ...
    'kGrid', 'omegaF', 'halfTr', 'inBand', 'inGap', ...
    'zakPhases', 'nBands', 'nGaps', 'bandStarts', 'bandEnds', ...
    'gapStarts', 'gapEnds', 'eps1', 'eps2', 'T', 't1', 't2');
fprintf('Saved data to fig1b_data.mat\n');
end

% =========================================================================
% Helper: compute left eigenvector for a given eigenvalue
% =========================================================================
function w = computeLeftEigenvector(U, lambda)
% Find left eigenvector: w' * U = lambda * w'
% This is the conjugate of the right eigenvector of U.'
[Vt, Dt] = eig(U.');
lamDiag = diag(Dt);
[~, idx] = min(abs(lamDiag - lambda));
w = conj(Vt(:, idx));
% Normalize
w = w / (norm(w) + eps);
end

% =========================================================================
% Helper: inline conditional
% =========================================================================
function result = iif(condition, trueVal, falseVal)
if condition
    result = trueVal;
else
    result = falseVal;
end
end
