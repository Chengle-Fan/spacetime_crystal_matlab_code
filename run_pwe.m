%RUN_PWE  Plane-wave-expansion (PWE) entry script for 1D spacetime crystals.
%
% V3 conventions: 1D, scalar, non-dispersive media; normalized units
% (c0 = eps0 = mu0 = 1, mu_r = 1). Spatial period Lambda, temporal period T,
% g = 2*pi/Lambda, Omega = 2*pi/T, and
%     eps(x,t) = sum_{m,n} eps_mn * exp(i*n*g*x - i*m*Omega*t).
%
% Two scenarios are computed and plotted:
%
%   Scenario 1 - a general spacetime crystal
%       eps(x,t) = eps0*(1 + delta*cos(g*x).*cos(Omega*t)).
%       Shows the double-folded first-spacetime-BZ spectrum (Re/Im
%       quasifrequency vs k, colored by the central m=0 harmonic weight),
%       swept over the FULL first spatial BZ k in [-g/2, g/2) (any material,
%       including non-reciprocal travelling-wave modulations).
%
%   Scenario 2 - a binary photonic time crystal (eps 4 -> 1, 50% duty),
%       cross-checked against the exact 2x2 D/B monodromy transfer matrix
%       (tmm_bands). The two principal PWE branches are selected
%       INDEPENDENTLY of TMM (by central m=0 weight plus across-k continuity),
%       then labelled against TMM's two branches; the real part, the
%       growth/decay branch (in-gap) and the momentum-gap boundaries are each
%       compared separately and threshold asserts guard the agreement
%       (acceptance section 7). A uniform-medium analytic check
%       (omega = k/sqrt(eps*mu)) is asserted for BOTH solvers as well.
%
% Plan: define material -> period/truncation params -> compute -> plot.
%
% Results stay in the workspace. Nothing is written to disk and no directory
% is created unless the doSave guard at the very end is enabled.

%% ===================== GLOBAL PARAMETERS =====================
% Normalized units: c0 = eps0 = mu0 = 1.

% ---- Scenario 1: general spacetime crystal ----
eps0    = 2;      % mean permittivity
delta   = 0.3;    % modulation depth
Lambda1 = 1;      % spatial period
T1      = 1;      % temporal period
Nspace1 = 3;      % spatial harmonic truncation (n = -3:3)
Mtime1  = 3;      % temporal harmonic truncation (m = -3:3)
nK1     = 120;    % number of k points across the full first BZ (even)

% ---- Scenario 2: binary photonic time crystal ----
epsA    = 4;      % permittivity for t in [0, dutyA*T)
epsB    = 1;      % permittivity for t in [dutyA*T, T)
dutyA   = 0.5;    % duty cycle of the high-eps layer
T2      = 1;      % temporal period
Mtime2  = 19;     % temporal harmonic truncation (converged)
Nspace2 = 0;      % no spatial modulation (pure time crystal)
kNorm   = linspace(0.02, 1.45, 181);   % k*T/(2*pi) = k/Omega sweep
matchTol = 1e-2;  % PWE-vs-TMM median |Delta(omega)|/Omega threshold (assert)

% ---- Output guard (off by default; set true to export figures) ----
doSave  = false;

%% ===================== SCENARIO 1 =====================
fprintf('=== Scenario 1: general spacetime crystal ===\n');

g1     = 2*pi/Lambda1;
Omega1 = 2*pi/T1;
epsFun1 = @(X,Time) eps0*(1 + delta*cos(g1*X).*cos(Omega1*Time));

cfg1 = struct('Lambda', Lambda1, 'T', T1, ...
              'Nspace', Nspace1, 'Mtime', Mtime1);
fourier1 = pwe_fourier(epsFun1, [], cfg1);      % mu_r = 1 by default

% Full first spatial BZ, half-open: k in [-g/2, g/2) with k = 0 included
% (even nK1). A zero-symmetric sweep is required for general (non-reciprocal)
% media; a 0..g/2 half-zone would silently drop the negative-k half.
kScan1 = -g1/2 + (0:nK1-1)*g1/nK1;
cfg1.kScan = kScan1;
bands1 = pwe_bands(fourier1, cfg1, 'omega');    % fixed-k -> quasifrequency

S1 = bands1.S;                                  % S = (2N+1)*(2M+1)
fprintf('  eps(x,t) = %.1f*(1 + %.2f*cos(g*x)*cos(Omega*t))\n', eps0, delta);
fprintf('  g = %.4f, Omega = %.4f (Lambda = %g, T = %g)\n', g1, Omega1, ...
    Lambda1, T1);
fprintf('  truncation Nspace=%d, Mtime=%d -> S=%d harmonics -> %d bands\n', ...
    Nspace1, Mtime1, S1, 2*S1);
fprintf('  full first BZ k in [-g/2, g/2) over %d points\n', nK1);

% Folded first-BZ spectrum, colored by central (m=0) harmonic weight.
wMat1 = bands1.m0Weight.';                      % 2S x nK (align with omega)
kGrid1 = repmat(kScan1, 2*S1, 1);
x1 = kGrid1(:)/(g1/2);                          % k / (g/2), spans [-1,1)
yRe1 = real(bands1.omega(:))/Omega1;
yIm1 = imag(bands1.omega(:))/Omega1;
w1 = wMat1(:);

fig1 = figure('Color','w','Position',[100 100 1100 460]);
tl1 = tiledlayout(fig1, 1, 2, 'TileSpacing','compact');

ax1a = nexttile(tl1);
scatter(ax1a, x1, yRe1, 12, w1, 'filled');
xlabel(ax1a, 'k/(g/2)');
ylabel(ax1a, 'Re(\omega)/\Omega');
title(ax1a, 'Folded first-BZ bands');
grid(ax1a,'on'); box(ax1a,'on');
caxis(ax1a, [0 1]);
cb1a = colorbar(ax1a); cb1a.Label.String = 'm=0 weight';

ax1b = nexttile(tl1);
scatter(ax1b, x1, yIm1, 12, w1, 'filled');
xlabel(ax1b, 'k/(g/2)');
ylabel(ax1b, 'Im(\omega)/\Omega');
title(ax1b, 'Momentum/energy gaps');
grid(ax1b,'on'); box(ax1b,'on');
caxis(ax1b, [0 1]);
cb1b = colorbar(ax1b); cb1b.Label.String = 'm=0 weight';

title(tl1, sprintf(['General spacetime crystal ' ...
    'eps(x,t)=%.1f(1+%.1f cos(gx)cos(\\Omegat)): folded first-BZ spectrum'], ...
    eps0, delta));

%% ===================== SCENARIO 2 =====================
fprintf('=== Scenario 2: binary PTC, PWE vs TMM ===\n');

Omega2  = 2*pi/T2;
kValues = kNorm*Omega2;                         % kNorm = k/Omega

% Analytic binary Fourier coefficient (inline, no legacy dependency).
%   eps(t) = epsA on [0, dutyA*T), epsB on [dutyA*T, T), expanded as
%   eps(t) = sum_m eps_m exp(-i*m*Omega*t)  (matches the V3 n=0 convention).
epsCoeff2 = @(m,n) double(n==0)*binaryCoeff(m, epsA, epsB, dutyA);
muCoeff2  = @(m,n) double(m==0 && n==0);

% Minimal fourier struct consumed by pwe_bands (Lambda irrelevant: Nspace=0).
fourier2 = struct('Lambda', 1, 'T', T2, 'g', 2*pi/1, 'Omega', Omega2, ...
                  'Nspace', Nspace2, 'Mtime', Mtime2, ...
                  'epsCoeff', epsCoeff2, 'muCoeff', muCoeff2);

cfg2 = struct('Lambda', 1, 'T', T2, ...
              'Nspace', Nspace2, 'Mtime', Mtime2, 'kScan', kValues);
bands2 = pwe_bands(fourier2, cfg2, 'omega');    % folded quasifrequencies

tmm2 = tmm_bands(kValues, [epsA epsB], 1, [dutyA 1-dutyA]*T2);

% Independent principal-band selection: the two highest m=0-weight PWE
% branches at each k, with a nearest-by-continuity relabel across k so the
% two paths stay continuous. TMM is NOT used to select these branches.
[pweSelected, selectedWeight] = select_principal_bands(bands2, Omega2);

% Label the two selected PWE paths against the two TMM branches by global
% assignment (the lowest summed distance); this is labelling, not selection.
assignCostAB = sum(sum(abs(pweSelected - tmm2.omega).^2, 1), 2);
swapped = [tmm2.omega(2,:); tmm2.omega(1,:)];
assignCostBA = sum(sum(abs(pweSelected - swapped).^2, 1), 2);
if assignCostAB <= assignCostBA
    tmmLabels = tmm2.omega;
else
    tmmLabels = swapped;
end

% Separate comparisons: real part, growth/decay branch, gap boundaries.
% The folded quasifrequency is periodic in omega with period Omega (the
% first temporal BZ is [-Omega/2, Omega/2]); +pi and -pi are the SAME folded
% point, so the real-part error must be wrapped onto that circle (a gap
% branch may legitimately be returned at either seam).
dRe = mod(real(pweSelected) - real(tmmLabels) + Omega2/2, Omega2) - Omega2/2;
realErr = abs(dRe)/Omega2;
imErr   = abs(imag(pweSelected) - imag(tmmLabels))/Omega2;
weightMask = selectedWeight >= 0.3;             % only well-resolved branches
compMask   = weightMask & isfinite(realErr) & isfinite(imErr);
medianRealErr = median(realErr(compMask), 'all');
maxRealErr    = max(realErr(compMask), [], 'all');
medianImErr   = median(imErr(compMask), 'all');

% In-gap-only growth/decay error (acceptance 7.308: the growth/decay branch
% is compared separately). The all-column median is dominated by out-of-gap
% columns where Im ~ 0 for BOTH solvers and understates the in-gap
% disagreement, so the in-gap error is the asserted quantity.
%
% The in-gap comparison is pairing-INVARIANT. Inside a momentum gap the PWE
% eig solver decides at each gap ENTRY which branch is growth vs decay by
% floating-point rounding: all four folded copies of the two modes are
% exactly equidistant from the real previous pair, so the growth/decay signs
% are set by the solver's arithmetic, and the passband-driven global
% assignment (tmmLabels) cannot consistently align those signs either (it
% minimizes total distance, which is dominated by passband Re). Comparing
% signed Im would therefore fail or pass spuriously depending on the LAPACK
% build. Instead, at each in-gap column take the SMALLER of the two row
% pairings (identity vs swapped TMM) of the elementwise |Im| error, then max
% over columns. A conjugate PWE pair then reports its rate mismatch; a
% same-sign (both growth / both decay) selection cannot be paired down and
% still fails loudly.
imTolGap = 0.05;   % pre-registered: PWE growth-rate at Mtime = 19 is
                   % truncation-limited; measured max in-gap error ~2.3e-2
                   % = ~23% of the peak |Im|/Omega ~ 0.096, concentrated at
                   % the gap EDGE (k/Omega ~ 1.228), so 0.05 is a ~2x margin
imTmm = imag(tmm2.omega);                            % 2 x nK
inGapMask = (max(abs(imTmm), [], 1) > 1e-3*Omega2);
elemOk = weightMask & isfinite(imErr) & repmat(inGapMask, 2, 1);
inGapCols = find(inGapMask & all(elemOk, 1));
if isempty(inGapCols)
    error(['No in-gap PWE/TMM columns in the comparison set; the momentum-', ...
        'gap growth/decay comparison cannot be evaluated (acceptance 7.308).']);
end
imErrId = abs(imag(pweSelected(:, inGapCols)) - imTmm(:, inGapCols));
imErrSw = abs(imag(pweSelected(:, inGapCols)) - imTmm([2 1], inGapCols));
colBest = min(max(imErrId, [], 1), max(imErrSw, [], 1));  % per-column best pairing
gapRate = colBest/Omega2;
medianImErrGap = median(gapRate, 'all');
maxImErrGap    = max(gapRate, [], 'all');

% Momentum-gap boundaries and centers from the exact TMM half-trace.
imMag = max(abs(imag(tmm2.omega)), [], 1);      % 1 x nK (max over branches)
above = imMag > 1e-3*Omega2;
dd = diff([0 above 0]);
starts = find(dd == 1);
ends   = find(dd == -1) - 1;
gapCenters = zeros(1, numel(starts));
gapPeaks   = zeros(1, numel(starts));
for i = 1:numel(starts)
    seg = imMag(starts(i):ends(i));
    [gapPeaks(i), loc] = max(seg);
    gapCenters(i) = kNorm(starts(i) + loc - 1);
end

% PWE-side gap boundaries from the SELECTED branches' imaginary parts, then
% compare the PWE and TMM gap centers (acceptance 7.308: the gap boundaries
% are compared separately, not extracted TMM-only).
imPweSel = max(abs(imag(pweSelected)), [], 1);
pweAbove = imPweSel > 1e-3*Omega2;
dp2 = diff([0 pweAbove 0]);
pStarts = find(dp2 == 1); pEnds = find(dp2 == -1) - 1;
pweGapCenters = zeros(1, numel(pStarts));
for i = 1:numel(pStarts)
    seg = imPweSel(pStarts(i):pEnds(i));
    [~, loc] = max(seg);
    pweGapCenters(i) = kNorm(pStarts(i) + loc - 1);
end
gapBoundTol = 0.02;     % in k/Omega (scan spacing ~ 7.9e-3, ~2.5 bins)
if numel(pweGapCenters) == numel(gapCenters)
    maxGapCenterErr = max(abs(gapCenters - pweGapCenters));
else
    dAll = abs(pweGapCenters - gapCenters.');
    maxGapCenterErr = max([max(min(dAll, [], 1)), max(min(dAll, [], 2).')]);
end

fprintf('  eps(t): %.1f -> %.1f, duty %.2f, T = %g, Omega = %.4f\n', ...
    epsA, epsB, dutyA, T2, Omega2);
fprintf('  Mtime=%d, Nspace=%d -> S=%d harmonics -> %d bands\n', ...
    Mtime2, Nspace2, bands2.S, 2*bands2.S);
fprintf('  principal bands selected by m0 weight + continuity, %d/%d columns valid\n', ...
    nnz(compMask), numel(compMask));
fprintf('  real part:     median |Re(PWE)-Re(TMM)|/Omega = %.4e (<= %.3g required)\n', ...
    medianRealErr, matchTol);
fprintf('                 max    |Re(PWE)-Re(TMM)|/Omega = %.4e\n', maxRealErr);
fprintf('  growth/decay:  all-column median |Im(PWE)-Im(TMM)|/Omega = %.4e\n', ...
    medianImErr);
fprintf('                 IN-GAP max |Im| error/Omega = %.4e (<= %.2f required)\n', ...
    maxImErrGap, imTolGap);
assert(medianRealErr <= matchTol, ...
    ['median PWE-vs-TMM |Delta Re(omega)|/Omega = %.4e exceeds %.3g; ', ...
    'the principal-band match is not converged.'], medianRealErr, matchTol);
assert(maxImErrGap <= imTolGap, ...
    ['in-gap PWE-vs-TMM growth/decay |Delta Im(omega)|/Omega = %.4e exceeds ', ...
    'the pre-registered %.2f (acceptance 7.308).'], maxImErrGap, imTolGap);
for i = 1:numel(gapCenters)
    fprintf('  momentum-gap center %d at k/Omega = %.3f (max |Im(w)/Omega| = %.4f)\n', ...
        i, gapCenters(i), gapPeaks(i)/Omega2);
end
fprintf('  gap boundaries: TMM centers %s; PWE centers %s; max center diff %.4f k/Omega (<= %.2f)\n', ...
    mat2str(gapCenters, 3), mat2str(pweGapCenters, 3), maxGapCenterErr, gapBoundTol);
assert(maxGapCenterErr <= gapBoundTol, ...
    ['PWE-vs-TMM momentum-gap boundary center diff %.4f exceeds %.2f k/Omega ', ...
    '(acceptance 7.308).'], maxGapCenterErr, gapBoundTol);

%% ===================== UNIFORM-MEDIUM ANALYTIC CHECK =====================
% repair.md 7.308: uniform-medium PWE/TMM analytic error within tolerance.
% For eps = const, mu = 1 both solvers must recover the folded dispersion
% omega = k/sqrt(eps) to machine precision (previously claimed-but-not-
% enforced; now asserted in the entry).
epsU  = 2.25;  nU = sqrt(epsU);
kU    = linspace(0.05, 1.2, 101)*Omega2;   % physical k in the first two BZ
wU    = mod(kU/nU + Omega2/2, Omega2) - Omega2/2;   % analytic folded Re(omega)
fourierU = struct('Lambda', 1, 'T', T2, 'g', 2*pi, 'Omega', Omega2, ...
                  'Nspace', 0, 'Mtime', 3, ...
                  'epsCoeff', @(m,n) double(n==0 && m==0)*epsU, ...
                  'muCoeff',  @(m,n) double(m==0 && n==0));
bandsU = pwe_bands(fourierU, struct('Lambda',1,'T',T2,'Nspace',0,'Mtime',3, ...
    'kScan', kU), 'omega');
dPweU = mod(real(bandsU.omega) - wU + Omega2/2, Omega2) - Omega2/2;
minPweU = min(abs(dPweU), [], 1)/Omega2;
tmmU = tmm_bands(kU, [epsU epsU], 1, [T2/2 T2/2]);
dTmmU = mod(real(tmmU.omega) - wU + Omega2/2, Omega2) - Omega2/2;
minTmmU = min(abs(dTmmU), [], 1)/Omega2;
uniformTol = 1e-6;
maxErrPweU = max(minPweU);
maxErrTmmU = max(minTmmU);
assert(maxErrPweU <= uniformTol, ...
    ['uniform-medium PWE error %.3e exceeds %.1e (acceptance 7.308).'], ...
    maxErrPweU, uniformTol);
assert(maxErrTmmU <= uniformTol, ...
    ['uniform-medium TMM error %.3e exceeds %.1e (acceptance 7.308).'], ...
    maxErrTmmU, uniformTol);
fprintf('  uniform-medium check (eps = %.2f, mu = 1; analytic omega = k/%.2f folded): ', ...
    epsU, nU);
fprintf('PWE max |err|/Omega = %.3e, TMM max = %.3e (<= %.1e)\n', ...
    maxErrPweU, maxErrTmmU, uniformTol);

fig2 = figure('Color','w','Position',[130 130 1100 460]);
tl2 = tiledlayout(fig2, 1, 2, 'TileSpacing','compact');

ax2a = nexttile(tl2);
hold(ax2a,'on');
hTmmRe = plot(ax2a, kNorm, real(tmmLabels.')/Omega2, 'k-', 'LineWidth',1.6);
kScat2 = repmat(kNorm, 2, 1);
hPweRe = scatter(ax2a, kScat2(:), real(pweSelected(:))/Omega2, 12, ...
    selectedWeight(:), 'filled');
xlabel(ax2a, 'k/\Omega');
ylabel(ax2a, 'Re(\omega_F)/\Omega');
title(ax2a, 'Principal Floquet bands');
grid(ax2a,'on'); box(ax2a,'on');
caxis(ax2a, [0 1]);
cb2a = colorbar(ax2a); cb2a.Label.String = 'PWE m=0 weight';
legend(ax2a, [hTmmRe(1), hPweRe], {'TMM (exact)','PWE (selected)'}, ...
    'Location','best','Box','off');

ax2b = nexttile(tl2);
hold(ax2b,'on');
hTmmIm = plot(ax2b, kNorm, imag(tmmLabels.')/Omega2, 'k-', 'LineWidth',1.6);
hPweIm = plot(ax2b, kNorm, imag(pweSelected.')/Omega2, 'o', ...
    'Color',[0.82 0.15 0.15], 'MarkerSize',3);
xlabel(ax2b, 'k/\Omega');
ylabel(ax2b, 'Im(\omega_F)/\Omega');
title(ax2b, 'Momentum-gap growth/decay');
grid(ax2b,'on'); box(ax2b,'on');
yline(ax2b, 0, 'k:');
legend(ax2b, [hTmmIm(1), hPweIm(1)], {'TMM (exact)','PWE (selected)'}, ...
    'Location','best','Box','off');

title(tl2, 'Binary photonic time crystal: PWE dots vs TMM lines');

fprintf('Base MATLAB only (no Toolbox); declared minimum R2020a (untested on this machine), verified R2026a.\n');
fprintf('DONE.\n');

%% ===================== OPTIONAL SAVE (off by default) =====================
if doSave
    outDir = fullfile(pwd, 'output');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    save_fig(fig1, fullfile(outDir, 'run_pwe_scenario1_spacetime.png'));
    save_fig(fig2, fullfile(outDir, 'run_pwe_scenario2_ptc_vs_tmm.png'));
    fprintf('Saved figures to %s\n', outDir);
end

% -------------------------------------------------------------------------
function [pweSelected, selectedWeight] = select_principal_bands(bands, Omega)
%SELECT_PRINCIPAL_BANDS Independently select the two physical Floquet bands.
%   At each scan column the two PWE branches with the largest central m=0
%   harmonic weight are chosen; a nearest-by-continuity relabel then keeps
%   the two selected paths continuous across the scan (branch identity).
%   Returns:
%     pweSelected   : 2 x nK complex folded quasifrequencies of the two paths.
%     selectedWeight: 2 x nK m0 weights of the selected branches.
nK = numel(bands.kScan);
omegaFold = bands.omega;                    % 2S x nK (already folded)
weights   = bands.m0Weight.';               % 2S x nK
pweSelected   = complex(nan(2, nK));
selectedWeight = nan(2, nK);

for ik = 1:nK
    wgt = weights(:, ik);
    [~, idSort] = sort(wgt, 'descend');
    if ik == 1
        chosen = idSort(1:2);               % two highest m0 weight
    else
        % Among the top-4 candidates, pick the distinct pair whose folded
        % omega lies closest (wrapped) to the previous two selected values.
        cand = idSort(1:min(4, numel(idSort)));
        prev = pweSelected(:, ik-1);
        bestCost = inf; chosen = [nan nan];
        bestIsConjugate = false;
        tieTol = 1e-9*Omega;      % degenerate-gap-edge tie tolerance
        for c1 = 1:numel(cand)
            for c2 = 1:numel(cand)
                if c1 == c2, continue; end
                cost = fold_dist(omegaFold(cand(c1), ik), prev(1), Omega) ...
                     + fold_dist(omegaFold(cand(c2), ik), prev(2), Omega);
                if cost < bestCost - tieTol
                    bestCost = cost; chosen = [cand(c1), cand(c2)];
                    bestIsConjugate = is_gap_conjugate(cand(c1), cand(c2), ...
                        omegaFold(:, ik), Omega);
                elseif abs(cost - bestCost) <= tieTol && ~bestIsConjugate && ...
                        is_gap_conjugate(cand(c1), cand(c2), omegaFold(:, ik), Omega)
                    % At a momentum-gap ENTRY column the previous pair is real
                    % and the two modes' folded copies cross at Re=0 (or the
                    % seam), so the [growth,decay], [growth,growth], and
                    % [decay,decay] continuations have IDENTICAL fold_dist to
                    % ~1e-15*Omega -- only eig-solver rounding decides between
                    % them. Prefer the physical conjugate pair (one growth,
                    % one decay) so the in-gap Im signs are deterministic
                    % across LAPACK builds. The strict-min pair keeps its win
                    % when it is already conjugate; a tied conjugate rescues
                    % only when the strict min is not.
                    bestCost = cost; chosen = [cand(c1), cand(c2)];
                    bestIsConjugate = true;
                end
            end
        end
    end
    pweSelected(:, ik) = omegaFold(chosen, ik);
    selectedWeight(:, ik) = weights(chosen, ik);
end
end

% -------------------------------------------------------------------------
function d = fold_dist(a, b, Omega)
%FOLD_DIST Distance between two folded complex quasifrequencies: the real
%parts compared on the wrapped circle of period Omega plus the imag parts.
d = abs(mod(real(a) - real(b) + Omega/2, Omega) - Omega/2) ...
    + abs(imag(a) - imag(b));
end

% -------------------------------------------------------------------------
function tf = is_gap_conjugate(i1, i2, wCol, Omega)
%IS_GAP_CONJUGATE True if PWE branches i1,i2 form an in-gap conjugate pair.
%   Inside a momentum gap the two physical Floquet branches are a
%   time-reversed (conjugate) pair: opposite Im signs and equal |Im|. At a
%   gap ENTRY column the folded copies of the two modes cross, so the
%   [growth,decay], [growth,growth] and [decay,decay] continuations have
%   identical fold_dist to ~1e-15*Omega and only eig-solver rounding decides
%   between them; preferring the conjugate pair makes the in-gap Im signs
%   deterministic across LAPACK builds. Branches only qualify if
%   |Im| > 1e-6*Omega: passband branches have Im ~ 1e-15*Omega and must NOT
%   satisfy this test.
imTol = 1e-6*Omega;
w1 = wCol(i1); w2 = wCol(i2);
tf = abs(imag(w1)) > imTol && abs(imag(w2)) > imTol ...
     && imag(w1)*imag(w2) < 0 ...
     && abs(abs(imag(w1)) - abs(imag(w2))) <= imTol;
end

% -------------------------------------------------------------------------
function coeff = binaryCoeff(m, epsA, epsB, dutyA)
%BINARYCOEFF Analytic Fourier coefficient eps_m of a binary time crystal.
%   eps(t) = epsA for t in [0, dutyA*T), epsB for t in [dutyA*T, T), with
%   eps(t) = sum_m eps_m exp(-i*m*Omega*t).
if m == 0
    coeff = dutyA*epsA + (1 - dutyA)*epsB;
else
    coeff = (epsA - epsB)*(exp(1i*2*pi*m*dutyA) - 1)/(1i*2*pi*m);
end
end

% -------------------------------------------------------------------------
function save_fig(fig, filePath)
%SAVE_FIG Export a figure to PNG with a fallback if exportgraphics fails.
try
    exportgraphics(fig, filePath, 'Resolution', 220);
catch
    print(fig, filePath, '-dpng', '-r220');
end
end
