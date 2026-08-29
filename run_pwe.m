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
%       quasifrequency vs k, colored by the central m=0 harmonic weight).
%
%   Scenario 2 - a binary photonic time crystal (eps 4 -> 1, 50% duty),
%       cross-checked against the exact 2x2 D/B monodromy transfer matrix
%       (tmm_bands).
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
nK1     = 121;    % number of k points swept across the first BZ

% ---- Scenario 2: binary photonic time crystal ----
epsA    = 4;      % permittivity for t in [0, dutyA*T)
epsB    = 1;      % permittivity for t in [dutyA*T, T)
dutyA   = 0.5;    % duty cycle of the high-eps layer
T2      = 1;      % temporal period
Mtime2  = 19;     % temporal harmonic truncation (converged)
Nspace2 = 0;      % no spatial modulation (pure time crystal)
kNorm   = linspace(0.02, 1.45, 181);   % k*T/(2*pi) = k/Omega sweep

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

kScan1 = linspace(0, g1/2, nK1);
cfg1.kScan = kScan1;
bands1 = pwe_bands(fourier1, cfg1, 'omega');    % fixed-k -> quasifrequency

S1 = bands1.S;                                  % S = (2N+1)*(2M+1)
fprintf('  eps(x,t) = %.1f*(1 + %.2f*cos(g*x)*cos(Omega*t))\n', eps0, delta);
fprintf('  g = %.4f, Omega = %.4f (Lambda = %g, T = %g)\n', g1, Omega1, ...
    Lambda1, T1);
fprintf('  truncation Nspace=%d, Mtime=%d -> S=%d harmonics -> %d bands\n', ...
    Nspace1, Mtime1, S1, 2*S1);
fprintf('  k sweep 0 .. g/2 over %d points (folded first BZ)\n', nK1);

% Folded first-BZ spectrum, colored by central (m=0) harmonic weight.
wMat1 = bands1.m0Weight.';                      % 2S x nK (align with omega)
kGrid1 = repmat(kScan1, 2*S1, 1);
x1 = kGrid1(:)/(g1/2);                          % k / (g/2)
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

% Match each of the two TMM branches to its closest (unused) PWE branch.
nK2 = numel(kNorm);
pweMatched    = complex(nan(2, nK2));
matchedWeight = nan(2, nK2);
for ik = 1:nK2
    folded = bands2.omega(:, ik);               % 2S x 1 (already folded)
    wgt    = bands2.m0Weight(ik, :);            % 1 x 2S
    valid  = isfinite(folded) & abs(imag(folded)) < 0.35*Omega2;
    ids = find(valid);
    used = false(size(ids));
    for branch = 1:2
        target = tmm2.omega(branch, ik);
        cost = abs(folded(ids) - target)/Omega2 ...
             + 2e-3*(1 - wgt(ids).');
        cost(used) = inf;
        [~, loc] = min(cost);
        chosen = ids(loc);
        used(loc) = true;
        pweMatched(branch, ik)    = folded(chosen);
        matchedWeight(branch, ik) = wgt(chosen);
    end
end

% Matching error and momentum-gap center(s) (from the imaginary panel).
err = abs(pweMatched - tmm2.omega)/Omega2;      % 2 x nK
medianErr = median(err(:));
imMag = max(abs(imag(tmm2.omega)), [], 1);      % 1 x nK (max over branches)

% Locate each momentum gap as a contiguous k-region where Im(w) is nonzero,
% and report its center (peak of |Im|) in units of k/Omega.
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

fprintf('  eps(t): %.1f -> %.1f, duty %.2f, T = %g, Omega = %.4f\n', ...
    epsA, epsB, dutyA, T2, Omega2);
fprintf('  Mtime=%d, Nspace=%d -> S=%d harmonics -> %d bands\n', ...
    Mtime2, Nspace2, bands2.S, 2*bands2.S);
fprintf('  median |PWE - TMM|/Omega = %.4e\n', medianErr);
for i = 1:numel(gapCenters)
    fprintf('  momentum-gap center %d at k/Omega = %.3f (max |Im(w)/Omega| = %.4f)\n', ...
        i, gapCenters(i), gapPeaks(i)/Omega2);
end

fig2 = figure('Color','w','Position',[130 130 1100 460]);
tl2 = tiledlayout(fig2, 1, 2, 'TileSpacing','compact');

ax2a = nexttile(tl2);
hold(ax2a,'on');
hTmmRe = plot(ax2a, kNorm, real(tmm2.omega.')/Omega2, 'k-', 'LineWidth',1.6);
kScat2 = repmat(kNorm, 2, 1);
hPweRe = scatter(ax2a, kScat2(:), real(pweMatched(:))/Omega2, 12, ...
    matchedWeight(:), 'filled');
xlabel(ax2a, 'k/\Omega');
ylabel(ax2a, 'Re(\omega_F)/\Omega');
title(ax2a, 'Principal Floquet bands');
grid(ax2a,'on'); box(ax2a,'on');
caxis(ax2a, [0 1]);
cb2a = colorbar(ax2a); cb2a.Label.String = 'PWE m=0 weight';
legend(ax2a, [hTmmRe(1), hPweRe], {'TMM (exact)','PWE (matched)'}, ...
    'Location','best','Box','off');

ax2b = nexttile(tl2);
hold(ax2b,'on');
hTmmIm = plot(ax2b, kNorm, imag(tmm2.omega.')/Omega2, 'k-', 'LineWidth',1.6);
hPweIm = plot(ax2b, kNorm, imag(pweMatched.')/Omega2, 'o', ...
    'Color',[0.82 0.15 0.15], 'MarkerSize',3);
xlabel(ax2b, 'k/\Omega');
ylabel(ax2b, 'Im(\omega_F)/\Omega');
title(ax2b, 'Momentum-gap growth/decay');
grid(ax2b,'on'); box(ax2b,'on');
yline(ax2b, 0, 'k:');
legend(ax2b, [hTmmIm(1), hPweIm(1)], {'TMM (exact)','PWE (matched)'}, ...
    'Location','best','Box','off');

title(tl2, 'Binary photonic time crystal: PWE dots vs TMM lines');

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
