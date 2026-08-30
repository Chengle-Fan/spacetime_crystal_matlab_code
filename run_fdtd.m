%RUN_FDTD  FDTD + FFT-band entry script (V3.0.1).
%
% V3 conventions: 1D, scalar, non-dispersive media; normalized units
% (c0 = eps0 = mu0 = 1, mu_r = 1). Spatial period Lambda, temporal period T,
% g = 2*pi/Lambda, Omega = 2*pi/T. The FDTD kernel advances D and B on a
% Yee grid (leapfrog); a soft broadband source excites many k and omega, and
% a space-time FFT of the E(x,t) history sampled INSIDE the finite-sample
% measurement ROI gives the Floquet band structure (real spectrum ridges
% only; Im(omega) inside momentum/energy gaps is not recovered by the FFT).
%
% Acceptance observable: the ELECTRIC FIELD E(x,t) sampled at the integer
% Yee nodes inside the finite-sample ROI (P0-02). D is retained only as an
% optional diagnostic; it is never the acceptance observable.
%
% All three cases use a TRUE FINITE material sample (P0-01):
%     sponge | background | broadband source | background | SAMPLE |
%     background | sponge
% The SAMPLE is the only region where the material is modulated; outside it
% the background permittivity epsBg is constant. This produces real sample
% interfaces, boundary reflections and finite-size resonances. The outer
% layers are a simplified absorbing sponge, NOT a PML.
%
%   Case A - field distribution: a Gaussian wave packet launched in the left
%            background propagates into the finite modulated sample and is
%            plotted as |E|(x,t). Display only: no quantitative acceptance
%            in 3.0.1.
%   Case B - finite binary time-crystal sample: E-FFT band structure inside
%            the sample-interior measurement ROI, cross-checked against the
%            exact 2x2 D/B monodromy (tmm_bands). This is the ONLY
%            quantitative gate in 3.0.1, supplemented by a no-sample
%            reference run (incident-spectrum coverage + sponge residual
%            reflection) and a 1/L ridge extrapolation onto a common
%            signed-k grid over three sample lengths, compared with the
%            pre-registered tolerance derived from the measurement windows
%            (P0-01 / acceptance section 7).
%   Case C - a finite spacetime crystal eps(x,t) = eps0*(1 + delta*
%            cos(g*x)*cos(Omega*t)); the double-folded first-spacetime-BZ
%            spectrum with a PWE theory overlay. Generated WITHOUT
%            quantitative acceptance in 3.0.1. This default demo prints a
%            REPORT-ONLY independent ridge-vs-PWE comparison (max-power
%            folded bin per active column vs the light-line-closest PWE
%            fundamental band); the QUANTITATIVE cross-length (L=960/1280)
%            comparison, the two record-duration classes, and the two
%            ROI/window convergence classes run in the one-time regression
%            plan (repair.md section 7 / 6.5), not in this default demo.
%
% The Case B eps(t) profile is the WEAK contrast eps 1.3 <-> 1 (see the
% in-section note): a strong-contrast binary crystal (e.g. 4 <-> 1) opens a
% momentum gap whose parametric growth Im(omega)*T ~ 0.5 dominates the
% finite-time-window spectrum before the bulk bands can be resolved. The
% weak-contrast benchmark is a meaningful TMM cross-check of the
% finite-sample E-FFT method; the strong-contrast growth is a documented
% physical limitation, not a suppressed effect.
%
% Plan: 定义材料 -> 有限样品几何 -> 计算 -> 作图
%
% Results stay in the workspace. Nothing is written to disk and no directory
% is created unless the doSave guard at the very end is enabled. The entry
% script is the ONLY place allowed to overlay the explicit TMM/PWE theory
% curves.

clear; clc; close all;

%% ===================== GLOBAL PARAMETERS =====================
% Normalized units: c0 = eps0 = mu0 = 1. Spatially-uniform media have mu_r=1.

% ---- shared material: binary photonic time crystal ----
%   eps(t) = epsHi for t in [0, T/2),  epsLo for t in [T/2, T).
epsHi = 1.3;           % high permittivity (first half-period)
epsLo = 1.0;           % low  permittivity (second half-period)
mu_r  = 1;             % relative permeability (V3 default)
Tper  = 1;             % temporal period T
Omega = 2*pi/Tper;     % temporal modulation frequency Omega = 2*pi/T = 2*pi
epsBg = 1;             % background permittivity (constant OUTSIDE the sample)

% scalar eps(t) -> material handles. epsFun(x,t) must return an array the
% same size as x (the FDTD kernel evaluates it with a scalar time t).
epsScalar = @(t) epsHi + (epsLo - epsHi)*double(mod(t, Tper) >= Tper/2);
muFun  = @(x,t) mu_r .* ones(size(x));

% CFL: a conservative global upper bound 1/sqrt(min(eps)*min(mu)) = 1 covers
% every case (background eps=1, sample eps in [1, 1.3], mu=1). Passing it
% makes the CFL decision exact instead of audit-estimated (P1-02/P1-03).
maxWaveSpeed = 1.0;

% ---- output guard (off by default; set true to export figures) ----
doSave = false;

%% ==================== CASE A: FIELD DISTRIBUTION ====================
% Wave packet in the left background incident on a finite modulated sample.
fprintf('=== Case A: wave packet incident on a finite modulated sample ===\n');

dxA    = 1/40;         % spatial step (Lambda = 1)
dtA    = 0.5*dxA;      % time step (CFL = 0.5 with maxWaveSpeed = 1)
nSpongeA = 120;        % sponge width in cells (3 Lambda)
nBgLA    = 320;        % left background cells (8 Lambda)
nSampleA = 800;        % sample cells (20 Lambda) -- binary time crystal
nBgRA    = 320;        % right background cells (8 Lambda)
nLeftA   = nSpongeA + nBgLA;
NxA      = 2*nLeftA + nSampleA;
xA       = (0:NxA-1)*dxA;               % E grid (no repeated endpoint)
iSample0A = nLeftA + 1;
iSample1A = nLeftA + nSampleA;
xSampleA  = [xA(iSample0A) xA(iSample1A)];   % sample slab x-extent

% Background + modulated sample material (mask-based; sample is the only
% modulated region, giving real interfaces at xSampleA).
inSampleA = @(X) double(X >= xSampleA(1) & X <= xSampleA(2));
epsFunA   = @(X,t) epsBg + (epsScalar(t) - epsBg).*inSampleA(X);

% Gaussian wave packet in the background (t = 0 state, P3-06 corrected).
% v = 1/sqrt(epsBg) in the launch region; the H half-grid gets ONLY the time
% half-step phase (the Yee numerical-dispersion frequency), not an extra
% k0*dx/2 spatial factor (the H grid is already the spatial half-grid).
k0A      = 2.0*Omega;        % packet center wavenumber (pass band)
sigmaA   = 2.0;              % packet Gaussian width (Lambda)
xcA      = 6.0;              % packet center at t = 0 (inside the left bg)
vBgA     = 1/sqrt(epsBg*mu_r);
omegaYeeA = (2/dtA)*asin(min(1, vBgA*dtA/dxA*sin(k0A*dxA/2)));
E0A = exp(-((xA - xcA)/sigmaA).^2) .* exp(1i*k0A*(xA - xcA));
xHA = xA(1:end-1) + dxA/2;                 % H grid: Nx-1 points (sponge)
Hhalf0A = sqrt(epsBg/mu_r) * exp(1i*omegaYeeA*dtA/2) ...
          .* exp(-((xHA - xcA)/sigmaA).^2) .* exp(1i*k0A*(xHA - xcA));

nStepsA   = round(50*Tper/dtA);           % ~50 temporal periods
tEndA     = nStepsA*dtA;
% Abrupt permittivity switches must sit on integer time nodes (centred D/B
% correction): t = T/2, T, 3T/2, ... < tEndA.
interfaceTimesA = (Tper/2 : Tper/2 : tEndA - Tper/2);

fieldCfg = struct();
fieldCfg.x = xA;
fieldCfg.dt = dtA;
fieldCfg.nSteps = nStepsA;
fieldCfg.epsFun = epsFunA;
fieldCfg.muFun = muFun;
fieldCfg.E0 = E0A;
fieldCfg.Hhalf0 = Hhalf0A;
fieldCfg.boundary = 'sponge';
fieldCfg.spongeCells = nSpongeA;
fieldCfg.spongeStrength = 0.08;
fieldCfg.recordEvery = 2;
fieldCfg.temporalInterfaces = interfaceTimesA;
fieldCfg.maxWaveSpeed = maxWaveSpeed;

fieldResult = fdtd1d(fieldCfg);

% E is already recovered by the kernel (E = D/eps on the integer nodes); the
% after-switch take-side convention is documented in fdtd_fft_bands.
E_A = fieldResult.E;                  % nRecords x NxA
x_A = fieldResult.x(:).';
t_A = fieldResult.t(:);
[~, probeIdA] = min(abs(x_A - xcA));
probeE_A = abs(E_A(:, probeIdA));

figA = figure('Color','w','Position',[60 60 1160 440]);
tlA = tiledlayout(figA, 1, 2, 'TileSpacing','compact', 'Padding','compact');

axA1 = nexttile(tlA);
imagesc(axA1, x_A, t_A/Tper, abs(E_A));
axis(axA1, 'xy');
colormap(axA1, parula(256));
cbA1 = colorbar(axA1); cbA1.Label.String = '|E|';
hold(axA1, 'on');
% Mark the sample slab (solid) and the inner sponge edges (dashed).
xline(axA1, xSampleA(1), 'k-', 'sample', 'LineWidth', 1.4);
xline(axA1, xSampleA(2), 'k-', 'sample', 'LineWidth', 1.4);
xline(axA1, nSpongeA*dxA, 'r--', 'sponge', 'LineWidth', 1.0);
xline(axA1, (NxA-nSpongeA)*dxA, 'r--', 'sponge', 'LineWidth', 1.0);
xlabel(axA1, 'x');
ylabel(axA1, 't/T');
title(axA1, '|E| packet incident on finite sample');
grid(axA1, 'on'); box(axA1, 'on');

axA2 = nexttile(tlA);
plot(axA2, t_A/Tper, probeE_A, '-', 'Color',[0.00 0.45 0.75], 'LineWidth', 1.4);
xlabel(axA2, 't/T');
ylabel(axA2, sprintf('|E| at x/\\Lambda = %.2f', xcA));
title(axA2, 'Probe field vs time');
grid(axA2, 'on'); box(axA2, 'on');

title(tlA, sprintf(['Finite binary time crystal \\epsilon(t)=[%g %g] in ', ...
    'background \\epsilon=%g: packet at k=%.2f\\Omega'], epsHi, epsLo, ...
    epsBg, k0A/Omega));
drawnow;

fprintf('  sample slab x in [%.2f, %.2f] (%.1f Lambda), Nx=%d\n', ...
    xSampleA(1), xSampleA(2), xSampleA(2)-xSampleA(1), NxA);
fprintf('  dx=%.4f, dt=%.4f, CFL=%.3f, %d steps (~%.1f T)\n', ...
    dxA, dtA, maxWaveSpeed*dtA/dxA, nStepsA, nStepsA*dtA/Tper);

%% ==================== CASE B: FDTD-FFT BANDS VS TMM ====================
% Finite binary time-crystal sample; E-FFT ridges inside the sample-interior
% measurement ROI, cross-checked against the exact 2x2 D/B monodromy.
%
% Design notes (V3.0.1 fixes):
%   * long sample + sample-interior ROI (minus 10 cells per face) + rect
%     spatial window: a short central ROI cannot separate the bulk plane-wave
%     modes, so momentum-gap/band-edge content leaks into every k column. The
%     ROI must span most of a LONG sample.
%   * weak contrast eps 1.3 <-> 1: a strong-contrast binary crystal opens a
%     momentum gap whose parametric growth Im(omega)*T dominates the
%     finite-time-window spectrum before the bands resolve.
%   * 12-period Hann time window starting at the MEASURED arrival + 10
%     periods: the finite-time window is the resolution bottleneck; 12 periods
%     give bin spacing 2*pi/12 = 0.52 with Hann localization.
%   * narrow source (sx = 0.4): the incident spectrum covers k/Omega to ~1.4
%     at -40 dB, above the maxKNorm = 1.2 statistics range with margin
%     (verified by the no-sample reference below, repair.md line 309).
%   * Method A ridge: the ridge is the MAX-POWER bin over ALL omega*T per
%     active k column (theory-independent), and the error is
%     min(|ridge - exact|, |ridge + exact|) over the TMM pass band. The
%     real-field spectrum is omega-symmetric and band-edge leakage can pollute
%     either half, so no frequency half is assumed to carry the band.
fprintf('=== Case B: finite time crystal, E-FFT bands vs TMM ===\n');

caseBCfg = struct();
caseBCfg.dxB = 1/40;
caseBCfg.dtB = 0.4*caseBCfg.dxB;
caseBCfg.nSpongeB = 120;
caseBCfg.nBgLB = 160;
caseBCfg.nBgRB = 160;
caseBCfg.nPerB = 44;              % record length (temporal periods)
caseBCfg.recordEveryB = 2;       % 50 samples per period
caseBCfg.maxKNormB = 1.2;        % |k/Omega| statistics range (inside support)
caseBCfg.xcSrcB = 4.0;           % source center (Lambda), inside the left bg
caseBCfg.sxSrcB = 0.4;           % source spatial width (broad k)
caseBCfg.stSrcB = 0.2*Tper;      % source temporal width (broad omega)
caseBCfg.tcSrcB = 1.0*Tper;      % source peak time
caseBCfg.spongeStrengthB = 0.08;

% 1/L ridge-stability convergence over three sample lengths (30/40/60 Lambda).
sampleLengthsB = [1200 1600 2400];
benchB = struct('medianErrB',[], 'maxErrB',[], 'nPassB',[], 'noiseFloorB',[], ...
    'bandsB',[], 'kSelB',[], 'passbandB',[], 'exactOmegaTB',[], 'ridgeA',[], ...
    'nSampleB',[], 'nRoiB',[], 'tWinStartT',[], 'arrivalT',[]);
for liB = 1:numel(sampleLengthsB)
    benchB(liB) = run_caseB_benchmark(sampleLengthsB(liB), caseBCfg, ...
        epsScalar, muFun, epsBg, Tper, Omega, maxWaveSpeed);
end

% No-sample reference: incident-spectrum coverage + sponge residual reflection.
refB = run_caseB_reference(caseBCfg, sampleLengthsB(2), epsBg, muFun, ...
    Tper, maxWaveSpeed);

% ---- acceptance gate on the main (40 Lambda) benchmark ----
bMain = benchB(2);
% Pre-registered finite-sample tolerance (repair.md 7.311: "有限样品容差不能
% 机械照搬旧周期 0.10/0.30"). The Method-A ridge is the max-power bin on the
% zero-padded folded omega grid, so the peak-picking quantization error of a
% resolved ridge is bounded by that grid spacing (zeroPaddedGridSpacing =
% 2*pi/nFold). Median ~ half a folded bin (typical), max <= one folded bin
% (worst case for a resolved single peak). The effective window linewidth
% (ENBW x bin) and |vg|*dk_eff only broaden the peak envelope; they do not
% bias its max-bin location, so the folded-grid spacing is the honest bound.
zeroPadGridB = bMain.bandsB.zeroPaddedGridSpacingOmegaT;
tolMedianB = zeroPadGridB/2;
tolMaxB    = zeroPadGridB;
assert(refB.supportKNorm >= caseBCfg.maxKNormB, ...
    ['incident spectrum covers k/Omega only to %.2f at -40 dB; the %.2f ', ...
    'statistics range is not excited (repair.md line 309).'], ...
    refB.supportKNorm, caseBCfg.maxKNormB);
assert(refB.reflectRatio < 1e-2, ...
    'no-sample reference residual backscatter %.3e exceeds 1e-2', ...
    refB.reflectRatio);
assert(bMain.medianErrB <= tolMedianB, ...
    ['median FDTD-FFT vs TMM error %.4f exceeds the pre-registered %.4f ', ...
    '(zero-padded folded bin / 2, NOT the old periodic 0.10)'], ...
    bMain.medianErrB, tolMedianB);
assert(bMain.maxErrB <= tolMaxB, ...
    ['maximum FDTD-FFT vs TMM error %.4f exceeds the pre-registered %.4f ', ...
    '(one zero-padded folded bin, NOT the old periodic 0.30)'], ...
    bMain.maxErrB, tolMaxB);

% ---- robustness: target response invariant under bg-buffer / sponge change ----
% repair.md section 7 line 309: fix sample/source/ROI, only double the
% background buffer or change the sponge; the target response must stay within
% the pre-registered tolerance (both variants must still pass the gate and the
% median drift must stay within driftTol = tolMedianB below).
caseCfgBg = caseBCfg;
caseCfgBg.nBgLB = 2*caseBCfg.nBgLB;
caseCfgBg.nBgRB = 2*caseBCfg.nBgRB;
caseCfgSp = caseBCfg;
caseCfgSp.spongeStrengthB = 0.5*caseBCfg.spongeStrengthB;
variantBg = run_caseB_benchmark(sampleLengthsB(2), caseCfgBg, ...
    epsScalar, muFun, epsBg, Tper, Omega, maxWaveSpeed);
variantSp = run_caseB_benchmark(sampleLengthsB(2), caseCfgSp, ...
    epsScalar, muFun, epsBg, Tper, Omega, maxWaveSpeed);
driftBg = abs(variantBg.medianErrB - bMain.medianErrB);
driftSp = abs(variantSp.medianErrB - bMain.medianErrB);
% Both variants must pass the SAME pre-registered gate, and the median drift
% must stay within the ridge-quantization scale (driftTol = tolMedianB).
driftTol = tolMedianB;
assert(variantBg.medianErrB <= tolMedianB && variantBg.maxErrB <= tolMaxB, ...
    'bg-buffer-doubled variant fails the gate (%.4f/%.4f)', ...
    variantBg.medianErrB, variantBg.maxErrB);
assert(variantSp.medianErrB <= tolMedianB && variantSp.maxErrB <= tolMaxB, ...
    'sponge-halved variant fails the gate (%.4f/%.4f)', ...
    variantSp.medianErrB, variantSp.maxErrB);
assert(driftBg <= driftTol && driftSp <= driftTol, ...
    ['target-response drift exceeds the pre-registered %.4f (bg: %.4f, ', ...
    'sponge: %.4f)'], driftTol, driftBg, driftSp);

% ---- 1/L ridge extrapolation onto a common signed-k grid (repair.md 7.311) ----
% Each sample length has its own ROI length, hence its own k grid; the
% per-column Method-A ridge error is interpolated onto a common grid, then
% e(k) = a(k) + b(k)/L is fit by least squares and the L -> inf residual
% a(k) (the large-sample ridge-vs-TMM agreement) is compared with the
% pre-registered tolerance. This is a genuine cross-length continuity +
% 1/L extrapolation, not a raw-median stability plot. The fixed-central-ROI
% convergence class and the two-record-duration class are exercised in the
% one-time regression battery (not shipped, per the flat-tree constraint);
% the expanding-ROI series shipped here is the class that is physically
% meaningful for bulk-mode separation (see the design notes above).
exB = run_caseB_1L_extrapolation(benchB, sampleLengthsB, caseBCfg.maxKNormB, ...
    epsHi, epsLo, Tper, Omega);
fprintf('  1/L extrapolation (common signed-k grid): L->inf residual median = %.4f, max = %.4f (1/L slope median = %.4f)\n', ...
    exB.eInfMed, exB.eInfMax, exB.slopeMed);
assert(exB.eInfMax <= tolMaxB, ...
    ['1/L extrapolated ridge residual max %.4f exceeds the pre-registered ', ...
    '%.4f; the finite-sample ridge does not converge to the TMM band as ', ...
    'L -> inf (repair.md 7.311).'], exB.eInfMax, tolMaxB);
assert(exB.eInfMed <= tolMedianB, ...
    ['1/L extrapolated ridge residual median %.4f exceeds the pre-registered ', ...
    '%.4f (repair.md 7.311).'], exB.eInfMed, tolMedianB);

% ---- Case B figure: spectrum + ridge overlay and 1/L convergence ----
bandsB = bMain.bandsB;
kN_B    = bandsB.kNormalized(:).';
kSelB   = bMain.kSelB;
exactOT = bMain.exactOmegaTB;
passBB  = bMain.passbandB;
ridgeA  = bMain.ridgeA;

figB = figure('Color','w','Position',[70 70 1180 640]);
tlB = tiledlayout(figB, 1, 2, 'TileSpacing','compact', 'Padding','compact');

% --- left tile: E-FFT spectrum + TMM pass band + Method A ridge markers ---
axB1 = nexttile(tlB);
imagesc(axB1, kN_B, bandsB.omegaT(:), bandsB.spectralDb);
axis(axB1, 'xy');
colormap(axB1, parula(256));
cbB = colorbar(axB1); cbB.Label.String = 'folded spectral power [dB]';
hold(axB1, 'on');
% TMM pass-band curves (solid) over the pass-band segments only.
kSelBNum = kN_B(kSelB);
d = diff([false passBB false]);
pbStart = find(d == 1); pbStop = find(d == -1) - 1;
for s = 1:numel(pbStart)
    ids = pbStart(s):pbStop(s);
    plot(axB1, kSelBNum(ids),  exactOT(ids), 'k-', 'LineWidth', 1.4, ...
        'HandleVisibility','off');
    plot(axB1, kSelBNum(ids), -exactOT(ids), 'k-', 'LineWidth', 1.4, ...
        'HandleVisibility','off');
end
% Method A ridge markers (max-power bin over all omega*T), thinned.
ridgeSel = ridgeA(kSelB);
mkIds = find(passBB);
mkIds = mkIds(1:2:end);
hRidge = scatter(axB1, kSelBNum(mkIds), ridgeSel(mkIds), 14, 'w', 'filled');
hTheory = plot(axB1, nan, nan, 'k-', 'LineWidth', 1.4, 'DisplayName','TMM');
xlabel(axB1, 'k/\Omega   (\Omega = 2\pi/T)');
ylabel(axB1, '\omega T');
title(axB1, sprintf('E-FFT ridges vs TMM (median |\\Delta\\omega T| = %.3f)', ...
    bMain.medianErrB));
xlim(axB1, [0 caseBCfg.maxKNormB]); ylim(axB1, [-pi pi]);
yticks(axB1, [-pi 0 pi]); yticklabels(axB1, {'-\pi','0','\pi'});
legend(axB1, [hRidge hTheory], {'E-FFT ridge (ROI)','TMM pass band'}, ...
    'Location','northeast','Box','off');
grid(axB1, 'on'); box(axB1, 'on');

% --- right tile: 1/L ridge error extrapolation (least-squares in 1/L) ---
axB2 = nexttile(tlB);
medB = [benchB.medianErrB];
maxB = [benchB.maxErrB];
invL = 1./sampleLengthsB;
cBlue = [0.00 0.45 0.75]; cOrange = [0.85 0.33 0.10];
plot(axB2, invL, medB, '-o', 'Color',cBlue, ...
    'LineWidth', 1.6, 'MarkerFaceColor',cBlue, 'DisplayName','median');
hold(axB2, 'on');
plot(axB2, invL, maxB, '-s', 'Color',cOrange, ...
    'LineWidth', 1.6, 'MarkerFaceColor',cOrange, 'DisplayName','max');
% least-squares fits in 1/L and the extrapolated L -> inf intercepts.
plot(axB2, [0 invL], polyval(exB.medFit, [0 invL]), '--', 'Color',cBlue, ...
    'LineWidth', 1.0, 'HandleVisibility','off');
plot(axB2, [0 invL], polyval(exB.maxFit, [0 invL]), '--', 'Color',cOrange, ...
    'LineWidth', 1.0, 'HandleVisibility','off');
scatter(axB2, 0, exB.medIntercept, 50, 'o', 'MarkerFaceColor',cBlue, ...
    'MarkerEdgeColor','k', 'DisplayName','median (L->inf)');
scatter(axB2, 0, exB.maxIntercept, 50, 's', 'MarkerFaceColor',cOrange, ...
    'MarkerEdgeColor','k', 'DisplayName','max (L->inf)');
% pre-registered tolerance (zero-padded folded grid, not old 0.10/0.30).
yline(axB2, tolMedianB, ':', 'pre-reg. tol/2', 'Color',[0.35 0.35 0.35], ...
    'HandleVisibility','off');
yline(axB2, tolMaxB, ':', 'pre-reg. tol', 'Color',[0.35 0.35 0.35], ...
    'HandleVisibility','off');
xlabel(axB2, '1/L  (1 / sample cells)');
ylabel(axB2, '|Delta(\omega T)| vs TMM');
legend(axB2, 'Location','northeast','Box','off');
grid(axB2, 'on'); box(axB2, 'on');
title(axB2, '1/L ridge error extrapolation (L->inf)');

title(tlB, sprintf(['Finite binary time crystal \\epsilon(t)=[%g %g] in ', ...
    '\\epsilon_{bg}=%g: E-FFT ROI ridges vs TMM pass band'], ...
    epsHi, epsLo, epsBg));
drawnow;

fprintf('  sample lengths (Lambda cells): %s\n', mat2str(sampleLengthsB));
fprintf('  main benchmark (%d cells, ROI %d cells): ', ...
    bMain.nSampleB, bMain.nRoiB);
fprintf('median |d(omega*T)| = %.4f, max = %.4f, %d pass-band columns\n', ...
    bMain.medianErrB, bMain.maxErrB, bMain.nPassB);
fprintf('  1/L series median: %s;  max: %s;  L->inf intercepts median %.4f, max %.4f\n', ...
    mat2str(medB, 3), mat2str(maxB, 3), exB.medIntercept, exB.maxIntercept);
fprintf('  no-sample reference: incident k/Omega support -40dB = %.2f (>= %.2f); ', ...
    refB.supportKNorm, caseBCfg.maxKNormB);
fprintf('sponge residual backscatter = %.2e of incident peak\n', ...
    refB.reflectRatio);

%% ==================== CASE C: DOUBLE-FOLDED BZ ====================
% Finite spacetime-crystal sample; double-folded E-FFT over the central ROI.
fprintf('=== Case C: finite spacetime crystal, double-folded BZ ===\n');

eps0C  = 2;             % mean permittivity
deltaC = 0.3;           % modulation depth
LambdaC = 1;            % spatial period
TC      = 1;            % temporal period
gC      = 2*pi/LambdaC;
OmegaC  = 2*pi/TC;
epsFunC_smooth = @(X,T) eps0C*(1 + deltaC*cos(gC*X).*cos(OmegaC*T));

dxC    = LambdaC/40;    % spatial step (40 points per cell)
dtC    = 0.4*dxC;       % time step (CFL = 0.4 with maxWaveSpeed = 1)
nSpongeC = 120;         % sponge width in cells (3 Lambda)
nBgLC    = 160;         % left background cells (4 Lambda)
nSampleC = 960;         % sample cells (24 Lambda, integer unit cells)
nBgRC    = 160;         % right background cells (4 Lambda)
nLeftC   = nSpongeC + nBgLC;
NxC      = 2*nLeftC + nSampleC;
xC       = (0:NxC-1)*dxC;             % E grid (no repeated endpoint)
iSample0C = nLeftC + 1;
iSample1C = nLeftC + nSampleC;
xSampleC  = [xC(iSample0C) xC(iSample1C)];
% Sample boundaries fall on integer Lambda (nLeftC is a multiple of 40), so
% the termination/modulation phase is fixed and grid-aligned.
inSampleC = @(X) double(X >= xSampleC(1) & X <= xSampleC(2));
epsFunC   = @(X,t) epsBg + (epsFunC_smooth(X,t) - epsBg).*inSampleC(X);

% Central ROI spanning an integer number of unit cells (12 cells) for the
% spatial Floquet fold (P1-10).
nRoiC    = 480;         % ROI cells (12 Lambda = 12 unit cells)
iROI0C   = iSample0C + floor((nSampleC - nRoiC)/2);
iROI1C   = iROI0C + nRoiC - 1;

% Soft broadband source pulse in the LEFT background.
xcSrcC  = 4.0;
sxSrcC  = 0.8;
stSrcC  = 0.2*TC;
tcSrcC  = 1.0*TC;
wcSrcC  = 1.0*OmegaC;
JampC   = 1.0;
JsrcC   = @(x,t) JampC*exp(-((x - xcSrcC)/sxSrcC).^2) ...
                 .* exp(-((t - tcSrcC)/stSrcC).^2) ...
                 .* cos(wcSrcC*(t - tcSrcC));
sourceDC = @(x,t,n) dtC*JsrcC(x,t);

nPerC    = 38;                          % temporal periods in the run
nStepsC  = round(nPerC*TC/dtC);
recordEveryC = 2;                       % 50 samples per period
samplesPerPeriodC = round(TC/(dtC*recordEveryC));

bandCfgC = struct();
bandCfgC.x = xC;
bandCfgC.dt = dtC;
bandCfgC.nSteps = nStepsC;
bandCfgC.epsFun = epsFunC;
bandCfgC.muFun = @(x,t) ones(size(x));
bandCfgC.boundary = 'sponge';
bandCfgC.spongeCells = nSpongeC;
bandCfgC.spongeStrength = 0.08;
bandCfgC.recordEvery = recordEveryC;
bandCfgC.maxWaveSpeed = maxWaveSpeed;
bandCfgC.sourceD = sourceDC;           % the source is actually injected

bandFieldC = fdtd1d(bandCfgC);

% Self-contained time gate: same measured-arrival window as Case B.
roiEnergyC = sum(abs(bandFieldC.E(:, iROI0C:iROI1C)).^2, 2);
roiPeakC   = max(roiEnergyC);
arrIdxC    = find(roiEnergyC > 1e-3*roiPeakC, 1, 'first');
tStartTC   = ceil(bandFieldC.t(arrIdxC)/TC) + 1;
timeROIC   = [tStartTC*samplesPerPeriodC + 1, nPerC*samplesPerPeriodC];

% Pass both temporalPeriod and spatialPeriod -> fold BOTH k and omega into the
% first spacetime Brillouin zone: k in [-g/2, g/2], omega*T in [-pi, pi].
% The spatial ROI spans an integer number of unit cells; the time ROI an
% integer number of modulation periods.
tROIC1 = bandFieldC.t(timeROIC(1));
fftCfgC = struct('temporalPeriod', TC, 'spatialPeriod', LambdaC, ...
    'fieldComponent', 'E', ...
    'spatialROI', [iROI0C iROI1C], ...
    'timeROI', timeROIC, ...
    'spatialWindow', 'tukey', 'spatialWindowTukeyR', 0.5, ...
    'temporalWindow', 'hann', ...
    'noiseFloor', 1e-9, ...
    'measurementROI', [xC(iROI0C) xC(iROI1C)], ...
    'timeWindowStart', tROIC1, ...
    'modulationPhase', mod(tROIC1, TC)/TC);
bandsC = fdtd_fft_bands(bandFieldC.E, bandFieldC.x, bandFieldC.t, fftCfgC);

% PWE theory overlay (entry script only): folded quasifrequencies of the
% INFINITE spacetime crystal at the ROI's folded k columns. Nspace = Mtime = 1
% is exact for the single-harmonic standing wave (m,n = (+-1,+-1)).
kScanC  = bandsC.kNormalized(:).' * OmegaC;
cfgC    = struct('Lambda', LambdaC, 'T', TC, 'Nspace', 1, 'Mtime', 1, ...
                 'kScan', kScanC);
fourierC = pwe_fourier(epsFunC_smooth, [], cfgC);
pweC     = pwe_bands(fourierC, cfgC, 'omega');
pweWeight = pweC.m0Weight.';                     % 2S x nK
pweShow   = pweWeight >= 0.3;                    % well-resolved branches
pweReT    = real(pweC.omega)*TC;                 % Re(omega)*T, folded to [-pi,pi]
                                                 % (matches the omegaT rad axis)
kRep      = repmat(bandsC.kNormalized(:).', 2*pweC.S, 1);

figC = figure('Color','w','Position',[80 80 920 580]);
axC = axes(figC);
imagesc(axC, bandsC.kNormalized(:), bandsC.omegaT(:), bandsC.spectralDb);
axis(axC, 'xy');
colormap(axC, parula(256));
cbC = colorbar(axC); cbC.Label.String = 'folded spectral power [dB]';
hold(axC, 'on');
hPWE = scatter(axC, kRep(pweShow), pweReT(pweShow), 6, [0.1 0.6 0.1], ...
    'filled', 'DisplayName', 'PWE');
xlabel(axC, 'folded k  (first spatial BZ, k/\Omega)');
ylabel(axC, '\omega T  (first temporal BZ)');
title(axC, sprintf(['Finite spacetime crystal: double-folded E-FFT ' ...
    '\\epsilon=%.1f(1+%.1f cos(gx)cos(\\Omegat))'], eps0C, deltaC));
xlim(axC, [min(bandsC.kNormalized(:)) max(bandsC.kNormalized(:))]);
ylim(axC, [-pi pi]);
yticks(axC, [-pi 0 pi]); yticklabels(axC, {'-\pi','0','\pi'});
legend(axC, hPWE, 'Location','northeast','Box','off');
grid(axC, 'on'); box(axC, 'on');
drawnow;

fprintf('  sample slab x in [%.2f, %.2f], ROI x in [%.2f, %.2f] (%d cells)\n', ...
    xSampleC(1), xSampleC(2), xC(iROI0C), xC(iROI1C), nRoiC);
fprintf('  spectrum %dx%d (native d(omega*T)=%.4f, dk/Omega=%.4f)\n', ...
    size(bandsC.spectralDb,1), size(bandsC.spectralDb,2), ...
    bandsC.dftBinSpacingOmegaT, bandsC.spatialDftBinSpacingKOverK0);
fprintf('  effective (ENBW) resolutions: d(omega*T)=%.4f, dk/Omega=%.4f\n', ...
    bandsC.effectiveFreqResolutionOmegaT, bandsC.spatialEffectiveResolutionKOverK0);

% ---- Case C independent ridge vs PWE fundamental band (report-only) ----
% The ridge is extracted INDEPENDENTLY (max-power folded bin per ACTIVE
% column, noiseFloor-masked); the target is the PWE fundamental band of the
% INFINITE crystal, identified per column as the branch closest (wrapped) to
% the mean-medium light line omega*T = k*T/sqrt(eps0) (theory-based,
% ridge-independent). The folded-BZ sign degeneracy at the spatial band edge
% is resolved with min(|r-f|, |r+f|) as in Case B. Report-only: the
% quantitative acceptance for Case C lives in the section-7 regression battery
% (repair.md 7.313, run per 6.5), which adds the cross-length L=960/1280
% continuity leg this single-length entry run cannot.
kCnorm = bandsC.kNormalized(:).';
selCRidge = bandsC.activeKMask(:).';
ridgeC = nan(1, numel(kCnorm));
for jj = find(selCRidge)
    [~, imax] = max(bandsC.power(:, jj), [], 1);
    ridgeC(jj) = bandsC.omegaT(imax);
end
lightC = kScanC/sqrt(eps0C);                    % omega*T on the light line
fundC = nan(1, numel(kCnorm));
for jj = 1:numel(kCnorm)
    d = wrap_omegaT(pweReT(:, jj) - lightC(jj));
    [~, bj] = min(abs(d));
    fundC(jj) = pweReT(bj, jj);
end
errCR = nan(1, numel(kCnorm));
for jj = find(selCRidge)
    d1 = wrap_omegaT(ridgeC(jj) - fundC(jj));
    d2 = wrap_omegaT(ridgeC(jj) + fundC(jj));   % BZ-edge sign degeneracy
    errCR(jj) = min(abs(d1), abs(d2));
end
absvgC = 1/sqrt(eps0C);
linewidthCT = sqrt(bandsC.effectiveFreqResolutionOmegaT^2 + ...
    (2*pi*absvgC*bandsC.spatialEffectiveResolutionKOverK0)^2);
tolMedC = 2*pi/bandsC.nFold;
tolMaxC = max(2*tolMedC, linewidthCT/2);
fprintf('  independent ridge (max-power bin per active column) vs PWE\n');
fprintf('    fundamental band (light-line-closest branch): %d active cols,\n', ...
    nnz(selCRidge));
fprintf('    median |d(omega*T)| = %.4f (<= %.4f), max = %.4f (<= %.4f) [report only]\n', ...
    median(errCR(selCRidge)), tolMedC, max(errCR(selCRidge)), tolMaxC);

%% ======================= CONSOLE SUMMARY =======================
fprintf('============================================================\n');
fprintf('FDTD finite-sample E-FFT band benchmark vs TMM (binary time crystal):\n');
fprintf('  sample lengths (Lambda cells) : %s\n', mat2str(sampleLengthsB));
fprintf('  main benchmark (%d cells, ROI %d cells):\n', bMain.nSampleB, bMain.nRoiB);
fprintf('    median  |Delta(omega*T)| = %.4f   (<= %.4f pre-registered)\n', ...
    bMain.medianErrB, tolMedianB);
fprintf('    maximum |Delta(omega*T)| = %.4f   (<= %.4f pre-registered)\n', ...
    bMain.maxErrB, tolMaxB);
fprintf('    pass-band k columns compared  : %d\n', bMain.nPassB);
fprintf('  1/L series median: %s;  max: %s;  L->inf intercepts median %.4f, max %.4f\n', ...
    mat2str([benchB.medianErrB], 3), mat2str([benchB.maxErrB], 3), ...
    exB.medIntercept, exB.maxIntercept);
fprintf('  1/L common-k extrapolated residual (L->inf): median %.4f, max %.4f (1/L slope median %.4f)\n', ...
    exB.eInfMed, exB.eInfMax, exB.slopeMed);
fprintf('  pre-registered tolerance from the zero-padded folded grid spacing %.4f (NOT the old periodic 0.10/0.30):\n', ...
    zeroPadGridB);
fprintf('    median <= %.4f (= spacing/2), max <= %.4f (= one spacing)\n', ...
    tolMedianB, tolMaxB);
fprintf('  robustness: bg-buffer-doubled median %.4f (drift %.4f); ', ...
    variantBg.medianErrB, driftBg);
fprintf('sponge-halved median %.4f (drift %.4f); drift tol %.4f\n', ...
    variantSp.medianErrB, driftSp, driftTol);
fprintf('  measurement config: sample x in [%.3f, %.3f] (%.0f cells, free termination),\n', ...
    (caseBCfg.nSpongeB + caseBCfg.nBgLB)*caseBCfg.dxB, ...
    (caseBCfg.nSpongeB + caseBCfg.nBgLB + bMain.nSampleB)*caseBCfg.dxB, ...
    bMain.nSampleB);
fprintf('    ROI x in [%.3f, %.3f] (%.1f cells), ', ...
    bandsB.roiXStart, bandsB.roiXEnd, bMain.nRoiB);
fprintf('window [%.0f, %.0f]T (12 periods, Hann), ', ...
    bMain.tWinStartT, bMain.tWinStartT + 12);
fprintf('modulation eps %g<->%g from t=0 phase 0, ', ...
    epsHi, epsLo);
fprintf('source xc=%.1f sx=%.2f (no carrier) tc=%.1fT, background eps=%g mu=1\n', ...
    caseBCfg.xcSrcB, caseBCfg.sxSrcB, caseBCfg.tcSrcB/Tper, epsBg);
fprintf('    native sampling: kernel dx = %.4f Lambda, kernel dt = %.4f T (%d steps/period);\n', ...
    caseBCfg.dxB, caseBCfg.dtB/Tper, round(Tper/caseBCfg.dtB));
fprintf('    recorded every %d steps -> dtRecord = %.4f T (%d recorded samples/period)\n', ...
    caseBCfg.recordEveryB, caseBCfg.dtB*caseBCfg.recordEveryB/Tper, ...
    round(Tper/(caseBCfg.dtB*caseBCfg.recordEveryB)));
fprintf('    spatial window rect, temporal window Hann, active-k noiseFloor=%g\n', ...
    bMain.noiseFloorB);
fprintf('  no-sample reference: incident k/Omega support at -40 dB = %.2f (>= %.2f)\n', ...
    refB.supportKNorm, caseBCfg.maxKNormB);
fprintf('  no-sample reference: sponge residual backscatter = %.2e of incident peak\n', ...
    refB.reflectRatio);
fprintf('FDTD->FFT resolution, main Case B ROI (native = DFT bin spacing;\n');
fprintf('  effective = window ENBW x bin spacing, the real resolution limit):\n');
fprintf('  temporal bin spacing d(omega*T) = %.4f (= 2*pi/%d periods native)\n', ...
    bandsB.dftBinSpacingOmegaT, round(2*pi/bandsB.dftBinSpacingOmegaT));
fprintf('  temporal effective resolution (Hann ENBW) = %.4f\n', ...
    bandsB.effectiveFreqResolutionOmegaT);
fprintf('  spatial  bin spacing dk/Omega = %.4f; effective (ENBW) = %.4f\n', ...
    bandsB.spatialDftBinSpacingKOverK0, bandsB.spatialEffectiveResolutionKOverK0);
fprintf('============================================================\n');
fprintf(['CASE B TMM BENCHMARK PASSED (median |d(omega*T)| = %.4f <= %.4f, ', ...
    'max = %.4f <= %.4f pre-registered; incident support + sponge reflection ', ...
    'verified; 1/L extrapolated residual within tolerance).\n'], ...
    bMain.medianErrB, tolMedianB, bMain.maxErrB, tolMaxB);
fprintf('Cases A and C were generated without quantitative acceptance in 3.0.1.\n');
fprintf('The acceptance-adjacent classes NOT gated here live in the one-time\n');
fprintf('  regression battery (repair.md section 7, run per 6.5; /tmp, not shipped):\n');
fprintf('  two record durations, two ROI/spatial-window classes, the dx/dt grid\n');
fprintf('  convergence, the Case C cross-length ridge-vs-PWE quantitative comparison,\n');
fprintf('  and the fixed-central-ROI class -- a DOCUMENTED PHYSICAL EXPECTED-FAIL\n');
fprintf('  (a short central ROI cannot hold a steady-state window longer than the\n');
fprintf('  pulse transit, so every k column leaks); it is recorded there, not gated.\n');
fprintf('Base MATLAB only (no Toolbox); declared minimum R2020a (untested on this machine), verified R2026a.\n');
fprintf('============================================================\n');

%% ===================== OPTIONAL SAVE (off by default) =====================
% By default nothing is written and no directory is created. Flip this guard
% to true to export the figures to a local 'output' dir.
if doSave
    outDir = fullfile(pwd, 'output');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    save_fig(figA, fullfile(outDir, 'run_fdtd_caseA_field.png'));
    save_fig(figB, fullfile(outDir, 'run_fdtd_caseB_bands_vs_tmm.png'));
    save_fig(figC, fullfile(outDir, 'run_fdtd_caseC_doublefolded_bz.png'));
    fprintf('Saved figures to %s\n', outDir);
end

% -------------------------------------------------------------------------
function res = run_caseB_benchmark(nSampleB, c, epsScalar, muFun, epsBg, ...
    Tper, Omega, maxWaveSpeed)
%RUN_CASEB_BENCHMARK Finite binary time crystal: E-FFT band ridges inside the
%sample-interior ROI vs the exact 2x2 D/B TMM pass bands. Returns one result.
%
% Geometry (P0-01): sponge | bg | source | bg | SAMPLE | bg | sponge. The
% sample is the only modulated region. ROI = sample minus 10 cells per face
% (endpoint-free open window).
%
% Method A ridge: the ridge is the max-power bin over ALL omega*T per active
% k column (theory-independent); the error is min(|ridge-exact|, |ridge+exact|)
% over the TMM pass band (|trace/2| <= 1), allowing the +/-omega symmetry of
% the real-field spectrum.

dxB = c.dxB; dtB = c.dtB;
nSpongeB = c.nSpongeB; nBgLB = c.nBgLB;
nLeftB = nSpongeB + nBgLB;
NxB = 2*nLeftB + nSampleB;
xB = (0:NxB-1)*dxB;                      % E grid (no repeated endpoint)
iSample0B = nLeftB + 1;
iSample1B = nLeftB + nSampleB;
xSampleB = [xB(iSample0B) xB(iSample1B)];
iROI0B = iSample0B + 10;                 % sample minus 10 cells per face
iROI1B = iSample1B - 10;
nRoiB = iROI1B - iROI0B + 1;

inSampleB = @(X) double(X >= xSampleB(1) & X <= xSampleB(2));
epsFunB = @(X,t) epsBg + (epsScalar(t) - epsBg).*inSampleB(X);

% Soft broadband source (no carrier) in the left background (P0-01/P2-08):
% sourceD = dt*J is injected into D after each Ampere update.
JsrcB = @(x,t) exp(-((x - c.xcSrcB)/c.sxSrcB).^2) ...
               .* exp(-((t - c.tcSrcB)/c.stSrcB).^2);
sourceDB = @(x,t,n) dtB*JsrcB(x,t);

nStepsB = round(c.nPerB*Tper/dtB);
samplesPerPeriodB = round(Tper/(dtB*c.recordEveryB));

bandCfg = struct();
bandCfg.x = xB;
bandCfg.dt = dtB;
bandCfg.nSteps = nStepsB;
bandCfg.epsFun = epsFunB;
bandCfg.muFun = muFun;
bandCfg.boundary = 'sponge';
bandCfg.spongeCells = nSpongeB;
bandCfg.spongeStrength = c.spongeStrengthB;
bandCfg.recordEvery = c.recordEveryB;
bandCfg.temporalInterfaces = (Tper/2 : Tper/2 : nStepsB*dtB - Tper/2);
bandCfg.maxWaveSpeed = maxWaveSpeed;
bandCfg.sourceD = sourceDB;              % the source is actually injected

bandField = fdtd1d(bandCfg);

% Measured arrival + self-contained window. arrival = first ROI-energy
% crossing of 1e-3*max, window = [ceil(arrival)+10, ceil(arrival)+22]
% periods (12 periods, Hann) -> robust to group velocity. The threshold is
% referenced to the max over the FIRST HALF of the record: the momentum-gap
% parametric growth (Im*T ~ 0.12 for eps 1.3<->1) makes the full-record max
% grow with record length, so the first-half max is the incident-pulse
% plateau before that growth dominates and the arrival is record-length
% independent (verified for 44 and 60 period records).
roiEnergyB = sum(abs(bandField.E(:, iROI0B:iROI1B)).^2, 2);
roiHalfB = floor(numel(roiEnergyB)/2);
roiPeakB = max(roiEnergyB(1:roiHalfB));
arrIdxB = find(roiEnergyB(1:roiHalfB) > 1e-3*roiPeakB, 1, 'first');
tArriveT = bandField.t(arrIdxB)/Tper;
tWinStartT = ceil(tArriveT) + 10;
tWinEndT = tWinStartT + 12;
timeROIB = [tWinStartT*samplesPerPeriodB + 1, tWinEndT*samplesPerPeriodB];

% E-FFT inside the sample-interior measurement ROI. No spatialPeriod (the
% medium inside the ROI is spatially uniform), so k folding is off; only
% omega is folded. rect spatial window + Hann temporal window (validated).
tROIB1 = bandField.t(timeROIB(1));
fftCfgB = struct('temporalPeriod', Tper, ...
    'fieldComponent', 'E', ...
    'spatialROI', [iROI0B iROI1B], ...
    'timeROI', timeROIB, ...
    'spatialWindow', 'rect', ...
    'temporalWindow', 'hann', ...
    'noiseFloor', 1e-9, ...
    'measurementROI', [xB(iROI0B) xB(iROI1B)], ...
    'timeWindowStart', tROIB1, ...
    'modulationPhase', mod(tROIB1, Tper)/Tper);
bandsB = fdtd_fft_bands(bandField.E, bandField.x, bandField.t, fftCfgB);

% TMM pass band at the active, excited columns: signed |k/Omega| <= maxKNormB.
kN_B = bandsB.kNormalized(:).';
kSelB = (abs(kN_B) <= c.maxKNormB) & bandsB.activeKMask(:).';
kValuesB = kN_B(kSelB)*Omega;
tmmB = tmm_bands(kValuesB, [epsScalar(0) epsScalar(Tper/2)], 1, ...
    [Tper/2 Tper/2]);
halfTraceB = real(tmmB.halfTrace);
passbandB = abs(halfTraceB) <= 1 + 1e-9;
exactOmegaTB = acos(min(1, max(-1, halfTraceB)));    % in [0, pi]

% Method A ridge: max-power bin over ALL omega*T per active column.
ridgeA = nan(1, numel(kN_B));
for jj = find(kSelB)
    [~, imax] = max(bandsB.power(:, jj), [], 1);
    ridgeA(jj) = bandsB.omegaT(imax);
end
errsB = min(abs(ridgeA(kSelB) - exactOmegaTB), ...
            abs(ridgeA(kSelB) + exactOmegaTB));
pbErr = errsB(passbandB);
res.medianErrB = median(pbErr);
res.maxErrB = max(pbErr);
res.nPassB = nnz(passbandB);
res.noiseFloorB = fftCfgB.noiseFloor;
res.bandsB = bandsB;
res.kSelB = kSelB;
res.passbandB = passbandB;
res.exactOmegaTB = exactOmegaTB;
res.ridgeA = ridgeA;
res.nSampleB = nSampleB;
res.nRoiB = nRoiB;
res.tWinStartT = tWinStartT;
res.arrivalT = tArriveT;
end

% -------------------------------------------------------------------------
function res = run_caseB_reference(c, nSampleB, epsBg, muFun, Tper, maxWaveSpeed)
%RUN_CASEB_REFERENCE No-sample reference (repair.md line 309): identical Case
%B geometry but with eps = epsBg everywhere and a LONGER left background so
%the incident spectrum can be measured in a wide, quiet, sponge-free window.
%Quantifies
%   (1) supportKNorm : the incident k/Omega support at -40 dB (must cover the
%       benchmark statistics range |k/Omega| <= maxKNormB), and
%   (2) reflectRatio : the sponge residual backscatter in the left background
%       late in the record, normalized by the incident peak |E|.

dxB = c.dxB; dtB = c.dtB;
nSpongeB = c.nSpongeB;
nBgLB = 640;                            % longer left background (16 Lambda)
nLeftRef = nSpongeB + nBgLB;
NxRef = 2*nLeftRef + nSampleB;
xRef = (0:NxRef-1)*dxB;

epsFunRef = @(X,t) epsBg*ones(size(X));
JsrcRef = @(x,t) exp(-((x - c.xcSrcB)/c.sxSrcB).^2) ...
                 .* exp(-((t - c.tcSrcB)/c.stSrcB).^2);
sourceDRef = @(x,t,n) dtB*JsrcRef(x,t);

nStepsRef = round(c.nPerB*Tper/dtB);
refCfg = struct('x', xRef, 'dt', dtB, 'nSteps', nStepsRef, ...
    'epsFun', epsFunRef, 'muFun', muFun, 'boundary', 'sponge', ...
    'spongeCells', nSpongeB, 'spongeStrength', c.spongeStrengthB, ...
    'recordEvery', c.recordEveryB, 'maxWaveSpeed', maxWaveSpeed, ...
    'sourceD', sourceDRef);
refField = fdtd1d(refCfg);
Eref = refField.E;
tref = refField.t(:);

% ---- (1) incident k/Omega support at -40 dB ----
% Wide incident window: x in [5, 15] (400 cells = 10 Lambda, dk/Omega = 0.1),
% t in [5, 8] periods (the pulse crosses the window centre x = 10 at t = 7).
ix0 = round(5.0/dxB) + 1;
ix1 = round(15.0/dxB);                   % last inclusive cell
it0 = find(tref >= 5.0*Tper, 1, 'first');
it1 = find(tref <= 8.0*Tper, 1, 'last');
winE = Eref(it0:it1, ix0:ix1);
winE = winE - mean(winE, 1);             % per-column DC removal
S = fft2(winE);
Pow = abs(fftshift(S)).^2;
nXS = size(winE, 2);
kC = (-floor(nXS/2):ceil(nXS/2)-1);      % centered spatial bin index
% k/Omega = k*T/(2*pi) with k = 2*pi*kC/(N*dx) -> kC*T/(N*dx) (2*pi cancels).
kNormGrid = kC*Tper/(nXS*dxB);
peakPow = max(Pow, [], 'all');
supportMask = Pow >= 1e-4*peakPow;       % -40 dB contour
res.supportKNorm = max(abs(kNormGrid(any(supportMask, 1))));

% ---- (2) sponge residual backscatter at a left-bg probe, late in record ----
iProbeL = round(4.5/dxB) + 1;            % x = 4.5
iProbeH = round(6.5/dxB) + 1;            % x = 6.5
itIncident = tref <= 6.0*Tper;           % pulse transit through the probes
itLate = tref >= 12.0*Tper;              % after the pulse has cleared them
res.incidentPeakE = max(abs(Eref(itIncident, iProbeL:iProbeH)), [], 'all');
res.reflectRatio = max(abs(Eref(itLate, iProbeL:iProbeH)), [], 'all') ...
                   / res.incidentPeakE;
end

% -------------------------------------------------------------------------
function ex = run_caseB_1L_extrapolation(benchB, sampleLengthsB, maxKNormB, ...
    epsHi, epsLo, Tper, Omega)
%RUN_CASEB_1L_EXTRAPOLATION Cross-length ridge continuity + 1/L -> 0
%extrapolation of the finite-sample E-FFT ridge error (repair.md 7.311).
%
% Each sample length has its own ROI length, hence its own k grid, so the
% per-column Method-A ridge error is first interpolated onto a common signed-k
% grid. At every common-k point with a valid three-length fit, e(k) is fit by
% least squares in 1/L: e(k) = a(k) + b(k)/L. a(k) is the L -> inf residual
% (the large-sample ridge-vs-TMM agreement) and b(k) is the 1/L trend. The
% common-k grid is additionally masked to the exact TMM passband (same
% modulation across all samples): the ridge error is only defined there, and
% the linear k-interpolation must not bridge across a momentum gap. The
% aggregate median/max errors are also fit in 1/L for the figure.
%
% Returns
%   ex.commonK        : common signed-k grid (k/Omega)
%   ex.eInf, ex.eInfK : extrapolated L->inf residual and its k values
%   ex.eInfMed/eInfMax: median/max of the extrapolated residual
%   ex.slopeMed       : median |b(k)| (1/L trend strength)
%   ex.medIntercept/ex.medSlope, ex.maxIntercept/ex.maxSlope : aggregate
%       median/max 1/L fits; ex.medFit = [medSlope medIntercept] (polyfit
%       convention, for polyval).
nL = numel(sampleLengthsB);
commonK = linspace(0, maxKNormB, 401).';
eInterp = nan(nL, numel(commonK));
spanK   = false(nL, numel(commonK));
for li = 1:nL
    bb = benchB(li);
    ksel = bb.kSelB;                                  % logical, length nK
    errFull = nan(size(ksel));
    errFull(ksel) = min(abs(bb.ridgeA(ksel) - bb.exactOmegaTB), ...
                        abs(bb.ridgeA(ksel) + bb.exactOmegaTB));
    pbFull = false(size(ksel));
    pbFull(ksel) = bb.passbandB;
    valid = pbFull & isfinite(errFull) & isfinite(bb.bandsB.kNormalized(:).');
    kk = bb.bandsB.kNormalized(valid); ee = errFull(valid);
    if numel(kk) >= 3
        [kk, isort] = sort(kk); ee = ee(isort);
        eInterp(li,:) = interp1(kk, ee, commonK, 'linear');
        spanK(li,:) = (commonK >= min(kk)) & (commonK <= max(kk));
    end
end
% Exact TMM passband on the common grid: the ridge error is only defined in
% the passband, and interp1's 'linear' would otherwise bridge across a
% momentum gap. The modulation is identical across samples, so one TMM call.
halfTraceC = real(tmm_bands(commonK*Omega, [epsHi epsLo], 1, ...
    [Tper/2 Tper/2]).halfTrace);
passbandCommon = abs(halfTraceC) <= 1 + 1e-9;
good = all(spanK, 1) & all(isfinite(eInterp), 1) & passbandCommon(:).';
invL = (1./sampleLengthsB(:));
a1L = nan(1, numel(commonK)); b1L = nan(1, numel(commonK));
for ik = find(good)
    p = polyfit(invL, eInterp(:, ik), 1);
    a1L(ik) = p(2); b1L(ik) = p(1);                  % a + b/L
end
pMed = polyfit(invL, [benchB.medianErrB], 1);
pMax = polyfit(invL, [benchB.maxErrB], 1);

ex.commonK = commonK;
ex.eInf    = a1L(good);
ex.eInfK   = commonK(good);
ex.eInfMed = median(a1L(good));
ex.eInfMax = max(a1L(good));
ex.slopeMed = median(abs(b1L(good)));
ex.medIntercept = pMed(2); ex.medSlope = pMed(1); ex.medFit = pMed;
ex.maxIntercept = pMax(2); ex.maxSlope = pMax(1); ex.maxFit = pMax;
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

% -------------------------------------------------------------------------
function d = wrap_omegaT(a)
%WRAP_OMEGAT Wrap an omega*T value onto the period-2*pi circle [-pi, pi].
d = a - 2*pi*round(a/(2*pi));
end
