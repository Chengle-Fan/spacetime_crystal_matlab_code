%RUN_FDTD  FDTD + FFT-band entry script (V3).
%
% V3 conventions: 1D, scalar, non-dispersive media; normalized units
% (c0 = eps0 = mu0 = 1, mu_r = 1). Spatial period Lambda, temporal period T,
% g = 2*pi/Lambda, Omega = 2*pi/T. The FDTD kernel advances D and B on a
% Yee grid (leapfrog); a soft/impulsive broadband source excites many k and
% omega, and a space-time FFT of the recorded D(x,t) history gives the
% Floquet band structure (real spectrum ridges only; Im(omega) inside gaps is
% not recovered by the FFT).
%
% Three cases share one material definition (a binary time crystal):
%
%   Case A - field distribution: a Gaussian wave packet (k ~ 0.6*Omega, in
%            the first momentum gap) on a sponge-bounded domain -> |E|(x,t).
%   Case B - FDTD-FFT band structure of the spatially-uniform binary time
%            crystal, cross-checked against the exact 2x2 D/B monodromy
%            (tmm_bands). This is the acceptance benchmark.
%   Case C - a general spacetime crystal eps(x,t) = eps0*(1 + delta*
%            cos(g*x)*cos(Omega*t)); the double-folded first-spacetime-BZ
%            spectrum (acceptance demo).
%
% Plan: 定义材料 -> 周期/截断参数 -> 计算 -> 作图
%
% Results stay in the workspace. Nothing is written to disk and no directory
% is created unless the doSave guard at the very end is enabled. The entry
% script is the ONLY place allowed to overlay the explicit TMM theory curve.

clear; clc; close all;

%% ===================== GLOBAL PARAMETERS =====================
% Normalized units: c0 = eps0 = mu0 = 1. Spatially-uniform media have mu_r=1.

% ---- shared material: binary photonic time crystal ----
%   eps(t) = epsHi for t in [0, T/2),  epsLo for t in [T/2, T).
epsHi = 4;             % high permittivity (first half-period)
epsLo = 1;             % low  permittivity (second half-period)
mu_r  = 1;             % relative permeability (V3 default)
Tper  = 1;             % temporal period T
Omega = 2*pi/Tper;     % temporal modulation frequency Omega = 2*pi/T = 2*pi

% scalar eps(t) -> material handles. epsFun(x,t) must return an array the
% same size as x (the FDTD kernel evaluates it with a scalar time t).
epsScalar = @(t) epsHi + (epsLo - epsHi)*double(mod(t, Tper) >= Tper/2);
epsFun = @(x,t) epsScalar(t) .* ones(size(x));
muFun  = @(x,t) mu_r .* ones(size(x));

% ---- Case A: field distribution (wave packet + sponge boundary) ----
Lambda = 1;            % length scale for the packet domain (eps depends on t only)
dxA    = Lambda/40;    % spatial step
dtA    = 0.5*dxA;      % time step (CFL = 0.5)
NxA    = 36*40;        % grid points -> domain ~36*Lambda (no repeated endpoint)
nStepsA   = round(3*Tper/dtA);   % cover ~3 temporal periods
k0A       = 0.6*Omega;           % packet center wavenumber (inside 1st momentum gap)
sigmaA    = 4*Lambda;            % packet Gaussian width
xcA       = 10*Lambda;           % initial packet center (away from the sponge)
spongeCellsA    = 120;           % sponge width in cells (3*Lambda)
spongeStrengthA = 0.08;          % sponge damping strength
recordEveryA    = 2;             % record every 2 steps

% ---- Case B: FDTD-FFT bands (broadband source + periodic) vs TMM ----
dxB = Lambda/40;       % spatial step (Lambda = 1 here is a spatial unit)
dtB = 0.4*dxB;         % time step (CFL = 0.4)
NxB = 2048;            % spatial cells (spatially uniform medium -> no k folding)
nPerB   = 24;          % number of temporal periods (>= 8; more = finer omega)
nStepsB = round(nPerB*Tper/dtB);
srcSigXB  = 0.10;      % broadband source spatial width (narrow -> broad k)
maxKNormB = 3.0;       % k/Omega upper limit for the TMM comparison (well excited)

% ---- Case C: general spacetime crystal, double-folded BZ ----
eps0C = 2;             % mean permittivity
deltaC = 0.3;          % modulation depth
LambdaC = 1;           % spatial period
TC      = 1;           % temporal period
gC     = 2*pi/LambdaC; % spatial reciprocal vector
OmegaC = 2*pi/TC;      % temporal frequency
nSpatPerC = 8;         % spatial periods in the domain
nTempPerC = 8;         % temporal periods in the run
dxC = LambdaC/40;      % spatial step
dtC = 0.4*dxC;         % time step (CFL = 0.4)
NxC     = nSpatPerC*40;              % grid points -> exactly 8*LambdaC
nStepsC = round(nTempPerC*TC/dtC);
srcSigXC = 0.25;       % broadband source spatial width

% ---- output guard (off by default; set true to export figures) ----
doSave = false;

%% ==================== CASE A: FIELD DISTRIBUTION ====================
fprintf('=== Case A: wave packet + sponge boundary ===\n');

xA     = (0:NxA-1)*dxA;                       % E grid (no repeated endpoint)
tEndA  = nStepsA*dtA;
% Abrupt permittivity switches must sit on integer time nodes (centred D/B
% correction): t = T/2, T, 3T/2, ... < nStepsA*dtA.
interfaceTimesA = (Tper/2 : Tper/2 : tEndA - Tper/2);

% One-way wave packet. In the locally-uniform medium eps(t=0)=epsHi the phase
% velocity is v = 1/sqrt(epsHi*mu) and the impedance H/E = sqrt(eps/mu). The
% H field is placed on the H grid (x + dx/2) at t = -dt/2 (the dt/2 shift).
vPhaseA = 1/sqrt(epsHi*mu_r);
omega0A = k0A*vPhaseA;
E0A = exp(-((xA - xcA)/sigmaA).^2) .* exp(1i*k0A*(xA - xcA));
xHA = xA(1:end-1) + dxA/2;                    % H grid: Nx-1 points (sponge)
Hhalf0A = sqrt(epsHi/mu_r) * exp(1i*(k0A*dxA/2 + omega0A*dtA/2)) ...
          .* exp(-((xHA - xcA)/sigmaA).^2) .* exp(1i*k0A*(xHA - xcA));

fieldCfg = struct();
fieldCfg.x = xA;
fieldCfg.dt = dtA;
fieldCfg.nSteps = nStepsA;
fieldCfg.epsFun = epsFun;
fieldCfg.muFun = muFun;
fieldCfg.E0 = E0A;
fieldCfg.Hhalf0 = Hhalf0A;
fieldCfg.boundary = 'sponge';
fieldCfg.spongeCells = spongeCellsA;
fieldCfg.spongeStrength = spongeStrengthA;
fieldCfg.recordEvery = recordEveryA;
fieldCfg.storeD = true;
fieldCfg.temporalInterfaces = interfaceTimesA;

fieldResult = fdtd1d(fieldCfg);

% Reconstruct |E| = |D|/eps(t); eps depends only on t, so this is exact away
% from the (measure-zero) switch instants.
D_A = fieldResult.D;                    % nRecords x NxA
x_A = fieldResult.x(:).';               % 1 x NxA
t_A = fieldResult.t(:);                 % nRecords x 1
E_A = D_A ./ epsScalar(t_A);            % broadcast the scalar eps(t) over x

% Line-out: |E| at the packet center vs time.
[~, probeIdA] = min(abs(x_A - xcA));
probeE_A = abs(E_A(:, probeIdA));

figA = figure('Color','w','Position',[60 60 1160 440]);
tlA = tiledlayout(figA, 1, 2, 'TileSpacing','compact', 'Padding','compact');

axA1 = nexttile(tlA);
imagesc(axA1, x_A/Lambda, t_A/Tper, abs(E_A));
axis(axA1, 'xy');
colormap(axA1, parula(256));
cbA1 = colorbar(axA1); cbA1.Label.String = '|E|';
hold(axA1, 'on');
% Mark the inner edges of the two sponge absorbing layers.
spongeEdgeLo = spongeCellsA*dxA/Lambda;
spongeEdgeHi = (NxA - spongeCellsA)*dxA/Lambda;
xline(axA1, spongeEdgeLo, 'r--', 'sponge', 'LineWidth', 1.2);
xline(axA1, spongeEdgeHi, 'r--', 'sponge', 'LineWidth', 1.2);
xlabel(axA1, 'x/\Lambda');
ylabel(axA1, 't/T');
title(axA1, '|E| wave packet (sponge edges dashed)');
grid(axA1, 'on'); box(axA1, 'on');

axA2 = nexttile(tlA);
plot(axA2, t_A/Tper, probeE_A, '-', 'Color',[0.00 0.45 0.75], 'LineWidth', 1.4);
xlabel(axA2, 't/T');
ylabel(axA2, sprintf('|E| at x/\\Lambda = %.2f', xcA/Lambda));
title(axA2, 'Probe field vs time');
grid(axA2, 'on'); box(axA2, 'on');

title(tlA, sprintf(['Binary time crystal \\epsilon(t)=[%g %g]: wave packet ' ...
    'at k=%.2f\\Omega'], epsHi, epsLo, k0A/Omega));
drawnow;

fprintf('  Nx=%d, dx=%.4f, dt=%.4f, CFL=%.3f, %d steps (~%.1f T)\n', ...
    NxA, dxA, dtA, dtA/dxA, nStepsA, nStepsA*dtA/Tper);
fprintf('  sponge: %d cells/side (%.1f*Lambda), strength %.2f\n', ...
    spongeCellsA, spongeCellsA*dxA/Lambda, spongeStrengthA);

%% ==================== CASE B: FDTD-FFT BANDS VS TMM ====================
fprintf('=== Case B: broadband FDTD-FFT bands vs TMM ===\n');

xB     = (0:NxB-1)*dxB;                 % periodic domain (no repeated endpoint)
tEndB  = nStepsB*dtB;
interfaceTimesB = (Tper/2 : Tper/2 : tEndB - Tper/2);

% Broadband localized source: a real Gaussian in x (narrow -> broad k) that is
% injected impulsively at t=0 (H = 0), so the temporal content is broadband in
% omega and both +/-k and +/-omega branches are excited. (Equivalently one may
% drive cfg.sourceD(x,t,n) over the first few steps as a soft source.)
E0B = exp(-(xB/srcSigXB).^2);
Hhalf0B = zeros(size(xB));              % periodic boundary: H grid has Nx cells

bandCfg = struct();
bandCfg.x = xB;
bandCfg.dt = dtB;
bandCfg.nSteps = nStepsB;
bandCfg.epsFun = epsFun;
bandCfg.muFun = muFun;
bandCfg.E0 = E0B;
bandCfg.Hhalf0 = Hhalf0B;
bandCfg.boundary = 'periodic';
bandCfg.recordEvery = 1;
bandCfg.storeD = true;
bandCfg.temporalInterfaces = interfaceTimesB;

bandField = fdtd1d(bandCfg);

% FFT bands. No spatialPeriod is passed: the medium is spatially uniform, so
% k folding is disabled (only omega is folded into the first temporal BZ).
fftCfgB = struct('temporalPeriod', Tper);
bandsB  = fdtd_fft_bands(bandField.D, bandField.x, bandField.t, fftCfgB);

% Explicit TMM theory overlay (entry script only). kNormalized = k/Omega.
kN_B    = bandsB.kNormalized(:).';          % 1 x nK (k/Omega)
kSelB   = kN_B <= maxKNormB;                % columns the source excites reliably
kValuesB = kN_B(kSelB) * Omega;             % physical k = (k/Omega)*Omega
tmmB     = tmm_bands(kValuesB, [epsHi epsLo], mu_r, [Tper/2 Tper/2]);
halfTraceB   = real(tmmB.halfTrace);
passbandB    = abs(halfTraceB) <= 1 + 1e-9;
exactOmegaTB = acos(min(1, max(-1, halfTraceB)));   % in [0, pi]

ridgePosB = bandsB.ridgePositive(:).';
ridgeNegB = bandsB.ridgeNegative(:).';
ridgePosB = ridgePosB(kSelB);
ridgeNegB = ridgeNegB(kSelB);

% Benchmark errors over the pass band: positive ridge ~ +exactOmegaT, negative
% ridge ~ -exactOmegaT (same audit as the fig2 reference).
posErrB = abs(ridgePosB - exactOmegaTB);
negErrB = abs(ridgeNegB + exactOmegaTB);
errsB   = [posErrB(passbandB), negErrB(passbandB)];
medianErrB = median(errsB);
maxErrB    = max(errsB);
nPassB     = nnz(passbandB);

% True FDTD->FFT resolution (zero-padding refines the plot grid only).
dkB      = 2*pi/(NxB*dxB);       % spatial k bin spacing
dOmegaTB = 2*pi/nPerB;           % native temporal resolution in omega*T

fprintf('  median |Delta(omega*T)| = %.4f over %d pass-band k columns\n', ...
    medianErrB, nPassB);
fprintf('  max    |Delta(omega*T)| = %.4f\n', maxErrB);
assert(medianErrB <= 0.10, ...
    'median FDTD-FFT vs TMM error %.4f exceeds 0.10', medianErrB);
assert(maxErrB <= 0.30, ...
    'maximum FDTD-FFT vs TMM error %.4f exceeds 0.30', maxErrB);
fprintf('  PASS: FDTD-FFT ridges agree with the TMM pass bands.\n');

figB = figure('Color','w','Position',[70 70 920 580]);
axB = axes(figB);
imagesc(axB, kN_B, bandsB.omegaT(:), bandsB.spectralDb);
axis(axB, 'xy');
colormap(axB, parula(256));
cbB = colorbar(axB); cbB.Label.String = 'folded spectral power [dB]';
hold(axB, 'on');

% TMM pass-band curves (solid) over the pass-band segments only.
kSelBNum = kN_B(kSelB);
d = diff([false passbandB false]);
pbStart = find(d == 1); pbStop = find(d == -1) - 1;
for s = 1:numel(pbStart)
    ids = pbStart(s):pbStop(s);
    plot(axB, kSelBNum(ids),  exactOmegaTB(ids), 'k-', 'LineWidth', 1.4, ...
        'HandleVisibility','off');
    plot(axB, kSelBNum(ids), -exactOmegaTB(ids), 'k-', 'LineWidth', 1.4, ...
        'HandleVisibility','off');
end
% FDTD-FFT ridge markers (thinned over the pass band).
mkIds = find(passbandB);
mkIds = mkIds(1:2:end);
hRidge = scatter(axB, kSelBNum(mkIds), ridgePosB(mkIds), 14, 'w', 'filled');
scatter(axB, kSelBNum(mkIds), ridgeNegB(mkIds), 14, 'w', 'filled', ...
    'HandleVisibility','off');
hTheory = plot(axB, nan, nan, 'k-', 'LineWidth', 1.4, 'DisplayName','TMM');

xlabel(axB, 'k/\Omega   (\Omega = 2\pi/T)');
ylabel(axB, '\omega T');
title(axB, sprintf(['Binary time crystal: FDTD-FFT ridges vs TMM ' ...
    '(median |\\Delta\\omega T| = %.3f)'], medianErrB));
xlim(axB, [0 maxKNormB]); ylim(axB, [-pi pi]);
yticks(axB, [-pi 0 pi]); yticklabels(axB, {'-\pi','0','\pi'});
legend(axB, [hRidge hTheory], {'FDTD-FFT ridge','TMM pass band'}, ...
    'Location','northeast','Box','off');
grid(axB, 'on'); box(axB, 'on');
drawnow;

%% ==================== CASE C: DOUBLE-FOLDED BZ ====================
fprintf('=== Case C: general spacetime crystal, double-folded BZ ===\n');

% General spacetime crystal: eps(x,t) = eps0*(1 + delta*cos(g*x)*cos(Omega*t)).
epsFunC = @(x,t) eps0C*(1 + deltaC*cos(gC*x).*cos(OmegaC*t));
muFunC  = @(x,t) ones(size(x));

xC = (0:NxC-1)*dxC;            % exactly 8*LambdaC (periodic, no repeated endpoint)
E0C = exp(-(xC/srcSigXC).^2);  % broadband localized source (broad k)
Hhalf0C = zeros(size(xC));

bandCfgC = struct();
bandCfgC.x = xC;
bandCfgC.dt = dtC;
bandCfgC.nSteps = nStepsC;
bandCfgC.epsFun = epsFunC;
bandCfgC.muFun = muFunC;
bandCfgC.E0 = E0C;
bandCfgC.Hhalf0 = Hhalf0C;
bandCfgC.boundary = 'periodic';
bandCfgC.recordEvery = 1;
bandCfgC.storeD = true;

bandFieldC = fdtd1d(bandCfgC);

% Pass both temporalPeriod and spatialPeriod -> fold BOTH k and omega into the
% first spacetime Brillouin zone: k in [-g/2, g/2], omega*T in [-pi, pi].
fftCfgC = struct('temporalPeriod', TC, 'spatialPeriod', LambdaC);
bandsC  = fdtd_fft_bands(bandFieldC.D, bandFieldC.x, bandFieldC.t, fftCfgC);

figC = figure('Color','w','Position',[80 80 920 580]);
axC = axes(figC);
imagesc(axC, bandsC.kNormalized(:), bandsC.omegaT(:), bandsC.spectralDb);
axis(axC, 'xy');
colormap(axC, parula(256));
cbC = colorbar(axC); cbC.Label.String = 'folded spectral power [dB]';
xlabel(axC, 'folded k  (first spatial BZ)');
ylabel(axC, '\omega T  (first temporal BZ)');
title(axC, sprintf(['Double-folded first spacetime BZ: ' ...
    '\\epsilon=%.1f(1+%.1f cos(gx)cos(\\Omegat))'], eps0C, deltaC));
xlim(axC, [min(bandsC.kNormalized(:)) max(bandsC.kNormalized(:))]);
ylim(axC, [-pi pi]);
yticks(axC, [-pi 0 pi]); yticklabels(axC, {'-\pi','0','\pi'});
grid(axC, 'on'); box(axC, 'on');
drawnow;

fprintf('  eps(x,t) = %.1f*(1 + %.2f*cos(g*x)*cos(Omega*t)), g=%.4f, Omega=%.4f\n', ...
    eps0C, deltaC, gC, OmegaC);
fprintf('  domain: %d spatial x %d temporal periods, Nx=%d, %d steps\n', ...
    nSpatPerC, nTempPerC, NxC, nStepsC);

%% ======================= CONSOLE SUMMARY =======================
fprintf('============================================================\n');
fprintf('FDTD-FFT band benchmark vs TMM (binary time crystal):\n');
fprintf('  median  |Delta(omega*T)| = %.4f   (<= 0.10 required)\n', medianErrB);
fprintf('  maximum |Delta(omega*T)| = %.4f   (<= 0.30 required)\n', maxErrB);
fprintf('  pass-band k columns compared  : %d\n', nPassB);
fprintf('FDTD->FFT resolution (native, not zero-padded):\n');
fprintf('  spatial  dk = %.6f   (dk/Omega = %.6f)\n', dkB, dkB/Omega);
fprintf('  temporal d(omega*T) = %.6f   (= 2*pi/%d periods)\n', dOmegaTB, nPerB);
fprintf('  spatial k columns (Case B)    : %d\n', numel(kN_B));
fprintf('============================================================\n');
fprintf('ALL CHECKS PASSED.\n');

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
function save_fig(fig, filePath)
%SAVE_FIG Export a figure to PNG with a fallback if exportgraphics fails.
try
    exportgraphics(fig, filePath, 'Resolution', 220);
catch
    print(fig, filePath, '-dpng', '-r220');
end
end
