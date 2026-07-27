function example_ptc_pwe_vs_tmm()
%EXAMPLE_PTC_PWE_VS_TMM Cross-validate PWE against exact TMM for a binary PTC.
%
% The same quasifrequency bands of a binary photonic time crystal are
% calculated by two independent methods:
%   1. Temporal Fourier PWE (plane-wave expansion)
%   2. Exact 2x2 D/B monodromy (temporal transfer matrix)
%
% This cross-validation is essential for verifying convergence of the
% Fourier truncation, especially near band edges.

rootDir = stm_init();
epsA = 1.0;
epsB = 4.0;
muA = 1.0;
muB = 1.0;
dutyA = 0.5;
T = 1;
durations = [dutyA, 1-dutyA]*T;
Omega = 2*pi/T;

kNorm = linspace(0.02, 1.45, 181); % k*c0/Omega, c0=1
kValues = kNorm*Omega;
tmmBands = temporal_crystal_bands(kValues, [epsA epsB], ...
    [muA muB], durations);

% A binary temporal waveform is discontinuous, so it converges much more
% slowly than the sinusoidal modulation in the Park-Min model.
Mtime = 19;
epsCoeff = @(m,n) double(n == 0)*temporal_binary_eps_coeff( ...
    m, epsA, epsB, dutyA);
muCoeff = @(m,n) double(m == 0 && n == 0);
sys = stpwe_build_system(epsCoeff, muCoeff, 0, Mtime, 1, Omega);

pweMatched = complex(nan(2,numel(kValues)));
matchedWeight = nan(2,numel(kValues));
for ik = 1:numel(kValues)
    sol = stpwe_solve_omega(sys, kValues(ik));
    folded = stpwe_fold_frequency(sol.omega, Omega);
    valid = isfinite(folded) & abs(imag(folded)) < 0.35*Omega;
    ids = find(valid);
    used = false(size(ids));
    for branch = 1:2
        target = tmmBands.omegaF(branch,ik);
        cost = abs(folded(ids)-target)/Omega ...
            + 2e-3*(1-sol.m0Weight(ids));
        cost(used) = inf;
        [~, loc] = min(cost);
        chosen = ids(loc);
        used(loc) = true;
        pweMatched(branch,ik) = folded(chosen);
        matchedWeight(branch,ik) = sol.m0Weight(chosen);
    end
end

fig = figure('Color','w','Position',[100 100 1100 460]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact');

ax1 = nexttile(tl);
hold(ax1,'on');
plot(ax1, kNorm, real(tmmBands.omegaF.')/Omega, 'k-', ...
    'LineWidth',1.6);
kScatter = repmat(kNorm,2,1);
scatter(ax1, kScatter(:), real(pweMatched(:))/Omega, ...
    12, matchedWeight(:), 'filled');
xlabel(ax1, 'kc/\Omega');
ylabel(ax1, 'Re(\omega_F)/\Omega');
title(ax1, 'Principal Floquet bands');
grid(ax1,'on'); box(ax1,'on');
caxis(ax1,[0 1]);
cb = colorbar(ax1); cb.Label.String = 'PWE m=0 weight';

ax2 = nexttile(tl);
hold(ax2,'on');
plot(ax2, kNorm, imag(tmmBands.omegaF.')/Omega, 'k-', ...
    'LineWidth',1.6);
plot(ax2, kNorm, imag(pweMatched.')/Omega, 'o', ...
    'Color',[0.82 0.15 0.15], 'MarkerSize',3);
xlabel(ax2, 'kc/\Omega');
ylabel(ax2, 'Im(\omega_F)/\Omega');
title(ax2, 'Momentum-gap growth/decay');
grid(ax2,'on'); box(ax2,'on');
yline(ax2,0,'k:');

title(tl, 'Binary photonic time crystal: PWE (dots) vs TMM (lines)');
outputFile = fullfile(rootDir, 'output', ...
    'example_ptc_pwe_vs_tmm.png');
try
    exportgraphics(fig, outputFile, 'Resolution',220);
catch
    print(fig, outputFile, '-dpng', '-r220');
end
fprintf('Saved %s\n', outputFile);
end
