function example_complex_gaps()
%EXAMPLE_COMPLEX_GAPS Compare momentum gaps and frequency gaps.
%
% This example demonstrates the two fundamental eigenvalue problems of ST-PWE:
%   - Fixed k -> complex omega reveals unstable momentum gaps (Eq. 6)
%   - Fixed omega -> complex k reveals evanescent frequency gaps (Eq. 7)
%
% See also STPWE_SOLVE_OMEGA, STPWE_SOLVE_K

rootDir = stm_init();
p = stm_preset_modulated_slab();
Nspace = 10;
Mtime = 1;
sys = stpwe_build_system(@(m,n) stm_fourier_modulated_slab(m,n,p), ...
    [], Nspace, Mtime, p.g, p.Omega);

kBarSweep = linspace(0.12, 0.23, 111);
kStore = cell(size(kBarSweep));
fStore = cell(size(kBarSweep));
fiStore = cell(size(kBarSweep));
for ik = 1:numel(kBarSweep)
    sol = stpwe_solve_omega(sys, p.g*kBarSweep(ik));
    f = sol.omega/(p.g*p.c0);
    keep = isfinite(f) & real(f) > 0.075 & real(f) < 0.125 ...
        & abs(imag(f)) < 0.03 & sol.m0Weight > 0.01;
    kStore{ik} = repmat(kBarSweep(ik), sum(keep), 1);
    fStore{ik} = real(f(keep));
    fiStore{ik} = imag(f(keep));
end

fBarSweep = linspace(0.385, 0.425, 121);
f2Store = cell(size(fBarSweep));
krStore = cell(size(fBarSweep));
kiStore = cell(size(fBarSweep));
for iw = 1:numel(fBarSweep)
    sol = stpwe_solve_k(sys, p.g*p.c0*fBarSweep(iw));
    kNormalized = sol.k/p.g;
    keep = isfinite(kNormalized) ...
        & real(kNormalized) > -0.5 & real(kNormalized) < 0.5 ...
        & abs(imag(kNormalized)) < 0.12;
    f2Store{iw} = repmat(fBarSweep(iw), sum(keep), 1);
    krStore{iw} = real(kNormalized(keep));
    kiStore{iw} = imag(kNormalized(keep));
end

fig = figure('Color','w','Position',[100 100 1150 720]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact');

nexttile(tl);
scatter(vertcat(kStore{:}), vertcat(fStore{:}), 8, ...
    abs(vertcat(fiStore{:})), 'filled');
xlabel('k\Lambda/(2\pi)');
ylabel('Re(\omega)\Lambda/(2\pi c)');
title('Fixed k: momentum gap');
grid on; box on;
cb = colorbar; cb.Label.String = '|Im(\omega)|\Lambda/(2\pi c)';

nexttile(tl);
scatter(vertcat(kStore{:}), vertcat(fiStore{:}), 8, ...
    vertcat(fStore{:}), 'filled');
xlabel('k\Lambda/(2\pi)');
ylabel('Im(\omega)\Lambda/(2\pi c)');
title('Parametric growth/decay pair');
grid on; box on;
yline(0, 'k:');

nexttile(tl);
scatter(vertcat(krStore{:}), vertcat(f2Store{:}), 8, ...
    abs(vertcat(kiStore{:})), 'filled');
xlabel('Re(k)\Lambda/(2\pi)');
ylabel('\omega\Lambda/(2\pi c)');
title('Fixed \omega: frequency gap');
grid on; box on;
cb = colorbar; cb.Label.String = '|Im(k)|\Lambda/(2\pi)';

nexttile(tl);
scatter(vertcat(f2Store{:}), vertcat(kiStore{:}), 8, ...
    vertcat(krStore{:}), 'filled');
xlabel('\omega\Lambda/(2\pi c)');
ylabel('Im(k)\Lambda/(2\pi)');
title('Evanescent Bloch solutions');
grid on; box on;
xline(0.4, 'k:');
yline(0, 'k:');

title(tl, 'Complex ST-PWE spectra: momentum and frequency gaps');
outputFile = fullfile(rootDir, 'output', ...
    'example_complex_gaps.png');
try
    exportgraphics(fig, outputFile, 'Resolution',220);
catch
    print(fig, outputFile, '-dpng', '-r220');
end
fprintf('Saved %s\n', outputFile);
end
