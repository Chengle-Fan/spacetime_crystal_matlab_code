function example_floquet_bands_and_fields(quality)
%EXAMPLE_FLOQUET_BANDS_AND_FIELDS Floquet band structure and eigenmode field patterns.
%
%   example_floquet_bands_and_fields('quick')
%   example_floquet_bands_and_fields('paper')
%
% Demonstrates the full ST-PWE workflow on a partially modulated unit cell:
%   1. Build the system matrix from analytic Fourier coefficients
%   2. Compute the Floquet band structure (fixed-k eigenvalue problem)
%   3. Select six representative eigenmodes across the Brillouin zone
%   4. Reconstruct and plot the E(x,t) field for each selected mode
%
% The color scale uses m=0 Floquet-sector participation as a transparent
% proxy for physical mode weight.

if nargin < 1
    quality = 'quick';
end
rootDir = stm_init();
p = stm_preset_modulated_slab();

switch lower(quality)
    case 'paper'
        Nspace = 20;
        Nk = 181;
    case 'quick'
        Nspace = 10;
        Nk = 101;
    otherwise
        error('quality must be ''quick'' or ''paper''.');
end
Mtime = 1;

epsCoeff = @(m,n) stm_fourier_modulated_slab(m,n,p);
sys = stpwe_build_system(epsCoeff, [], Nspace, Mtime, p.g, p.Omega);

kBar = linspace(-0.5, 0.5, Nk);
kValues = p.g*kBar;
fMax = 0.82;
weightThreshold = 0.035;
imagTolerance = 3e-3;

allK = cell(Nk,1);
allF = cell(Nk,1);
allWeight = cell(Nk,1);

fprintf('Floquet bands: %d k points, matrix size %d x %d.\n', ...
    Nk, 2*sys.S, 2*sys.S);
for ik = 1:Nk
    sol = stpwe_solve_omega(sys, kValues(ik));
    fBar = sol.omega/(p.g*p.c0);
    keep = isfinite(fBar) ...
        & real(fBar) >= 0 & real(fBar) <= fMax ...
        & abs(imag(fBar)) <= imagTolerance ...
        & sol.m0Weight >= weightThreshold;
    allK{ik} = repmat(kBar(ik), sum(keep), 1);
    allF{ik} = real(fBar(keep));
    allWeight{ik} = sol.m0Weight(keep);
end
allK = vertcat(allK{:});
allF = vertcat(allF{:});
allWeight = vertcat(allWeight{:});

epsCoeff0 = @(n) stm_fourier_modulated_slab(0,n,p);
fStatic = stpwe_static_bands(kValues, epsCoeff0, Nspace, ...
    p.g, p.c0, 8);

labels = {'c','d','e','f','g','h'};
kTargetBar = [-0.175, 0.175, 0.10, 0.29, 0.10, 0.29];
fTargetBar = [0.10, 0.10, 0.45, 0.43, 0.35, 0.32];
modes = cell(1,numel(labels));
for q = 1:numel(labels)
    modes{q} = stpwe_select_mode(sys, p.g*kTargetBar(q), ...
        p.g*p.c0*fTargetBar(q), [0, 0.9*p.g*p.c0]);
    fprintf('(%s): kBar=% .4f, fBar=% .5f%+.2ei, m0 weight=%.3f\n', ...
        labels{q}, modes{q}.k/p.g, ...
        real(modes{q}.omega/(p.g*p.c0)), ...
        imag(modes{q}.omega/(p.g*p.c0)), modes{q}.m0Weight);
end

fig = figure('Color','w','Position',[60 60 1450 850]);
tl = tiledlayout(fig, 3, 4, 'TileSpacing','compact', ...
    'Padding','compact');

axA = nexttile(tl, 1, [1 2]);
xPlot = linspace(0, 3*p.Lambda, 601);
tauPlot = linspace(-1, 1, 161);
[XX, Tau] = meshgrid(xPlot, tauPlot);
epsPlot = stm_permittivity_modulated_slab(XX, Tau*p.T, p);
surf(axA, XX/p.Lambda, Tau, epsPlot, 'EdgeColor','none');
view(axA, -38, 28);
xlabel(axA, 'Position x/\Lambda');
ylabel(axA, 'Time t/T');
zlabel(axA, '\epsilon_r');
title(axA, '(a) Space-time permittivity');
axis(axA, 'tight');
colormap(axA, hot(256));

axB = nexttile(tl, 5, [2 2]);
hold(axB, 'on');
for ib = 1:size(fStatic,1)
    plot(axB, kBar, fStatic(ib,:), 'Color',[0.25 0.25 0.25], ...
        'LineWidth',0.8);
    for shift = [-2 -1 1 2]
        plot(axB, kBar, fStatic(ib,:) + shift*p.OmegaBar, '--', ...
            'Color',[0.65 0.65 0.65], 'LineWidth',0.45);
    end
end
scatter(axB, allK, allF, 8, allWeight, 'filled');
for q = 1:numel(labels)
    kb = modes{q}.k/p.g;
    fb = real(modes{q}.omega/(p.g*p.c0));
    plot(axB, kb, fb, 'ko', 'MarkerFaceColor','w', 'MarkerSize',4);
    text(axB, kb + 0.012, fb + 0.012, ['(' labels{q} ')'], ...
        'FontSize',9);
end
xlim(axB, [-0.5 0.5]);
ylim(axB, [0 fMax]);
xlabel(axB, 'k\Lambda/(2\pi)');
ylabel(axB, '\omega\Lambda/(2\pi c)');
title(axB, '(b) Floquet band structure');
grid(axB, 'on');
box(axB, 'on');
cmapBand = [linspace(1,0,256).', zeros(256,1), ones(256,1)];
colormap(axB, cmapBand);
caxis(axB, [0 1]);
cb = colorbar(axB, 'northoutside');
cb.Label.String = 'm=0 sector participation';

xField = linspace(0, 3*p.Lambda, 451);
tauField = linspace(0, 3, 361);
fieldTiles = [3 4 7 8 11 12];
for q = 1:numel(labels)
    ax = nexttile(tl, fieldTiles(q));
    field = stpwe_reconstruct_field(modes{q}, xField, ...
        tauField*p.T, false);
    imagesc(ax, xField/p.Lambda, tauField, field);
    set(ax, 'YDir','normal');
    caxis(ax, [-1 1]);
    colormap(ax, stm_redblue(256));
    hold(ax, 'on');
    for cellNo = 0:2
        xline(ax, cellNo + p.xModStart/p.Lambda, 'k-', ...
            'LineWidth',0.45);
    end
    title(ax, ['(' labels{q} ')']);
    xlabel(ax, 'x/\Lambda');
    if mod(q,2) == 1
        ylabel(ax, 't/T');
    else
        set(ax, 'YTickLabel',[]);
    end
end

title(tl, 'Floquet bands and eigenmode fields (modulated slab)');
outputFile = fullfile(rootDir, 'output', ...
    ['example_floquet_bands_and_fields_' lower(quality) '.png']);
save_example_figure(fig, outputFile);
fprintf('Saved %s\n', outputFile);
end

function save_example_figure(fig, outputFile)
try
    exportgraphics(fig, outputFile, 'Resolution',220);
catch
    print(fig, outputFile, '-dpng', '-r220');
end
end
