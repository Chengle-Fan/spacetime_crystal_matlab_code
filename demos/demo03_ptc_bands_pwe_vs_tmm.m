function demo03_ptc_bands_pwe_vs_tmm()
%DEMO03_PTC_BANDS_PWE_VS_TMM Cross-check a binary photonic time crystal.
%
% The same quasifrequency bands are calculated in two independent ways:
% temporal Fourier PWE and exact 2x2 D/B monodromy (TMM).

rootDir = startup_stm();
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
% slowly than the sinusoidal modulation in Fig. 2.
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
hTmm1 = plot(ax1, kNorm, real(tmmBands.omegaF.')/Omega, 'k-', ...
    'LineWidth',1.6);
kScatter = repmat(kNorm,2,1);
hPwe1 = scatter(ax1, kScatter(:), real(pweMatched(:))/Omega, ...
    12, matchedWeight(:), 'filled');
xlabel(ax1, 'kc/\Omega');
ylabel(ax1, 'Re(\omega_F)/\Omega');
title(ax1, 'Principal Floquet bands');
grid(ax1,'on'); box(ax1,'on');
clim(ax1,[0 1]);
cb = colorbar(ax1); cb.Label.String = 'PWE m=0 weight';
legend(ax1, [hTmm1(1), hPwe1], {'TMM (exact)','PWE (matched)'}, ...
    'Location','best','Box','off');

ax2 = nexttile(tl);
hold(ax2,'on');
hTmm2 = plot(ax2, kNorm, imag(tmmBands.omegaF.')/Omega, 'k-', ...
    'LineWidth',1.6);
hPwe2 = plot(ax2, kNorm, imag(pweMatched.')/Omega, 'o', ...
    'Color',[0.82 0.15 0.15], 'MarkerSize',3);
xlabel(ax2, 'kc/\Omega');
ylabel(ax2, 'Im(\omega_F)/\Omega');
title(ax2, 'Momentum-gap growth/decay');
grid(ax2,'on'); box(ax2,'on');
yline(ax2,0,'k:');
legend(ax2, [hTmm2(1), hPwe2(1)], {'TMM (exact)','PWE (matched)'}, ...
    'Location','best','Box','off');

title(tl, 'Binary photonic time crystal: PWE dots vs TMM lines');
outputFile = fullfile(rootDir, 'output', ...
    'demo03_ptc_pwe_vs_tmm.png');
try
    exportgraphics(fig, outputFile, 'Resolution',220);
catch
    print(fig, outputFile, '-dpng', '-r220');
end
fprintf('Saved %s\n', outputFile);

%% ---- Figure 2: how the binary PTC is modulated ----
nPer = 6;   % show two modulation periods
fig2 = figure('Color','w','Position',[130 620 1000 460]);
tl2 = tiledlayout(fig2, 1, 1, 'TileSpacing','compact');

% Left: binary epsilon(t) square wave over two periods.
ax3 = nexttile(tl2);
hold(ax3,'on');
xlim(ax3, [0 nPer]); ylim(ax3, [0.4 4.8]);
for p = 0:nPer-1   % lightly shade the epsA (duty) half of each period
    patch(ax3, [p, p+dutyA, p+dutyA, p], ...
        [0.4 0.4 4.8 4.8], [0.90 0.95 1.00], 'EdgeColor','none');
end
tPer = [0, dutyA*T, dutyA*T, T];       % corners of one period
epsPer = [epsA, epsA, epsB, epsB];
tFull = []; eFull = [];
for p = 0:nPer-1
    tFull = [tFull, tPer + p*T];
    eFull = [eFull, epsPer];
end
plot(ax3, tFull/T, eFull, 'b-', 'LineWidth', 2);
yline(ax3, epsA, ':', 'Color',[0 0.45 0.75]);
yline(ax3, epsB, ':', 'Color',[0.85 0.33 0.1]);
text(ax3, 0.25, epsA+0.35, '\epsilon_A = 1.0', ...
    'HorizontalAlignment','center','Color',[0 0.45 0.75]);
text(ax3, 0.75, epsB-0.35, '\epsilon_B = 4.0', ...
    'HorizontalAlignment','center','Color',[0.85 0.33 0.1]);
text(ax3, 1.0, 4.55, sprintf('one period T = %.2f', T), ...
    'HorizontalAlignment','center','FontSize',9);
xlabel(ax3, 't/T');
ylabel(ax3, '\epsilon(t)');
title(ax3, 'Temporal modulation \epsilon(t)');
grid(ax3,'on'); box(ax3,'on');

% % Right: |eps_m| of the binary waveform (only odd m for 50% duty, ~1/|m| decay).
% ax4 = nexttile(tl2);
% mList2 = (-25:25).';
% cAll = arrayfun(@(mm) temporal_binary_eps_coeff(mm, epsA, epsB, dutyA), ...
%     mList2);
% nz = abs(cAll) > 0;   % drop the exactly-zero even harmonics
% stem(ax4, mList2(nz), abs(cAll(nz)), 'b', 'filled', 'MarkerSize', 4);
% hold(ax4,'on');
% xline(ax4, Mtime, 'r--', 'LineWidth', 1);
% xline(ax4, -Mtime, 'r--', 'LineWidth', 1);
% text(ax4, Mtime+0.5, 0.95*max(abs(cAll)), 'Mtime = 19', ...
%     'HorizontalAlignment','left','Color','r','FontSize',8);
% xlabel(ax4, 'temporal harmonic m');
% ylabel(ax4, '|\epsilon_m|');
% title(ax4, 'Fourier spectrum: only odd m, \propto 1/|m|');
% grid(ax4,'on'); box(ax4,'on');

% title(tl2, sprintf(['How the binary PTC is modulated: \epsilon(t) switches ' ...
%     '%.1f <-> %.1f with duty %.0f%%'], epsA, epsB, 100*dutyA));
% outputFile2 = fullfile(rootDir, 'output', ...
%     'demo03_modulation_waveform.png');
% try
%     exportgraphics(fig2, outputFile2, 'Resolution',220);
% catch
%     print(fig2, outputFile2, '-dpng', '-r220');
% end
% fprintf('Saved %s\n', outputFile2);
% end
