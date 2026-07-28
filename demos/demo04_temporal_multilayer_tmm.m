function demo04_temporal_multilayer_tmm()
%DEMO04_TEMPORAL_MULTILAYER_TMM Reproduce the four APL design classes.
%
% The first panel uses Table I of Ramaccia et al. The remaining panels show
% transparent, periodic transparent, and amplifying equal-travel-distance
% temporal multilayers.

rootDir = startup_stm();
fNorm = linspace(0.5, 1.5, 401);
omega0 = 2*pi*fNorm; % T0=1

cases(1).title = 'Example 1: arbitrary stack';
cases(1).na = 1;
cases(1).nb = 3;
cases(1).n = [2.5 4.4 2.0 1.6];
cases(1).dt = [0.5 1.0 3.5 0.25];

cases(2).title = 'Example 2: transparent at f_0';
cases(2).na = 1;
cases(2).nb = 2;
cases(2).n = [2 7 1.5 5];
cases(2).dt = cases(2).n/(2*cases(2).na); % delta=pi

cases(3).title = 'Example 3: periodic transparent stack';
cases(3).na = 1;
cases(3).nb = 1;
cases(3).n = [6 1.5 6 1.5 6];
cases(3).dt = cases(3).n/(2*cases(3).na); % delta=pi

cases(4).title = 'Example 4: temporal amplification';
cases(4).na = 1;
cases(4).nb = 1;
cases(4).n = [6 1.5 6 1.5 6];
cases(4).dt = cases(4).n/(4*cases(4).na); % delta=pi/2

fig = figure('Color','w','Position',[80 80 1120 760]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact');
for q = 1:4
    epsInitial = cases(q).na^2;
    epsFinal = cases(q).nb^2;
    epsSlabs = cases(q).n.^2;
    spectrum = temporal_tmm_spectrum(omega0, epsInitial, 1, ...
        epsSlabs, ones(size(epsSlabs)), cases(q).dt, epsFinal, 1);

    ax = nexttile(tl);
    plot(ax, fNorm, abs(spectrum.forward), 'b-', 'LineWidth',1.5);
    hold(ax,'on');
    plot(ax, fNorm, abs(spectrum.backward), 'r--', 'LineWidth',1.5);
    xline(ax,1,'k:');
    xlabel(ax,'f/f_0');
    ylabel(ax,'Electric-field amplitude');
    title(ax,cases(q).title);
    legend(ax,'|E_b^+/E_a^+|','|E_b^-/E_a^+|', ...
        'Location','best');
    grid(ax,'on'); box(ax,'on');
end
title(tl, 'Temporal multilayer transfer functions');

outputFile = fullfile(rootDir, 'output', ...
    'demo04_temporal_multilayer_tmm.png');
try
    exportgraphics(fig, outputFile, 'Resolution',220);
catch
    print(fig, outputFile, '-dpng', '-r220');
end
fprintf('Saved %s\n', outputFile);
end
