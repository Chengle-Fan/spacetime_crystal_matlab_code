function fig1a_ptc_schematic()
%FIG1A_PTC_SCHEMATIC  Reproduce Fig. 1(a) of Lustig et al. (arXiv:1803.08731).
%
%   A schematic of the binary photonic time-crystal (PTC):
%   The permittivity epsilon(t) alternates between eps1 = 3 and eps2 = 1
%   with equal segment durations. The time origin t=0 is placed at the
%   midpoint of segment 1, preserving time-reversal symmetry.
%
%   Reference:
%     E. Lustig, Y. Sharabi, and M. Segev,
%     "Topology of photonic time-crystals," arXiv:1803.08731v1 (2018).

% --- Add parent toolbox to path ---
rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(rootDir, 'startup_stm.m'));

% --- Parameters (matching the paper) ---
eps1 = 3;
eps2 = 1;
T    = 2*pi;          % modulation period
t1   = 0.5*T;         % duration of segment 1
t2   = 0.5*T;         % duration of segment 2

% --- Build epsilon(t) over multiple periods ---
nPeriodsPlot = 4;     % how many periods to show
dtFine = T/200;       % fine time sampling
tPlot = linspace(-nPeriodsPlot*T/2, nPeriodsPlot*T/2, 2*nPeriodsPlot*200+1);
epsPlot = zeros(size(tPlot));

for i = 1:length(tPlot)
    tVal = tPlot(i);
    % Shift to [0, T) frame
    tMod = mod(tVal + T/2, T);  % t=0 center of seg1 → tMod goes 0 to T
    if tMod < t1
        epsPlot(i) = eps1;
    else
        epsPlot(i) = eps2;
    end
end

% --- Create figure ---
fig = figure('Color', 'w', 'Position', [100 100 800 400]);

% Plot epsilon(t)
fill([tPlot(1) tPlot, tPlot(end)], [0 epsPlot, 0], ...
    [0.7 0.85 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;
stairs(tPlot, epsPlot, 'b-', 'LineWidth', 2);

% Mark the time-reversal symmetry point
xline(0, 'r--', 'LineWidth', 1.5);
text(0.05, eps1+0.15, 't=0 (TR center)', 'Color', 'r', ...
    'FontSize', 11, 'FontWeight', 'bold');

% Annotate
ylim([0, eps1+0.8]);
xlabel('Time  t', 'FontSize', 13);
ylabel('Permittivity  \epsilon(t)', 'FontSize', 13);
title('Fig. 1(a): Binary photonic time-crystal – \epsilon(t)', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Annotate segment values
yMax = eps1 + 0.8;
plot([-T/4, -T/4], [0, yMax], 'k:', 'LineWidth', 1);
text(-T/4-0.1, eps1+0.15, sprintf('\\epsilon_1=%d', eps1), ...
    'FontSize', 12, 'HorizontalAlignment', 'right');
plot([T/4, T/4], [0, yMax], 'k:', 'LineWidth', 1);
text(T/4+0.1, eps2+0.15, sprintf('\\epsilon_2=%d', eps2), ...
    'FontSize', 12, 'HorizontalAlignment', 'left');

% Annotate period
plot([-T/2, T/2], [-0.25, -0.25], 'k-', 'LineWidth', 1.5);
plot([-T/2, -T/2], [-0.35, -0.15], 'k-', 'LineWidth', 1);
plot([T/2, T/2], [-0.35, -0.15], 'k-', 'LineWidth', 1);
text(0, -0.45, sprintf('T = %.4g  (\\Omega = %.4g)', T, 2*pi/T), ...
    'FontSize', 11, 'HorizontalAlignment', 'center');

% Axis formatting
set(gca, 'FontSize', 12);
box on;

% --- Save ---
outputDir = fullfile(fileparts(mfilename('fullpath')), 'output');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
outputFile = fullfile(outputDir, 'fig1a_ptc_schematic.png');
try
    exportgraphics(fig, outputFile, 'Resolution', 200);
catch
    print(fig, outputFile, '-dpng', '-r200');
end
fprintf('Saved: %s\n', outputFile);
end
