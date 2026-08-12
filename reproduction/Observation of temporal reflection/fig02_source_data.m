function result = fig02_source_data()
%FIG02_SOURCE_DATA Reproduce Fig. 2c--f from official source data.
%   RESULT = FIG02_SOURCE_DATA() plots measured red/blue-shift wave packets
%   in time and frequency and compares the frequency-sweep measurements to
%   the theory column supplied by the authors.  The gates used for each
%   component are explicit in RESULT.gates; processing is limited to the
%   documented baseline/taper operation in OTR_GATE_SIGNAL.
%
%   Outputs: output/fig02_source_data.png and .mat (MATLAB R2019b+).

    scriptDir = fileparts(mfilename('fullpath'));
    addpath(scriptDir);
    dataDir = local_find_data_dir(scriptDir);
    outputDir = fullfile(scriptDir, 'output');
    if exist(outputDir, 'dir') ~= 7
        mkdir(outputDir);
    end

    files.redshiftScope = fullfile(dataDir, '2c_raw', 'tek0002.csv');
    files.blueshiftScope = fullfile(dataDir, '2d_raw', 'tek0000.csv');
    files.redshiftSweep = fullfile(dataDir, '2e', 'freq_conv.csv');
    files.blueshiftSweep = fullfile(dataDir, '2f', 'freq_conv_reverse.csv');
    names = fieldnames(files);
    for ii = 1:numel(names)
        if exist(files.(names{ii}), 'file') ~= 2
            error('otr:MissingSourceData', 'Missing Fig. 2 source file: %s', files.(names{ii}));
        end
    end

    redScope = otr_read_scope_csv(files.redshiftScope);
    blueScope = otr_read_scope_csv(files.blueshiftScope);
    redSweepTable = otr_read_numeric_table(files.redshiftSweep);
    blueSweepTable = otr_read_numeric_table(files.blueshiftSweep);
    if size(redScope.channels, 2) < 2 || size(blueScope.channels, 2) < 2
        error('otr:ScopeChannels', 'Fig. 2 scope files must have CH1 and CH2.');
    end
    if size(redSweepTable.data, 2) < 4 || size(blueSweepTable.data, 2) < 4
        error('otr:SweepColumns', 'Fig. 2 sweep files must have four numeric columns.');
    end

    % Gates are expressed on the raw oscilloscope TIME coordinate.  CH1 is
    % the input-port trace (incident followed by TR); CH2 is the output-port
    % trace (time refraction).  They were selected to isolate the components
    % visible in the source waveforms and are saved verbatim for audit.
    gates = struct();
    gates.redshift.incident = [-30 40] * 1e-9;
    gates.redshift.reflected = [40 150] * 1e-9;
    gates.redshift.transmitted = [60 160] * 1e-9;
    gates.blueshift.incident = [-20 45] * 1e-9;
    gates.blueshift.reflected = [40 160] * 1e-9;
    gates.blueshift.transmitted = [35 130] * 1e-9;
    taperFraction = 0.10;

    red = local_analyse_scope(redScope, gates.redshift, taperFraction);
    blue = local_analyse_scope(blueScope, gates.blueshift, taperFraction);
    redSweep = redSweepTable.data(:, 1:4);
    blueSweep = blueSweepTable.data(:, 1:4);
    redErrors = [local_rmse(redSweep(:, 2), redSweep(:, 4)), ...
        local_rmse(redSweep(:, 3), redSweep(:, 4))];
    blueErrors = [local_rmse(blueSweep(:, 2), blueSweep(:, 4)), ...
        local_rmse(blueSweep(:, 3), blueSweep(:, 4))];

    cIncident = [0.12 0.42 0.80];
    cReflected = [0.56 0.20 0.67];
    cTransmitted = [0.95 0.65 0.08];
    cTheory = [0.20 0.20 0.20];
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [60 60 1280 780]);

    ax1 = subplot(2, 3, 1);
    local_plot_time(ax1, red, cIncident, cReflected, cTransmitted);
    title(ax1, '(c) Impedance decrease: time');
    ax2 = subplot(2, 3, 2);
    local_plot_spectra(ax2, red, cIncident, cReflected, cTransmitted);
    title(ax2, sprintf('Redshift: %.1f -> %.1f / %.1f MHz', ...
        red.peakMHz.incident, red.peakMHz.reflected, red.peakMHz.transmitted));

    ax3 = subplot(2, 3, 4);
    local_plot_time(ax3, blue, cIncident, cReflected, cTransmitted);
    title(ax3, '(d) Impedance increase: time');
    ax4 = subplot(2, 3, 5);
    local_plot_spectra(ax4, blue, cIncident, cReflected, cTransmitted);
    title(ax4, sprintf('Blueshift: %.1f -> %.1f / %.1f MHz', ...
        blue.peakMHz.incident, blue.peakMHz.reflected, blue.peakMHz.transmitted));

    ax5 = subplot(2, 3, 3);
    hold(ax5, 'on');
    plot(ax5, redSweep(:, 1), redSweep(:, 4), '--', 'Color', cTheory, ...
        'LineWidth', 1.5, 'DisplayName', 'theory CSV');
    plot(ax5, redSweep(:, 1), redSweep(:, 2), 'o-.', 'Color', cReflected, ...
        'MarkerFaceColor', 'w', 'LineWidth', 1.1, 'DisplayName', 'TR measured');
    plot(ax5, redSweep(:, 1), redSweep(:, 3), 's-', 'Color', cTransmitted, ...
        'MarkerFaceColor', 'w', 'LineWidth', 1.1, 'DisplayName', 'refracted measured');
    xlabel(ax5, 'input f_1 (MHz)');
    ylabel(ax5, 'output f_2 (MHz)');
    title(ax5, '(e) Broadband redshift sweep');
    grid(ax5, 'on');
    box(ax5, 'on');
    legend(ax5, 'Location', 'northwest');

    ax6 = subplot(2, 3, 6);
    hold(ax6, 'on');
    plot(ax6, blueSweep(:, 1), blueSweep(:, 4), '--', 'Color', cTheory, ...
        'LineWidth', 1.5, 'DisplayName', 'theory CSV');
    plot(ax6, blueSweep(:, 1), blueSweep(:, 2), 'o-.', 'Color', cReflected, ...
        'MarkerFaceColor', 'w', 'LineWidth', 1.1, 'DisplayName', 'TR measured');
    plot(ax6, blueSweep(:, 1), blueSweep(:, 3), 's-', 'Color', cTransmitted, ...
        'MarkerFaceColor', 'w', 'LineWidth', 1.1, 'DisplayName', 'refracted measured');
    xlabel(ax6, 'input f_2 (MHz)');
    ylabel(ax6, 'output f_1 (MHz)');
    title(ax6, '(f) Broadband blueshift sweep');
    grid(ax6, 'on');
    box(ax6, 'on');
    legend(ax6, 'Location', 'northwest');

    local_set_fonts(fig);
    pngFile = fullfile(outputDir, 'fig02_source_data.png');
    matFile = fullfile(outputDir, 'fig02_source_data.mat');
    print(fig, pngFile, '-dpng', '-r180');
    close(fig);

    result = struct();
    result.description = ['Official source-data reproduction of Fig. 2c--f. ' ...
        'FFT inputs use only the explicit gates below, linear endpoint baseline removal, ' ...
        'and a 10% cosine taper.'];
    result.dataDirectory = dataDir;
    result.files = files;
    result.gates = gates;
    result.taperFraction = taperFraction;
    result.redshift = red;
    result.blueshift = blue;
    result.redshiftSweep.headers = redSweepTable.headers;
    result.redshiftSweep.data = redSweep;
    result.redshiftSweep.rmseMHz = redErrors;
    result.blueshiftSweep.headers = blueSweepTable.headers;
    result.blueshiftSweep.data = blueSweep;
    result.blueshiftSweep.rmseMHz = blueErrors;
    result.outputFiles = {pngFile, matFile};
    save(matFile, 'result');

    fprintf('\n[Fig. 2 source data] Explicit time-gated FFT (no toolbox/filter)\n');
    fprintf('  Redshift peaks [incident, TR, refracted] = [%.2f, %.2f, %.2f] MHz\n', ...
        red.peakMHz.incident, red.peakMHz.reflected, red.peakMHz.transmitted);
    fprintf('  Redshift ratios [TR, refracted]/incident = [%.3f, %.3f]\n', ...
        red.peakMHz.reflected / red.peakMHz.incident, ...
        red.peakMHz.transmitted / red.peakMHz.incident);
    fprintf('  Blueshift peaks [incident, TR, refracted] = [%.2f, %.2f, %.2f] MHz\n', ...
        blue.peakMHz.incident, blue.peakMHz.reflected, blue.peakMHz.transmitted);
    fprintf('  Blueshift ratios [TR, refracted]/incident = [%.3f, %.3f]\n', ...
        blue.peakMHz.reflected / blue.peakMHz.incident, ...
        blue.peakMHz.transmitted / blue.peakMHz.incident);
    fprintf('  Sweep RMSE vs theory CSV [TR, refracted]: red [%.2f, %.2f] MHz, blue [%.2f, %.2f] MHz\n', ...
        redErrors(1), redErrors(2), blueErrors(1), blueErrors(2));
    fprintf('  Wrote %s and %s\n', pngFile, matFile);
end

function analysed = local_analyse_scope(scope, gateLimits, taperFraction)
    incident = otr_gate_signal(scope.time, scope.channels(:, 1), gateLimits.incident, taperFraction);
    reflected = otr_gate_signal(scope.time, scope.channels(:, 1), gateLimits.reflected, taperFraction);
    transmitted = otr_gate_signal(scope.time, scope.channels(:, 2), gateLimits.transmitted, taperFraction);
    [incident.frequencyHz, incident.spectrum, incident.peakHz] = local_spectrum(incident);
    [reflected.frequencyHz, reflected.spectrum, reflected.peakHz] = local_spectrum(reflected);
    [transmitted.frequencyHz, transmitted.spectrum, transmitted.peakHz] = local_spectrum(transmitted);

    analysed = struct();
    analysed.scope = scope;
    analysed.incident = incident;
    analysed.reflected = reflected;
    analysed.transmitted = transmitted;
    analysed.peakMHz = struct('incident', incident.peakHz / 1e6, ...
        'reflected', reflected.peakHz / 1e6, 'transmitted', transmitted.peakHz / 1e6);
end

function [frequency, amplitude, peakFrequency] = local_spectrum(gate)
    x = gate.signal(:);
    n = numel(x);
    nfft = 2 ^ nextpow2(max(n, 16384));
    spectrum = fft(x, nfft);
    count = floor(nfft / 2) + 1;
    frequency = (0:count - 1).' / (nfft * gate.dt);
    amplitude = abs(spectrum(1:count)) * gate.dt;
    search = frequency >= 5e6 & frequency <= 100e6;
    indices = find(search);
    [~, localIndex] = max(amplitude(search));
    peakIndex = indices(localIndex);
    peakFrequency = local_parabolic_peak(frequency, amplitude, peakIndex);
end

function peak = local_parabolic_peak(frequency, amplitude, index)
    peak = frequency(index);
    if index <= 1 || index >= numel(amplitude)
        return;
    end
    y = log(max(realmin, amplitude(index - 1:index + 1)));
    denominator = y(1) - 2 * y(2) + y(3);
    if abs(denominator) > eps
        offset = 0.5 * (y(1) - y(3)) / denominator;
        offset = max(-1, min(1, offset));
        peak = frequency(index) + offset * (frequency(index + 1) - frequency(index));
    end
end

function local_plot_time(ax, analysed, cIncident, cReflected, cTransmitted)
    hold(ax, 'on');
    plot(ax, analysed.incident.time * 1e9, analysed.incident.detrended, ...
        'Color', cIncident, 'LineWidth', 1.1, 'DisplayName', 'incident (CH1)');
    plot(ax, analysed.reflected.time * 1e9, analysed.reflected.detrended, ...
        'Color', cReflected, 'LineWidth', 1.1, 'DisplayName', 'TR (CH1)');
    plot(ax, analysed.transmitted.time * 1e9, analysed.transmitted.detrended, ...
        'Color', cTransmitted, 'LineWidth', 1.1, 'DisplayName', 'refracted (CH2)');
    xlabel(ax, 'raw scope time (ns)');
    ylabel(ax, 'baseline-corrected V (V)');
    grid(ax, 'on');
    box(ax, 'on');
    legend(ax, 'Location', 'best');
end

function local_plot_spectra(ax, analysed, cIncident, cReflected, cTransmitted)
    hold(ax, 'on');
    local_spectrum_line(ax, analysed.incident, cIncident, 'incident');
    local_spectrum_line(ax, analysed.reflected, cReflected, 'TR');
    local_spectrum_line(ax, analysed.transmitted, cTransmitted, 'refracted');
    xlabel(ax, 'f (MHz)');
    ylabel(ax, 'normalized |FFT|');
    xlim(ax, [0 100]);
    ylim(ax, [0 1.05]);
    grid(ax, 'on');
    box(ax, 'on');
    legend(ax, 'Location', 'northeast');
end

function local_spectrum_line(ax, gate, color, name)
    use = gate.frequencyHz <= 100e6;
    amplitude = gate.spectrum;
    amplitude = amplitude / max(eps, max(amplitude(use)));
    plot(ax, gate.frequencyHz(use) / 1e6, amplitude(use), 'Color', color, ...
        'LineWidth', 1.3, 'DisplayName', name);
end

function value = local_rmse(measured, theory)
    use = isfinite(measured) & isfinite(theory);
    value = sqrt(mean((measured(use) - theory(use)).^2));
end

function dataDir = local_find_data_dir(scriptDir)
    candidates = {fullfile(scriptDir, 'data', 'Time_Interface_Source_Data'), ...
        fullfile(scriptDir, 'Time_Interface_Source_Data'), ...
        fullfile(pwd, 'data', 'Time_Interface_Source_Data'), ...
        fullfile(pwd, 'Time_Interface_Source_Data')};
    for ii = 1:numel(candidates)
        if exist(candidates{ii}, 'dir') == 7
            dataDir = candidates{ii};
            return;
        end
    end
    roots = unique({scriptDir, pwd});
    for rr = 1:numel(roots)
        paths = regexp(genpath(roots{rr}), pathsep, 'split');
        for ii = 1:numel(paths)
            [~, leaf] = fileparts(paths{ii});
            if strcmp(leaf, 'Time_Interface_Source_Data')
                dataDir = paths{ii};
                return;
            end
        end
    end
    error('otr:DataDirectoryMissing', 'Could not locate data/Time_Interface_Source_Data.');
end

function local_set_fonts(fig)
    axesList = findall(fig, 'Type', 'axes');
    for ii = 1:numel(axesList)
        set(axesList(ii), 'FontName', 'Helvetica', 'FontSize', 9, 'LineWidth', 0.8);
    end
end
