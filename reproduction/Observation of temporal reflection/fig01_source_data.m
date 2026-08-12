function result = fig01_source_data()
%FIG01_SOURCE_DATA Reproduce Fig. 1e/f from the official source-data files.
%   RESULT = FIG01_SOURCE_DATA() reads the measured/simulated phase and
%   oscilloscope/ADS waveforms, writes output/fig01_source_data.png and .mat,
%   and prints quantitative checks.  No project-core function or toolbox is
%   required.  The total phase length is 30*0.208 m, following Methods.

    scriptDir = fileparts(mfilename('fullpath'));
    addpath(scriptDir);
    dataDir = local_find_data_dir(scriptDir);
    outputDir = fullfile(scriptDir, 'output');
    if exist(outputDir, 'dir') ~= 7
        mkdir(outputDir);
    end

    phaseDir = fullfile(dataDir, '1e_raw');
    waveDir = fullfile(dataDir, '1f_raw');
    files = struct();
    files.offMeasured = fullfile(phaseDir, 'switch_off_Measured.txt');
    files.offSimulation = fullfile(phaseDir, 'switch_off_Sim.txt');
    files.onMeasured = fullfile(phaseDir, 'switch_on_Measured.txt');
    files.onSimulation = fullfile(phaseDir, 'switch_on_Sim.txt');
    files.scope = fullfile(waveDir, 'tek0000.csv');
    files.adsInput = fullfile(waveDir, 'V_in_sim.txt');
    files.adsOutput = fullfile(waveDir, 'V_out_sim.txt');

    names = fieldnames(files);
    for ii = 1:numel(names)
        if exist(files.(names{ii}), 'file') ~= 2
            error('otr:MissingSourceData', 'Missing Fig. 1 source file: %s', files.(names{ii}));
        end
    end

    tables = struct();
    tables.offMeasured = otr_read_numeric_table(files.offMeasured);
    tables.offSimulation = otr_read_numeric_table(files.offSimulation);
    tables.onMeasured = otr_read_numeric_table(files.onMeasured);
    tables.onSimulation = otr_read_numeric_table(files.onSimulation);

    cellCount = 30;
    cellLength = 0.208; % m, Methods section of the paper
    totalLength = cellCount * cellLength;
    dispersion = struct();
    dispersion.offMeasured = local_phase_to_beta(tables.offMeasured, totalLength);
    dispersion.offSimulation = local_phase_to_beta(tables.offSimulation, totalLength);
    dispersion.onMeasured = local_phase_to_beta(tables.onMeasured, totalLength);
    dispersion.onSimulation = local_phase_to_beta(tables.onSimulation, totalLength);

    scope = otr_read_scope_csv(files.scope);
    adsInput = otr_read_numeric_table(files.adsInput);
    adsOutput = otr_read_numeric_table(files.adsOutput);
    if size(scope.channels, 2) < 3
        error('otr:ScopeChannels', 'Fig. 1f CSV must contain CH1, CH2 and CH3.');
    end

    [~, measuredPeakIndex] = max(abs(scope.channels(:, 1)));
    [~, adsPeakIndex] = max(abs(adsInput.data(:, 2)));
    timeShift = adsInput.data(adsPeakIndex, 1) - scope.time(measuredPeakIndex);
    measuredTime = scope.time + timeShift;

    fitBandMHz = [10 70];
    velocities = struct();
    velocities.offMeasured = local_group_velocity(dispersion.offMeasured, fitBandMHz);
    velocities.offSimulation = local_group_velocity(dispersion.offSimulation, fitBandMHz);
    velocities.onMeasured = local_group_velocity(dispersion.onMeasured, fitBandMHz);
    velocities.onSimulation = local_group_velocity(dispersion.onSimulation, fitBandMHz);

    inputFrequencyMHz = 60;
    betaAtInput = interp1(dispersion.offMeasured.frequencyHz, ...
        dispersion.offMeasured.beta, inputFrequencyMHz * 1e6, 'linear');
    translatedFrequencyHz = local_interp_unique(dispersion.onMeasured.beta, ...
        dispersion.onMeasured.frequencyHz, betaAtInput);
    measuredTranslationRatio = translatedFrequencyHz / (inputFrequencyMHz * 1e6);

    blue = [0.10 0.40 0.80];
    red = [0.85 0.18 0.16];
    measuredGray = [0.15 0.15 0.15];
    green = [0.12 0.58 0.28];
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1120 760]);

    ax1 = subplot(2, 2, [1 3]);
    hold(ax1, 'on');
    local_plot_dispersion(ax1, dispersion.offSimulation, blue, '-', false);
    local_plot_dispersion(ax1, dispersion.onSimulation, red, '-', false);
    local_plot_dispersion(ax1, dispersion.offMeasured, blue, 'o', true);
    local_plot_dispersion(ax1, dispersion.onMeasured, red, 'o', true);
    plot(ax1, NaN, NaN, '-', 'Color', blue, 'LineWidth', 1.8, ...
        'DisplayName', 'switch OFF');
    plot(ax1, NaN, NaN, '-', 'Color', red, 'LineWidth', 1.8, ...
        'DisplayName', 'switch ON');
    plot(ax1, NaN, NaN, '-', 'Color', measuredGray, 'LineWidth', 1.6, ...
        'DisplayName', 'ADS simulation');
    plot(ax1, NaN, NaN, 'o', 'Color', measuredGray, 'MarkerFaceColor', 'w', ...
        'DisplayName', 'Measurement');
    xlabel(ax1, '\beta (rad m^{-1})');
    ylabel(ax1, 'f (MHz)');
    title(ax1, '(e) Dispersion extracted from phase / total length');
    xlim(ax1, [-11 11]);
    ylim(ax1, [-105 105]);
    grid(ax1, 'on');
    box(ax1, 'on');
    legend(ax1, 'Location', 'northwest');
    text(ax1, -10.3, -96, sprintf('L = 30 \\times 0.208 m = %.2f m', totalLength), ...
        'FontSize', 9, 'BackgroundColor', 'w');

    ax2 = subplot(2, 2, 2);
    hold(ax2, 'on');
    useMeasured = measuredTime >= 0 & measuredTime <= 300e-9;
    useAds = adsInput.data(:, 1) >= 0 & adsInput.data(:, 1) <= 300e-9;
    plot(ax2, measuredTime(useMeasured) * 1e9, scope.channels(useMeasured, 1), ...
        'Color', [0.10 0.10 0.10], 'LineWidth', 1.0, 'DisplayName', 'V_1 measured');
    plot(ax2, adsInput.data(useAds, 1) * 1e9, adsInput.data(useAds, 2), '--', ...
        'Color', red, 'LineWidth', 1.2, 'DisplayName', 'V_1 ADS');
    ylabel(ax2, 'V_1 (V)');
    title(ax2, '(f) Input port: measured and ADS');
    xlim(ax2, [0 300]);
    grid(ax2, 'on');
    box(ax2, 'on');
    legend(ax2, 'Location', 'northeast');

    % Put the switch monitor on a compact inset-like right axis.  The
    % measured TIME axis was shifted only to align the incident peak to ADS.
    yyaxis(ax2, 'right');
    logic = scope.channels(:, 3);
    logicLow = local_quantile(logic(useMeasured), 0.05);
    logicHigh = local_quantile(logic(useMeasured), 0.95);
    logicNorm = (logic - logicLow) ./ max(eps, logicHigh - logicLow);
    plot(ax2, measuredTime(useMeasured) * 1e9, logicNorm(useMeasured), ...
        'Color', green, 'LineWidth', 0.8, 'DisplayName', 'switch logic');
    ylabel(ax2, 'logic (normalized)');
    ylim(ax2, [-0.15 1.15]);
    yyaxis(ax2, 'left');

    ax3 = subplot(2, 2, 4);
    hold(ax3, 'on');
    useAdsOut = adsOutput.data(:, 1) >= 0 & adsOutput.data(:, 1) <= 300e-9;
    plot(ax3, measuredTime(useMeasured) * 1e9, scope.channels(useMeasured, 2), ...
        'Color', [0.10 0.10 0.10], 'LineWidth', 1.0, 'DisplayName', 'V_2 measured');
    plot(ax3, adsOutput.data(useAdsOut, 1) * 1e9, adsOutput.data(useAdsOut, 2), '--', ...
        'Color', red, 'LineWidth', 1.2, 'DisplayName', 'V_2 ADS');
    xlabel(ax3, 'aligned time (ns)');
    ylabel(ax3, 'V_2 (V)');
    title(ax3, 'Time-refracted output');
    xlim(ax3, [0 300]);
    grid(ax3, 'on');
    box(ax3, 'on');
    legend(ax3, 'Location', 'northeast');

    local_set_fonts(fig);
    pngFile = fullfile(outputDir, 'fig01_source_data.png');
    matFile = fullfile(outputDir, 'fig01_source_data.mat');
    print(fig, pngFile, '-dpng', '-r180');
    close(fig);

    result = struct();
    result.description = ['Official source-data reproduction of Fig. 1e/f. ' ...
        'beta=-phase/L; measured TIME shifted only to align its incident peak with ADS.'];
    result.dataDirectory = dataDir;
    result.files = files;
    result.totalLengthM = totalLength;
    result.dispersion = dispersion;
    result.scope = scope;
    result.adsInput = adsInput;
    result.adsOutput = adsOutput;
    result.measuredTimeShiftS = timeShift;
    result.fitBandMHz = fitBandMHz;
    result.groupVelocityMPerS = velocities;
    result.translationCheck.inputFrequencyMHz = inputFrequencyMHz;
    result.translationCheck.betaRadPerM = betaAtInput;
    result.translationCheck.outputFrequencyMHz = translatedFrequencyHz / 1e6;
    result.translationCheck.ratio = measuredTranslationRatio;
    result.outputFiles = {pngFile, matFile};
    save(matFile, 'result');

    fprintf('\n[Fig. 1 source data] beta = -phase(rad) / %.3f m\n', totalLength);
    fprintf('  Low-frequency group velocity (10--70 MHz), Mm/s:\n');
    fprintf('    OFF measured %.3f, OFF ADS %.3f\n', ...
        velocities.offMeasured / 1e6, velocities.offSimulation / 1e6);
    fprintf('    ON  measured %.3f, ON  ADS %.3f\n', ...
        velocities.onMeasured / 1e6, velocities.onSimulation / 1e6);
    fprintf('  Conserved-beta check: %.1f MHz OFF -> %.2f MHz ON (ratio %.3f).\n', ...
        inputFrequencyMHz, translatedFrequencyHz / 1e6, measuredTranslationRatio);
    fprintf('  Waveform alignment: measured TIME shifted by %+.3f ns to align input peaks.\n', ...
        timeShift * 1e9);
    fprintf('  Wrote %s and %s\n', pngFile, matFile);
end

function out = local_phase_to_beta(tbl, totalLength)
    data = tbl.data;
    valid = size(data, 2) >= 2 & isfinite(data(:, 1)) & isfinite(data(:, 2));
    frequency = data(valid, 1);
    phaseDeg = data(valid, 2);
    [frequency, order] = sort(frequency);
    phaseDeg = phaseDeg(order);
    out = struct('frequencyHz', frequency, 'phaseDeg', phaseDeg, ...
        'beta', -phaseDeg * pi / 180 / totalLength, 'source', tbl.filename);
end

function velocity = local_group_velocity(dispersion, bandMHz)
    use = dispersion.frequencyHz >= bandMHz(1) * 1e6 & ...
        dispersion.frequencyHz <= bandMHz(2) * 1e6 & dispersion.beta > 0;
    beta = dispersion.beta(use);
    omega = 2 * pi * dispersion.frequencyHz(use);
    velocity = (beta' * omega) / (beta' * beta);
end

function value = local_interp_unique(x, y, query)
    [x, order] = sort(x(:));
    y = y(order);
    [x, uniqueIndex] = unique(x, 'stable');
    y = y(uniqueIndex);
    value = interp1(x, y, query, 'linear');
end

function local_plot_dispersion(ax, d, color, style, measured)
    use = d.frequencyHz >= 0 & d.frequencyHz <= 100e6 & isfinite(d.beta);
    beta = d.beta(use);
    fMHz = d.frequencyHz(use) / 1e6;
    if measured
        stride = max(1, round(numel(beta) / 32));
        select = 1:stride:numel(beta);
        plot(ax, beta(select), fMHz(select), style, 'Color', color, ...
            'MarkerSize', 4.0, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
        plot(ax, -beta(select), -fMHz(select), style, 'Color', color, ...
            'MarkerSize', 4.0, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
    else
        plot(ax, beta, fMHz, style, 'Color', color, 'LineWidth', 1.7, ...
            'HandleVisibility', 'off');
        plot(ax, -beta, -fMHz, style, 'Color', color, 'LineWidth', 1.7, ...
            'HandleVisibility', 'off');
    end
end

function value = local_quantile(x, q)
    x = sort(x(isfinite(x)));
    if isempty(x)
        value = NaN;
        return;
    end
    position = 1 + (numel(x) - 1) * q;
    lo = floor(position);
    hi = ceil(position);
    value = x(lo) + (position - lo) * (x(hi) - x(lo));
end

function dataDir = local_find_data_dir(scriptDir)
    direct = {fullfile(scriptDir, 'data', 'Time_Interface_Source_Data'), ...
        fullfile(scriptDir, 'Time_Interface_Source_Data'), ...
        fullfile(pwd, 'data', 'Time_Interface_Source_Data'), ...
        fullfile(pwd, 'Time_Interface_Source_Data')};
    for ii = 1:numel(direct)
        if exist(direct{ii}, 'dir') == 7
            dataDir = direct{ii};
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
    error('otr:DataDirectoryMissing', ...
        'Could not locate data/Time_Interface_Source_Data below %s or %s.', scriptDir, pwd);
end

function local_set_fonts(fig)
    axesList = findall(fig, 'Type', 'axes');
    for ii = 1:numel(axesList)
        set(axesList(ii), 'FontName', 'Helvetica', 'FontSize', 9, 'LineWidth', 0.8);
    end
end
