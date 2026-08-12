function result = fig03_source_data()
%FIG03_SOURCE_DATA Plot the official raw Fig. 3c/d source waveforms.
%   RESULT = FIG03_SOURCE_DATA() loads tau_15/25/35.csv and displays CH1,
%   CH2 and the switch monitor on their native oscilloscope time coordinate.
%   The supplied archive does not state unique gates or the separate control
%   acquisition needed to reconstruct the processed reflection spectra in
%   Fig. 3d.  This routine therefore labels its output RAW DATA and performs
%   only robust offset removal for visualization; it does not claim a gated
%   Fig. 3d reconstruction.
%
%   Outputs: output/fig03_source_data.png and .mat (MATLAB R2019b+).

    scriptDir = fileparts(mfilename('fullpath'));
    addpath(scriptDir);
    dataDir = local_find_data_dir(scriptDir);
    outputDir = fullfile(scriptDir, 'output');
    if exist(outputDir, 'dir') ~= 7
        mkdir(outputDir);
    end

    tauNs = [15 25 35];
    fileNames = {'tau_15.csv', 'tau_25.csv', 'tau_35.csv'};
    files = cell(size(fileNames));
    scopes = cell(size(fileNames));
    traces = cell(size(fileNames));
    edgeSummary = repmat(struct('riseTimeNs', NaN, 'fallTimeNs', NaN, ...
        'highDurationNs', NaN, 'ch1PeakToPeakV', NaN, 'ch2PeakToPeakV', NaN), ...
        size(fileNames));

    for ii = 1:numel(fileNames)
        files{ii} = fullfile(dataDir, '3cd_raw', fileNames{ii});
        if exist(files{ii}, 'file') ~= 2
            error('otr:MissingSourceData', 'Missing Fig. 3 source file: %s', files{ii});
        end
        scopes{ii} = otr_read_scope_csv(files{ii});
        if size(scopes{ii}.channels, 2) < 3
            error('otr:ScopeChannels', '%s must contain CH1, CH2 and CH3.', files{ii});
        end
        traces{ii} = local_raw_trace(scopes{ii}, [-40 230] * 1e-9);
        edgeSummary(ii) = local_edge_summary(traces{ii});
    end

    c1 = [0.10 0.10 0.10];
    c2 = [0.12 0.43 0.80];
    cLogic = [0.12 0.58 0.28];
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [70 70 1080 850]);
    for ii = 1:numel(tauNs)
        ax = subplot(3, 1, ii);
        hold(ax, 'on');
        trace = traces{ii};
        plot(ax, trace.time * 1e9, trace.ch1, 'Color', c1, 'LineWidth', 1.0, ...
            'DisplayName', 'CH1 raw - robust offset');
        plot(ax, trace.time * 1e9, trace.ch2, 'Color', c2, 'LineWidth', 1.0, ...
            'DisplayName', 'CH2 raw - robust offset');
        ylabel(ax, 'voltage (V)');
        xlim(ax, [-40 230]);
        grid(ax, 'on');
        box(ax, 'on');
        title(ax, sprintf('RAW DATA: nominal \\tau = %d ns (scope logic high %.1f ns)', ...
            tauNs(ii), edgeSummary(ii).highDurationNs));
        yyaxis(ax, 'right');
        plot(ax, trace.time * 1e9, trace.logicNormalized, 'Color', cLogic, ...
            'LineWidth', 0.9, 'DisplayName', 'CH3 switch monitor');
        ylabel(ax, 'logic (normalized)');
        ylim(ax, [-0.15 1.15]);
        yyaxis(ax, 'left');
        if ii == 1
            legend(ax, 'Location', 'northeast');
        end
        if ii == numel(tauNs)
            xlabel(ax, 'raw scope time (ns)');
        end
    end
    annotation(fig, 'textbox', [0.13 0.002 0.75 0.035], 'String', ...
        ['Only robust DC-offset removal and logic normalization are shown. ' ...
         'No unique Fig. 3d gate is inferable from the archive/paper.'], ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9);
    local_set_fonts(fig);

    pngFile = fullfile(outputDir, 'fig03_source_data.png');
    matFile = fullfile(outputDir, 'fig03_source_data.mat');
    print(fig, pngFile, '-dpng', '-r180');
    close(fig);

    result = struct();
    result.description = ['RAW official tau_15/25/35 acquisitions. Only a robust median ' ...
        'offset is removed from CH1/CH2 and CH3 is normalized by robust levels.'];
    result.limitation = ['The archive does not specify unique incident/reflection gates, ' ...
        'video-leakage reference periods, connector/loss compensation, or a mapping that ' ...
        'would uniquely reproduce processed Fig. 3d. No such spectrum is claimed here.'];
    result.dataDirectory = dataDir;
    result.files = files;
    result.nominalTauNs = tauNs;
    result.scope = scopes;
    result.rawDisplayTraces = traces;
    result.edgeSummary = edgeSummary;
    result.displayWindowNs = [-40 230];
    result.outputFiles = {pngFile, matFile};
    save(matFile, 'result');

    fprintf('\n[Fig. 3 source data] RAW DATA (no inferred spectral gate)\n');
    for ii = 1:numel(tauNs)
        fprintf(['  nominal tau %2d ns: switch monitor rise %.2f ns, fall %.2f ns, ' ...
            'threshold-high %.2f ns; CH1/CH2 p-p %.4f/%.4f V\n'], ...
            tauNs(ii), edgeSummary(ii).riseTimeNs, edgeSummary(ii).fallTimeNs, ...
            edgeSummary(ii).highDurationNs, edgeSummary(ii).ch1PeakToPeakV, ...
            edgeSummary(ii).ch2PeakToPeakV);
    end
    fprintf('  Limitation: %s\n', result.limitation);
    fprintf('  Wrote %s and %s\n', pngFile, matFile);
end

function trace = local_raw_trace(scope, limits)
    use = scope.time >= limits(1) & scope.time <= limits(2) & ...
        all(isfinite(scope.channels(:, 1:3)), 2);
    time = scope.time(use);
    channel = scope.channels(use, 1:3);
    if numel(time) < 4
        error('otr:EmptyDisplayWindow', 'Raw display window contains insufficient samples.');
    end
    quiet = time <= min(-5e-9, time(1) + 0.12 * (time(end) - time(1)));
    if nnz(quiet) < 4
        quiet = false(size(time));
        quiet(1:max(4, round(0.10 * numel(time)))) = true;
    end
    ch1Offset = median(channel(quiet, 1));
    ch2Offset = median(channel(quiet, 2));
    logicLow = local_quantile(channel(:, 3), 0.05);
    logicHigh = local_quantile(channel(:, 3), 0.95);

    trace = struct();
    trace.time = time;
    trace.ch1Raw = channel(:, 1);
    trace.ch2Raw = channel(:, 2);
    trace.logicRaw = channel(:, 3);
    trace.ch1OffsetV = ch1Offset;
    trace.ch2OffsetV = ch2Offset;
    trace.ch1 = channel(:, 1) - ch1Offset;
    trace.ch2 = channel(:, 2) - ch2Offset;
    trace.logicLevelsV = [logicLow logicHigh];
    trace.logicNormalized = (channel(:, 3) - logicLow) / max(eps, logicHigh - logicLow);
end

function summary = local_edge_summary(trace)
    logic = trace.logicNormalized;
    time = trace.time;
    rise = find(logic(1:end - 1) < 0.5 & logic(2:end) >= 0.5, 1, 'first') + 1;
    fallCandidates = find(logic(1:end - 1) >= 0.5 & logic(2:end) < 0.5) + 1;
    fallCandidates = fallCandidates(fallCandidates > rise);
    if isempty(rise) || isempty(fallCandidates)
        riseTime = NaN;
        fallTime = NaN;
        highDuration = NaN;
    else
        fall = fallCandidates(1);
        riseTime = local_crossing(time, logic, rise, 0.5);
        fallTime = local_crossing(time, logic, fall, 0.5);
        highDuration = fallTime - riseTime;
    end
    summary = struct('riseTimeNs', riseTime * 1e9, 'fallTimeNs', fallTime * 1e9, ...
        'highDurationNs', highDuration * 1e9, ...
        'ch1PeakToPeakV', max(trace.ch1) - min(trace.ch1), ...
        'ch2PeakToPeakV', max(trace.ch2) - min(trace.ch2));
end

function tCross = local_crossing(time, value, rightIndex, level)
    leftIndex = max(1, rightIndex - 1);
    denominator = value(rightIndex) - value(leftIndex);
    if abs(denominator) < eps
        fraction = 0;
    else
        fraction = (level - value(leftIndex)) / denominator;
    end
    fraction = max(0, min(1, fraction));
    tCross = time(leftIndex) + fraction * (time(rightIndex) - time(leftIndex));
end

function value = local_quantile(x, q)
    x = sort(x(isfinite(x)));
    position = 1 + (numel(x) - 1) * q;
    lo = floor(position);
    hi = ceil(position);
    value = x(lo) + (position - lo) * (x(hi) - x(lo));
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
