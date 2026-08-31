function results = reproduce_figures_2_3(options)
%REPRODUCE_FIGURES_2_3 Reproduce computable panels of paper Figs. 2 and 3.
%
% results = reproduce_figures_2_3(options)
%
% The function deliberately reuses the numerical kernels in the sibling
% "Transmission line" directory.  It adds only paper parameters, scans,
% finite-chain excitation assumptions, validation, and plotting.

if nargin < 1 || isempty(options)
    options = struct();
end
if ~isstruct(options) || ~isscalar(options)
    error('options must be a scalar struct.');
end
options = parse_options(options);

reproductionDirectory = fileparts(mfilename('fullpath'));
repositoryDirectory = fileparts(reproductionDirectory);
templateDirectory = fullfile(repositoryDirectory,'Transmission line');
if ~isfolder(templateDirectory)
    error('Transmission-line template directory not found: %s', ...
        templateDirectory);
end
oldPath = path;
pathCleanup = onCleanup(@() path(oldPath));
addpath(templateDirectory,'-begin');
templateFiles = verify_templates(templateDirectory);

if options.saveOutputs && ~isfolder(options.outputDirectory)
    mkdir(options.outputDirectory);
end

fprintf('Reproducing Fig. 2 SSPP bands with the transmission-line templates.\n');
[figure2Bands,figure2BandFigure] = calculate_figure2_bands(options);

fprintf('Reproducing Fig. 3 CROW bands with the transmission-line templates.\n');
[figure3Bands,figure3BandFigure] = calculate_figure3_bands(options);

if options.runFftBands
    fprintf(['Reconstructing finite-chain Gaussian-packet FDTD-FFT ' ...
        'responses for Figs. 2 and 3.\n']);
    [figure2FftBands,figure2FftFigure] = ...
        calculate_figure2_fft_bands(options);
    [figure3FftBands,figure3FftFigure] = ...
        calculate_figure3_fft_bands(options);
else
    figure2FftBands = struct('status','skipped by options.runFftBands');
    figure3FftBands = struct('status','skipped by options.runFftBands');
    figure2FftFigure = [];
    figure3FftFigure = [];
end

if options.runFields
    fprintf('Running finite-chain FDTD fields through t = 70 ns.\n');
    [figure2Fields,figure2FieldFigure] = ...
        calculate_figure2_fields(options);
    [figure3Fields,figure3FieldFigure] = ...
        calculate_figure3_fields(options);
else
    figure2Fields = struct('status','skipped by options.runFields');
    figure3Fields = struct('status','skipped by options.runFields');
    figure2FieldFigure = [];
    figure3FieldFigure = [];
end

results = struct();
results.paper = paper_snapshot();
results.assumptions = assumption_snapshot(options);
results.templateFiles = templateFiles;
results.options = options;
results.figure2 = struct('bands',figure2Bands, ...
    'fdtdFft',figure2FftBands,'fields',figure2Fields);
results.figure3 = struct('bands',figure3Bands, ...
    'fdtdFft',figure3FftBands,'fields',figure3Fields);
results.generatedAt = char(datetime('now','TimeZone','local', ...
    'Format','yyyy-MM-dd HH:mm:ss Z'));

if options.saveOutputs
    exportgraphics(figure2BandFigure, ...
        fullfile(options.outputDirectory,'figure2_bands.png'), ...
        'Resolution',options.imageResolution);
    exportgraphics(figure3BandFigure, ...
        fullfile(options.outputDirectory,'figure3_bands.png'), ...
        'Resolution',options.imageResolution);
    if options.runFftBands
        exportgraphics(figure2FftFigure, ...
            fullfile(options.outputDirectory,'figure2_fdtd_fft.png'), ...
            'Resolution',options.imageResolution);
        exportgraphics(figure3FftFigure, ...
            fullfile(options.outputDirectory,'figure3_fdtd_fft.png'), ...
            'Resolution',options.imageResolution);
    end
    if options.runFields
        exportgraphics(figure2FieldFigure, ...
            fullfile(options.outputDirectory,'figure2_fields.png'), ...
            'Resolution',options.imageResolution);
        exportgraphics(figure3FieldFigure, ...
            fullfile(options.outputDirectory,'figure3_fields.png'), ...
            'Resolution',options.imageResolution);
    end
    save(fullfile(options.outputDirectory,'reproduction_results.mat'), ...
        'results','-v7.3');
    fprintf('Saved reproduction outputs in %s\n',options.outputDirectory);
end

print_summary(results);
clear pathCleanup;
end

% -------------------------------------------------------------------------
function options = parse_options(input)
here = fileparts(mfilename('fullpath'));
options = struct();
options.quick = read_logical(input,'quick',false);
options.runFields = read_logical(input,'runFields',true);
options.runPhaseMap = read_logical(input,'runPhaseMap',true);
options.runFftBands = read_logical(input,'runFftBands',true);
options.saveOutputs = read_logical(input,'saveOutputs',true);
options.outputDirectory = read_text(input,'outputDirectory', ...
    fullfile(here,'outputs'));
options.imageResolution = read_positive_integer(input, ...
    'imageResolution',220);
if options.quick
    options.kCount = 61;
    options.pweOrder = 4;
    options.pweSamples = 128;
    options.tmmSlices = 128;
    options.phaseCount = 11;
    options.phaseKCount = 21;
    options.stepsPerPeriod = 64;
    options.recordEvery = 2;
    fftKCountDefault = 21;
    fftAnalysisPeriodsDefault = 6;
    fftCellCountDefault = 128;
else
    options.kCount = 121;
    options.pweOrder = 5;
    options.pweSamples = 256;
    options.tmmSlices = 256;
    options.phaseCount = 21;
    options.phaseKCount = 31;
    options.stepsPerPeriod = 128;
    options.recordEvery = 2;
    fftKCountDefault = 61;
    fftAnalysisPeriodsDefault = 12;
    fftCellCountDefault = 160;
end
options.fftKCount = read_positive_integer( ...
    input,'fftKCount',fftKCountDefault);
options.fftAnalysisPeriodCount = read_positive_integer( ...
    input,'fftAnalysisPeriodCount',fftAnalysisPeriodsDefault);
options.fftCellCount = read_positive_integer( ...
    input,'fftCellCount',fftCellCountDefault);
options.fftPulseIntensityFwhmCells = read_positive( ...
    input,'fftPulseIntensityFwhmCells',14);
options.fftZeroPaddingFactor = read_positive_integer( ...
    input,'fftZeroPaddingFactor',4);
options.fftProbeOffsets = read_integer_vector( ...
    input,'fftProbeOffsets',[-6 -3 0 3 6]);
options.fftInitialTailTolerance = read_fraction( ...
    input,'fftInitialTailTolerance',1e-4);
options.fftDynamicRangeDb = read_positive( ...
    input,'fftDynamicRangeDb',60);
if options.fftKCount < 3 || options.fftAnalysisPeriodCount < 2 || ...
        options.fftCellCount < 16
    error(['FFT reproduction requires fftKCount>=3, ' ...
        'fftAnalysisPeriodCount>=2, and fftCellCount>=16.']);
end
end

% -------------------------------------------------------------------------
function files = verify_templates(templateDirectory)
names = {'tl_build_model','tl_pwe_fourier','tl_pwe_bands', ...
    'tl_tmm_bands','tl_fdtd1d','tl_fdtd_gaussian_k_scan', ...
    'tl_fdtd_fft_bands'};
files = struct();
for index = 1:numel(names)
    resolved = which(names{index});
    expected = fullfile(templateDirectory,[names{index} '.m']);
    if ~strcmp(resolved,expected)
        error('%s resolved to %s instead of %s.',names{index}, ...
            resolved,expected);
    end
    files.(names{index}) = resolved;
end
end

% -------------------------------------------------------------------------
function [data,targetFigure] = calculate_figure2_bands(options)
paper = paper_snapshot();
baseCfg = struct('topology','sspp','allowAssumptions',true, ...
    'a',paper.fig2.a,'C0',paper.fig2.C0, ...
    'deltaC',paper.fig2.deltaC,'Ls',paper.fig2.Ls, ...
    'mutualS',0,'Rs',0,'Gp',0);
staticModel = tl_build_model(merge_struct(baseCfg, ...
    struct('deltaC',0,'fmHz',paper.fig2.fmHz(2))));
k = linspace(0,pi/staticModel.cell.a,options.kCount);
staticFrequency = static_positive_frequency(staticModel,k,1);

targetFigure = figure('Color','w','Position',[60 80 1380 390]);
layout = tiledlayout(targetFigure,1,4,'TileSpacing','compact', ...
    'Padding','compact');
staticAxes = nexttile(layout);
plot(staticAxes,k*staticModel.cell.a/pi,staticFrequency/1e6, ...
    'w:','LineWidth',2.2);
set_dark_axes(staticAxes);
xlabel(staticAxes,'k (\pi/a)');
ylabel(staticAxes,'Frequency (MHz)');
title(staticAxes,'(d) Static SSPP dispersion');
xlim(staticAxes,[0 1]);
ylim(staticAxes,[0 430]);
yline(staticAxes,paper.fig2.fr0Hz/1e6,'--', ...
    'f_{r0}=371 MHz','Color',[0.3 0.9 0.75]);

caseCells = cell(1,numel(paper.fig2.fmHz));
for index = 1:numel(paper.fig2.fmHz)
    fmHz = paper.fig2.fmHz(index);
    model = tl_build_model(merge_struct(baseCfg,struct('fmHz',fmHz)));
    [band,check] = dynamic_band(model,k,options);
    caseCells{index} = band_summary(model,band,check);

    targetAxes = nexttile(layout);
    plot_dynamic_band(targetAxes,band,model,[0 1]);
    title(targetAxes,sprintf('(%s) f_m = %.0f MHz', ...
        char('d'+index),fmHz/1e6));
    xlabel(targetAxes,'k (\pi/a)');
    if index == 1
        ylabel(targetAxes,'Re(\omega)/\Omega');
    end
end
cases = [caseCells{:}];
title(layout,'Fig. 2(d-g): calculated SSPP dispersion and expanding k-gap');

data = struct();
data.paperParameters = paper.fig2;
data.static = struct('kaOverPi',k*staticModel.cell.a/pi, ...
    'frequencyHz',staticFrequency,'model',staticModel.snapshot);
data.dynamic = cases;
data.observableNote = ['PWE lines and sparse TMM markers are theoretical; ' ...
    'the experimental FFT colormap is unavailable in the paper data.'];
end

% -------------------------------------------------------------------------
function [data,targetFigure] = calculate_figure3_bands(options)
paper = paper_snapshot();
baseCfg = struct('topology','crow','allowAssumptions',true, ...
    'a',paper.fig3.a,'Ls',paper.fig3.Ls,'mutualS',0, ...
    'C0',paper.fig3.C0,'deltaC',paper.fig3.deltaC, ...
    'L0',paper.fig3.L0,'Cblock',paper.fig3.Cblock, ...
    'fmHz',paper.fig3.fmHz,'Rs',0,'Gp',0,'R0',0);
finiteStaticModel = tl_build_model(merge_struct(baseCfg, ...
    struct('deltaC',0)));
idealStaticModel = tl_build_model(merge_struct(baseCfg, ...
    struct('deltaC',0,'Cblock',Inf)));
dynamicModel = tl_build_model(baseCfg);
 k = linspace(0,pi/dynamicModel.cell.a,options.kCount);
finiteFrequency = static_positive_frequency(finiteStaticModel,k,1);
idealFrequency = static_positive_frequency(idealStaticModel,k,1);
[band,check] = dynamic_band(dynamicModel,k,options);

if options.runPhaseMap
    fprintf('Calculating the Fig. 3(e) Q-kappa phase map.\n');
    phaseMap = crow_phase_map(options,paper.fig3);
else
    phaseMap = struct('status','skipped by options.runPhaseMap');
end

targetFigure = figure('Color','w','Position',[60 80 1160 390]);
layout = tiledlayout(targetFigure,1,3,'TileSpacing','compact', ...
    'Padding','compact');

staticAxes = nexttile(layout);
plot(staticAxes,k*dynamicModel.cell.a/pi,idealFrequency/1e6, ...
    'w:','LineWidth',2.2);
hold(staticAxes,'on');
plot(staticAxes,k*dynamicModel.cell.a/pi,finiteFrequency/1e6, ...
    '--','Color',[0.3 0.9 0.75],'LineWidth',1.6);
set_dark_axes(staticAxes);
xlabel(staticAxes,'k (\pi/a)');
ylabel(staticAxes,'Frequency (MHz)');
title(staticAxes,'(d) Static CROW dispersion');
xlim(staticAxes,[0 1]);
ylim(staticAxes,[280 420]);
legend(staticAxes,{'C_{block}=\infty (paper formula)', ...
    'C_{block}=200 pF (circuit)'},'Location','southeast', ...
    'TextColor','w','Color',[0.08 0.08 0.12]);

phaseAxes = nexttile(layout);
if options.runPhaseMap
    imagesc(phaseAxes,phaseMap.kappaRatio,phaseMap.QRatio, ...
        phaseMap.gapFraction);
    axis(phaseAxes,'xy');
    hold(phaseAxes,'on');
    contour(phaseAxes,phaseMap.kappaRatio,phaseMap.QRatio, ...
        phaseMap.gapFraction,[0.995 0.995],'w-','LineWidth',1.4);
    plot(phaseAxes,1,1,'ro','MarkerFaceColor','r','MarkerSize',6);
    colorbar(phaseAxes);
    clim(phaseAxes,[0 1]);
    xlabel(phaseAxes,'\kappa/\kappa_0');
    ylabel(phaseAxes,'Q/Q_0');
else
    axis(phaseAxes,'off');
    text(phaseAxes,0.5,0.5,'Phase map skipped', ...
        'HorizontalAlignment','center');
end
title(phaseAxes,'(e) Unstable k-span / (\pi/a)');

dynamicAxes = nexttile(layout);
plot_dynamic_band(dynamicAxes,band,dynamicModel,[0.4 0.6]);
xlabel(dynamicAxes,'k (\pi/a)');
ylabel(dynamicAxes,'Re(\omega)/\Omega');
title(dynamicAxes,'(f) f_m = 700 MHz');
title(layout,['Fig. 3(d-f): CROW dispersion, phase map, and reported ' ...
    'full-k-gap test']);

data = struct();
data.paperParameters = paper.fig3;
data.static = struct('kaOverPi',k*dynamicModel.cell.a/pi, ...
    'idealFrequencyHz',idealFrequency, ...
    'finiteBlockFrequencyHz',finiteFrequency, ...
    'idealModel',idealStaticModel.snapshot, ...
    'finiteBlockModel',finiteStaticModel.snapshot);
data.dynamic = band_summary(dynamicModel,band,check);
data.phaseMap = phaseMap;
data.staticBlockShiftHz = finiteFrequency(1)-idealFrequency(1);
end

% -------------------------------------------------------------------------
function [data,targetFigure] = calculate_figure2_fft_bands(options)
paper = paper_snapshot();
fmList = paper.fig2.fmHz;
caseCells = cell(1,numel(fmList));
targetFigure = figure('Color','w','Position',[60 80 1180 410]);
layout = tiledlayout(targetFigure,1,numel(fmList), ...
    'TileSpacing','compact','Padding','compact');

for index = 1:numel(fmList)
    cfg = struct('topology','sspp','allowAssumptions',true, ...
        'a',paper.fig2.a,'C0',paper.fig2.C0, ...
        'deltaC',paper.fig2.deltaC,'Ls',paper.fig2.Ls, ...
        'mutualS',0,'Rs',0,'Gp',0,'fmHz',fmList(index), ...
        'cellCount',options.fftCellCount);
    model = tl_build_model(cfg);
    caseCells{index} = finite_packet_fft_case( ...
        model,options,fmList(index)/2);
    targetAxes = nexttile(layout);
    plot_fdtd_fft_response(targetAxes,caseCells{index},options);
    title(targetAxes,sprintf('f_m = %.0f MHz',fmList(index)/1e6));
    if index == 1
        ylabel(targetAxes,'Re(\omega)/\Omega');
    end
end
cases = [caseCells{:}];
title(layout,['Fig. 2 finite-chain voltage-probe response; horizontal ' ...
    'axis is Gaussian source centre k_c']);
data = struct('cases',cases,'paperParameters',paper.fig2, ...
    'interpretation',['Finite-chain, finite-packet, finite-window response; ' ...
    'not the unavailable experimental colormap or exact bulk eigenvalues.']);
end

% -------------------------------------------------------------------------
function [data,targetFigure] = calculate_figure3_fft_bands(options)
paper = paper_snapshot();
cfg = struct('topology','crow','allowAssumptions',true, ...
    'a',paper.fig3.a,'Ls',paper.fig3.Ls,'mutualS',0, ...
    'C0',paper.fig3.C0,'deltaC',paper.fig3.deltaC, ...
    'L0',paper.fig3.L0,'Cblock',paper.fig3.Cblock, ...
    'fmHz',paper.fig3.fmHz,'Rs',0,'Gp',0,'R0',0, ...
    'cellCount',options.fftCellCount);
model = tl_build_model(cfg);
caseData = finite_packet_fft_case(model,options,paper.fig3.fcolHz);

targetFigure = figure('Color','w','Position',[60 80 1120 430]);
layout = tiledlayout(targetFigure,1,2,'TileSpacing','compact', ...
    'Padding','compact');
bandAxes = nexttile(layout);
plot_fdtd_fft_response(bandAxes,caseData,options);
ylabel(bandAxes,'Re(\omega)/\Omega');
title(bandAxes,'Finite-chain multi-probe V(t) FFT');

fieldAxes = nexttile(layout);
representative = caseData.scan.representativeField;
imagesc(fieldAxes,representative.node.x/model.cell.a, ...
    representative.node.t/model.modulation.period, ...
    real(representative.node.V));
axis(fieldAxes,'xy');
colorbar(fieldAxes);
xlabel(fieldAxes,'node position x/a');
ylabel(fieldAxes,'t/T');
title(fieldAxes,sprintf('Representative V(x,t), k_c a/\\pi=%.2f', ...
    representative.k*model.cell.a/pi));
title(layout,['Fig. 3 finite-chain Gaussian-packet response using the ' ...
    'reported CROW circuit parameters']);

data = struct('case',caseData,'paperParameters',paper.fig3, ...
    'interpretation',['Finite-chain voltage response. Spectral linewidth ' ...
    'does not measure Im(omega).']);
end

% -------------------------------------------------------------------------
function caseData = finite_packet_fft_case(model,options,targetFrequencyHz)
k = linspace(0,pi/model.cell.a,options.fftKCount);
T = model.modulation.period;
dt = T/options.stepsPerPeriod;
centerNode = round((model.finite.cellCount+1)/2);
probeNodes = centerNode+options.fftProbeOffsets;
if any(probeNodes < 1) || any(probeNodes > model.finite.cellCount) || ...
        numel(unique(probeNodes)) ~= numel(probeNodes)
    error('options.fftProbeOffsets do not define unique in-chain probes.');
end

scanCfg = struct();
scanCfg.kScan = k;
scanCfg.dt = dt;
scanCfg.nSteps = options.fftAnalysisPeriodCount*options.stepsPerPeriod;
scanCfg.recordEvery = 1;
scanCfg.pulseIntensityFwhmCells = ...
    options.fftPulseIntensityFwhmCells;
scanCfg.voltageAmplitude = 1;
scanCfg.centerNode = centerNode;
scanCfg.probeNodeIndices = probeNodes;
scanCfg.targetInitialFrequencyHz = targetFrequencyHz;
scanCfg.representativeK = 0.55*pi/model.cell.a;
scanCfg.precision = 'single';
scanCfg.boundaryType = 'open';
scanCfg.modulationEnabled = true;
scanCfg.modulationStart = 0;
scanCfg.modulationEnd = Inf;
scanCfg.initialTailTolerance = options.fftInitialTailTolerance;
scanCfg.requireNoBoundaryArrival = true;
scanCfg.progressEvery = max(1,ceil(options.fftKCount/5));
scan = tl_fdtd_gaussian_k_scan(model,scanCfg);

fftCfg = struct();
fftCfg.temporalPeriod = T;
fftCfg.analysisTimeRange = [0 options.fftAnalysisPeriodCount*T];
fftCfg.zeroPaddingFactor = options.fftZeroPaddingFactor;
fftCfg.dynamicRangeDb = options.fftDynamicRangeDb;
fftCfg.activeColumnRelativeThreshold = 1e-12;
fftCfg.returnProbePower = false;
fftCfg.returnProbeSignals = false;
bands = tl_fdtd_fft_bands(scan.probeSignals,scan.time,k,fftCfg);
bands.kaOverPi = k*model.cell.a/pi;

tmm = tl_tmm_bands(model,k,struct('temporalSlices',options.tmmSlices));
[~,dominantIndices] = max(bands.foldedPower,[],1);
dominantOmega = bands.foldedOmega(dominantIndices).';
dominantError = nan(1,numel(k));
for index = find(bands.activeKMask)
    distance = circular_distance(dominantOmega(index), ...
        real(tmm.omegaFolded(:,index)),model.modulation.OmegaRadPerSec);
    dominantError(index) = min(distance)/model.modulation.OmegaRadPerSec;
end
validError = dominantError(isfinite(dominantError));
validation = struct();
validation.dominantOmega = dominantOmega;
validation.dominantCircularErrorOverOmega = dominantError;
validation.medianDominantErrorOverOmega = median(validError);
validation.maximumDominantErrorOverOmega = max(validError);
validation.nativeResolutionOverOmega = ...
    bands.nativeOmegaResolutionNormalized;
validation.method = ['The dominant response peak is selected without a ' ...
    'theory target, then compared with the nearest TMM real frequency.'];

scan = rmfield(scan,'probeSignals');
caseData = struct('bands',bands,'scan',scan,'validation',validation, ...
    'model',model.snapshot,'targetInitialFrequencyHz',targetFrequencyHz);
fprintf(['  %s f_m=%.0f MHz FDTD-FFT: median/max dominant error ' ...
    '%.3g/%.3g Omega, native resolution %.3g Omega.\n'], ...
    upper(model.kind),model.modulation.fmHz/1e6, ...
    validation.medianDominantErrorOverOmega, ...
    validation.maximumDominantErrorOverOmega, ...
    validation.nativeResolutionOverOmega);
end

% -------------------------------------------------------------------------
function plot_fdtd_fft_response(targetAxes,caseData,options)
imagesc(targetAxes,caseData.bands.kaOverPi, ...
    caseData.bands.omegaOverOmega,caseData.bands.spectrumDb);
axis(targetAxes,'xy');
xlim(targetAxes,[0 1]);
ylim(targetAxes,[-0.5 0.5]);
clim(targetAxes,[-options.fftDynamicRangeDb 0]);
colormap(targetAxes,parula(256));
colorbar(targetAxes);
xlabel(targetAxes,'Gaussian source centre k_c a/\pi');
end

% -------------------------------------------------------------------------
function [band,check] = dynamic_band(model,k,options)
pweCfg = struct('Mtime',options.pweOrder,'Nt',options.pweSamples);
fourier = tl_pwe_fourier(model,pweCfg);
checkIndices = unique(round(linspace(1,numel(k),13)));
kCheck = k(checkIndices);
pweCheck = tl_pwe_bands(fourier,model,kCheck,pweCfg);
tmmCfg = struct('temporalSlices',options.tmmSlices);
tmmCheck = tl_tmm_bands(model,kCheck,tmmCfg);

if strcmp(model.kind,'crow')
    % The finite-Cblock CROW has a low-frequency parasitic pair in addition
    % to the main CROW pair.  The PWE central-harmonic heuristic can swap
    % these identities near a degeneracy, so use the direct monodromy TMM
    % for the plotted physical set and use the complete PWE replica pool
    % only as an independent coverage check.
    primary = tl_tmm_bands(model,k,tmmCfg);
    omegaPrimary = primary.omegaFolded;
    [omegaCheck,checkError] = match_pwe_to_tmm( ...
        pweCheck.omegaFoldedRaw,tmmCheck.omegaFolded, ...
        model.modulation.OmegaRadPerSec);
    primaryMethod = 'temporal TMM';
    checkMethod = 'raw PWE replica coverage';
else
    primary = tl_pwe_bands(fourier,model,k,pweCfg);
    omegaPrimary = primary.omegaSelected;
    omegaCheck = tmmCheck.omegaFolded;
    [~,checkError] = match_pwe_to_tmm(pweCheck.omegaFoldedRaw, ...
        tmmCheck.omegaFolded,model.modulation.OmegaRadPerSec);
    primaryMethod = 'PWE';
    checkMethod = 'temporal TMM with raw PWE replica coverage error';
end

band = struct();
band.kaOverPi = primary.kaOverPi;
band.omega = omegaPrimary;
band.frequencyNormalized = mod(real(omegaPrimary), ...
    model.modulation.OmegaRadPerSec)/model.modulation.OmegaRadPerSec;
band.growthHz = imag(omegaPrimary)/(2*pi);
band.growthMaximumHz = max(band.growthHz,[],1);
band.unstableMask = band.growthMaximumHz > ...
    1e-5*model.modulation.fmHz;
band.kGap = mask_extent(band.kaOverPi,band.unstableMask);
band.pweConfig = pweCfg;
band.primaryMethod = primaryMethod;

check = struct();
check.kaOverPi = tmmCheck.kaOverPi;
check.frequencyNormalized = mod(real(omegaCheck), ...
    model.modulation.OmegaRadPerSec)/model.modulation.OmegaRadPerSec;
check.growthHz = imag(omegaCheck)/(2*pi);
check.maximumCircularComplexErrorOverOmega = checkError;
check.tmmConfig = tmmCfg;
check.method = checkMethod;
band.check = check;
end

% -------------------------------------------------------------------------
function plot_dynamic_band(targetAxes,band,model,yLimits)
hold(targetAxes,'on');
set_dark_axes(targetAxes);
for mode = 1:size(band.frequencyNormalized,1)
    stable = ~band.unstableMask;
    plot(targetAxes,band.kaOverPi(stable), ...
        band.frequencyNormalized(mode,stable),'.', ...
        'Color',[0.45 0.75 1.0],'MarkerSize',8);
    plot(targetAxes,band.kaOverPi(band.unstableMask), ...
        band.frequencyNormalized(mode,band.unstableMask),'.', ...
        'Color',[1.0 0.35 0.2],'MarkerSize',9);
end
for mode = 1:size(band.check.frequencyNormalized,1)
    plot(targetAxes,band.check.kaOverPi, ...
        band.check.frequencyNormalized(mode,:),'wo', ...
        'MarkerSize',3.5,'LineWidth',0.7);
end
yline(targetAxes,0.5,':','Color',[0.85 0.85 0.85]);
xlabel(targetAxes,'k (\pi/a)');
xlim(targetAxes,[0 1]);
ylim(targetAxes,yLimits);
box(targetAxes,'on');
if all(band.unstableMask)
    label = 'full k-gap';
elseif any(band.unstableMask)
    label = sprintf('k-gap %.2f to %.2f \\pi/a', ...
        band.kGap.minimum,band.kGap.maximum);
else
    label = 'no resolved unstable k-gap';
end
text(targetAxes,0.04,yLimits(2)-0.08*diff(yLimits),label, ...
    'Color',[1 0.8 0.25],'FontSize',9,'VerticalAlignment','top');
text(targetAxes,0.04,yLimits(1)+0.06*diff(yLimits), ...
    sprintf('PWE/TMM err %.2g', ...
    band.check.maximumCircularComplexErrorOverOmega), ...
    'Color',[0.85 0.85 0.85],'FontSize',8);
if strcmp(model.kind,'crow')
    ylim(targetAxes,yLimits);
end
end

% -------------------------------------------------------------------------
function summary = band_summary(model,band,check)
summary = band;
summary.model = model.snapshot;
summary.check = check;
end

% -------------------------------------------------------------------------
function frequency = static_positive_frequency(model,k,branchFromTop)
frequency = zeros(size(k));
for index = 1:numel(k)
    A = model.functions.bulkStateMatrix(k(index),0);
    omega = 1i*eig(A);
    positive = sort(real(omega(real(omega) >= -1e-8)),'descend');
    positive = positive(positive > 1e3 | abs(k(index)) < 1e-12);
    if isempty(positive)
        frequency(index) = 0;
    else
        selected = min(branchFromTop,numel(positive));
        frequency(index) = positive(selected)/(2*pi);
    end
end
end

% -------------------------------------------------------------------------
function phaseMap = crow_phase_map(options,paper)
QRatio = linspace(0.05,1.5,options.phaseCount);
kappaRatio = linspace(0.05,3,options.phaseCount);
k = linspace(0,pi/paper.a,options.phaseKCount);
gapFraction = zeros(numel(QRatio),numel(kappaRatio));
maxGrowthHz = zeros(size(gapFraction));
omegaCol = 2*pi*paper.fcolHz;
pweCfg = struct('Mtime',2,'Nt',64);

for iq = 1:numel(QRatio)
    C0 = paper.C0*QRatio(iq);
    L0 = 1/(omegaCol^2*C0);
    for ikappa = 1:numel(kappaRatio)
        Ls = paper.Ls/kappaRatio(ikappa);
        cfg = struct('topology','crow','allowAssumptions',true, ...
            'a',paper.a,'Ls',Ls,'mutualS',0,'C0',C0, ...
            'deltaC',paper.phaseRelativeDeltaC*C0,'L0',L0, ...
            'Cblock',paper.Cblock,'fmHz',paper.fmHz, ...
            'Rs',0,'Gp',0,'R0',0);
        model = tl_build_model(cfg);
        fourier = tl_pwe_fourier(model,pweCfg);
        bands = tl_pwe_bands(fourier,model,k,pweCfg);
        growth = max(imag(bands.omegaFoldedRaw),[],1)/(2*pi);
        unstable = growth > 1e-4*paper.fmHz;
        gapFraction(iq,ikappa) = sum(unstable)/(numel(k)-1);
        if unstable(1), gapFraction(iq,ikappa) = gapFraction(iq,ikappa)-0.5/(numel(k)-1); end
        if unstable(end), gapFraction(iq,ikappa) = gapFraction(iq,ikappa)-0.5/(numel(k)-1); end
        maxGrowthHz(iq,ikappa) = max(growth);
    end
end
gapFraction = max(0,min(1,gapFraction));
phaseMap = struct('QRatio',QRatio,'kappaRatio',kappaRatio, ...
    'gapFraction',gapFraction,'maximumGrowthHz',maxGrowthHz, ...
    'pweConfig',pweCfg,'Cblock',paper.Cblock, ...
    'definition',['fraction of sampled 0 <= k <= pi/a with a positive ' ...
    'PWE growth rate above 1e-4 f_m']);
end

% -------------------------------------------------------------------------
function [data,targetFigure] = calculate_figure2_fields(options)
paper = paper_snapshot();
fmList = paper.fig2.fmHz;
caseCells = cell(1,numel(fmList));
for index = 1:numel(fmList)
    cfg = struct('topology','sspp','allowAssumptions',true, ...
        'a',paper.fig2.a,'C0',paper.fig2.C0, ...
        'deltaC',paper.fig2.deltaC,'Ls',paper.fig2.Ls, ...
        'mutualS',0,'Rs',0,'Gp',0,'fmHz',fmList(index), ...
        'cellCount',97);
    model = tl_build_model(cfg);
    sourceNode = 49;
    sourceCenter = 20e-9;
    sourceWidth = 0.75e-9;
    sourceFcn = @(time) 1e-3*exp(-0.5*((time-sourceCenter)/ ...
        sourceWidth).^2).*cos(pi*fmList(index)*(time-sourceCenter));
    field = run_finite_field(model,sourceNode,sourceFcn,options);
    xCentered = field.branch.x/model.cell.a-mean(field.branch.x/model.cell.a);
    caseCells{index} = reduce_field(field.branch.IAtNodeTime,xCentered, ...
        field.node.t,field.modulation.start,model,'series current');
end
cases = [caseCells{:}];

targetFigure = plot_field_cases(cases, ...
    arrayfun(@(f) sprintf('f_m = %.0f MHz',f/1e6),fmList, ...
    'UniformOutput',false),[-48 48], ...
    'Fig. 2(h-j): SSPP |I|^2 proxy for |H_z|^2');
data = struct('cases',cases,'paperParameters',paper.fig2, ...
    'sourceAssumption',['1 mA Gaussian current pulse at x/a=0, t=20 ns, ' ...
    'sigma=0.75 ns, carrier=f_m/2; matched 50 ohm ends']);
end

% -------------------------------------------------------------------------
function [data,targetFigure] = calculate_figure3_fields(options)
paper = paper_snapshot();
sourceCoordinates = [-7 0 9];
sourceNodes = sourceCoordinates+17;
caseCells = cell(1,numel(sourceNodes));
for index = 1:numel(sourceNodes)
    cfg = struct('topology','crow','allowAssumptions',true, ...
        'a',paper.fig3.a,'Ls',paper.fig3.Ls,'mutualS',0, ...
        'C0',paper.fig3.C0,'deltaC',paper.fig3.deltaC, ...
        'L0',paper.fig3.L0,'Cblock',paper.fig3.Cblock, ...
        'fmHz',paper.fig3.fmHz,'Rs',0,'Gp',0,'R0',0, ...
        'cellCount',33);
    model = tl_build_model(cfg);
    sourceCenter = 12e-9;
    sourceWidth = 0.75e-9;
    sourceFcn = @(time) 1e-3*exp(-0.5*((time-sourceCenter)/ ...
        sourceWidth).^2).*cos(pi*paper.fig3.fmHz* ...
        (time-sourceCenter));
    field = run_finite_field(model,sourceNodes(index),sourceFcn,options);
    xCentered = field.resonator.x/model.cell.a-16;
    caseCells{index} = reduce_field(field.resonator.IHalf,xCentered, ...
        field.node.t,field.modulation.start,model,'resonator current');
    caseCells{index}.sourceCoordinate = sourceCoordinates(index);
end
cases = [caseCells{:}];

labels = arrayfun(@(x) sprintf('source x/a = %d',x), ...
    sourceCoordinates,'UniformOutput',false);
targetFigure = plot_field_cases(cases,labels,[-16 16], ...
    'Fig. 3(g-i): CROW |I_0|^2 proxy for |H_z|^2');
data = struct('cases',cases,'paperParameters',paper.fig3, ...
    'sourceCoordinates',sourceCoordinates, ...
    'sourceAssumption',['1 mA Gaussian current pulse at t=12 ns, ' ...
    'sigma=0.75 ns, carrier=f_m/2; matched 50 ohm ends']);
end

% -------------------------------------------------------------------------
function field = run_finite_field(model,sourceNode,sourceFcn,options)
dt = model.modulation.period/options.stepsPerPeriod;
modulationStart = round(32e-9/dt)*dt;
nSteps = round(70e-9/dt);
cfg = struct();
cfg.dt = dt;
cfg.nSteps = nSteps;
cfg.recordEvery = options.recordEvery;
cfg.boundaryType = 'matched';
cfg.modulationEnabled = true;
cfg.modulationStart = modulationStart;
cfg.modulationEnd = Inf;
cfg.precision = 'single';
cfg.source = struct('type','current','node',sourceNode, ...
    'waveformFcn',sourceFcn);
field = tl_fdtd1d(model,cfg);
end

% -------------------------------------------------------------------------
function reduced = reduce_field(signal,x,time,modulationStart,model,label)
signal = double(signal);
power = abs(signal).^2;
referenceRows = time < modulationStart & time > modulationStart-8e-9;
referencePower = max(power(referenceRows,:),[],'all');
if isempty(referencePower) || referencePower <= realmin
    referencePower = max(power,[],'all');
end
relativePower = power/max(referencePower,realmin);
displayPower = log10(max(relativePower,1e-8));

finalPower = power(end,:);
if max(finalPower) > 0
    finalPower = finalPower/max(finalPower);
end
weightSum = sum(finalPower);
if weightSum > 0
    center = sum(x.*finalPower)/weightSum;
    sigma = sqrt(sum((x-center).^2.*finalPower)/weightSum);
else
    center = 0;
    sigma = 1;
end
sigma = max(sigma,0.35);
gaussianFit = exp(-0.5*((x-center)/sigma).^2);
gaussianFit = gaussianFit/max(gaussianFit);

fitRows = time >= max(modulationStart+8e-9,time(end)-20e-9);
amplitude = sqrt(sum(power,2));
valid = fitRows(:) & amplitude > 0;
fit = polyfit(time(valid),log(amplitude(valid)),1);
growthHz = fit(1)/(2*pi);

reduced = struct();
reduced.xOverA = x;
reduced.timeSeconds = time;
reduced.log10PowerRelativeToPrePump = single(displayPower);
reduced.finalNormalizedPower = finalPower;
reduced.gaussianFit = gaussianFit;
reduced.gaussianCenterOverA = center;
reduced.gaussianSigmaOverA = sigma;
reduced.fittedLateTimeGrowthHz = growthHz;
reduced.referencePower = referencePower;
reduced.modulationStartSeconds = modulationStart;
reduced.observable = label;
reduced.model = model.snapshot;
end

% -------------------------------------------------------------------------
function targetFigure = plot_field_cases(cases,labels,xLimits,figureTitle)
targetFigure = figure('Color','w','Position',[60 80 1120 620]);
layout = tiledlayout(targetFigure,2,numel(cases), ...
    'TileSpacing','compact','Padding','compact');
allMaximum = -Inf;
for index = 1:numel(cases)
    allMaximum = max(allMaximum,max( ...
        cases(index).log10PowerRelativeToPrePump,[],'all'));
end
colorLimits = [-6,max(0.5,allMaximum)];

for index = 1:numel(cases)
    lineAxes = nexttile(layout,index);
    plot(lineAxes,cases(index).xOverA, ...
        cases(index).finalNormalizedPower,'Color',[1 0.55 0.1], ...
        'LineWidth',1.0);
    hold(lineAxes,'on');
    plot(lineAxes,cases(index).xOverA,cases(index).gaussianFit, ...
        'r-','LineWidth',1.5);
    xlim(lineAxes,xLimits);
    ylim(lineAxes,[0 1.05]);
    grid(lineAxes,'on');
    title(lineAxes,sprintf('%s, growth %.2f MHz',labels{index}, ...
        cases(index).fittedLateTimeGrowthHz/1e6));
    if index == 1
        ylabel(lineAxes,'normalized |I|^2 at 70 ns');
    end

    fieldAxes = nexttile(layout,numel(cases)+index);
    imagesc(fieldAxes,cases(index).xOverA, ...
        cases(index).timeSeconds/1e-9, ...
        cases(index).log10PowerRelativeToPrePump);
    axis(fieldAxes,'xy');
    xlim(fieldAxes,xLimits);
    ylim(fieldAxes,[5 70]);
    clim(fieldAxes,colorLimits);
    colormap(fieldAxes,hot(256));
    hold(fieldAxes,'on');
    yline(fieldAxes,cases(index).modulationStartSeconds/1e-9, ...
        'w--','LineWidth',1.0);
    xlabel(fieldAxes,'x/a');
    if index == 1
        ylabel(fieldAxes,'t (ns)');
    end
end
colorbar(fieldAxes,'Location','eastoutside');
title(layout,figureTitle);
end

% -------------------------------------------------------------------------
function extent = mask_extent(coordinate,mask)
if any(mask)
    extent = struct('minimum',min(coordinate(mask)), ...
        'maximum',max(coordinate(mask)), ...
        'fraction',sum(mask)/numel(mask));
else
    extent = struct('minimum',NaN,'maximum',NaN,'fraction',0);
end
end

% -------------------------------------------------------------------------
function [matched,errorValue] = match_pwe_to_tmm(pwePool,tmmSet,Omega)
matched = complex(zeros(size(tmmSet)));
errorValue = 0;
for ik = 1:size(tmmSet,2)
    candidates = pwePool(:,ik);
    for mode = 1:size(tmmSet,1)
        target = tmmSet(mode,ik);
        realError = abs(mod(real(candidates)-real(target)+Omega/2, ...
            Omega)-Omega/2);
        imagError = abs(imag(candidates)-imag(target));
        distance = hypot(realError,imagError);
        [minimum,index] = min(distance);
        matched(mode,ik) = candidates(index);
        errorValue = max(errorValue,minimum/Omega);
    end
end
end

% -------------------------------------------------------------------------
function set_dark_axes(targetAxes)
targetAxes.Color = [0.04 0.04 0.08];
targetAxes.XColor = [0.85 0.85 0.85];
targetAxes.YColor = [0.85 0.85 0.85];
grid(targetAxes,'on');
targetAxes.GridColor = [0.45 0.45 0.5];
targetAxes.GridAlpha = 0.2;
box(targetAxes,'on');
end

% -------------------------------------------------------------------------
function paper = paper_snapshot()
paper.source = 'arXiv:2604.17408v1, pp. 5-8 and Figs. 2-3';
paper.fig2 = struct();
paper.fig2.a = 4e-3;
paper.fig2.reportedStripSpacing = 4e-3;
paper.fig2.C0 = 19.8e-12;
paper.fig2.deltaC = 2.38e-12;
paper.fig2.fmHz = [575 675 725]*1e6;
paper.fig2.fr0Hz = 371e6;
paper.fig2.Vdc = 14;
paper.fig2.Cvar = 12e-12;
paper.fig2.Ls = 4/((2*pi*paper.fig2.fr0Hz)^2*paper.fig2.C0);
paper.fig2.mutualInductance = 0;

paper.fig3 = struct();
paper.fig3.a = 12e-3;
paper.fig3.reportedStripSpacing = 12e-3;
paper.fig3.C0 = 27.8e-12;
paper.fig3.deltaC = 5e-12;
paper.fig3.phaseRelativeDeltaC = 0.18;
paper.fig3.fmHz = 700e6;
paper.fig3.fcolHz = 315e6;
paper.fig3.L0 = 9.2e-9;
paper.fig3.Ls = 165e-9;
paper.fig3.Cblock = 200e-12;
paper.fig3.Vdc = 10.8;
paper.fig3.Cvar = 17e-12;
paper.fig3.Q0 = 20.7;
end

% -------------------------------------------------------------------------
function assumptions = assumption_snapshot(options)
assumptions = struct();
assumptions.aEqualsReportedStripSpacing = true;
assumptions.mutualInductanceH = 0;
assumptions.seriesResistanceOhm = 0;
assumptions.shuntConductanceSiemens = 0;
assumptions.resonatorResistanceOhm = 0;
assumptions.portImpedanceOhm = 50;
assumptions.fig2NodeCount = 97;
assumptions.fig3NodeCount = 33;
assumptions.fieldNodeCounts = struct('figure2',97,'figure3',33);
assumptions.fftNodeCount = options.fftCellCount;
assumptions.fftPulseIntensityFwhmCells = ...
    options.fftPulseIntensityFwhmCells;
assumptions.fftProbeOffsets = options.fftProbeOffsets;
assumptions.fftAnalysisPeriodCount = options.fftAnalysisPeriodCount;
assumptions.modulationStartSeconds = 32e-9;
assumptions.finalTimeSeconds = 70e-9;
assumptions.note = ['The PDF does not report these quantities completely; ' ...
    'they are simulation assumptions, not experimental parameters.'];
end

% -------------------------------------------------------------------------
function print_summary(results)
fprintf('\nReproduction diagnostics:\n');
for index = 1:numel(results.figure2.bands.dynamic)
    item = results.figure2.bands.dynamic(index);
    fprintf(['  Fig. 2, f_m=%.0f MHz: k-gap %.3f...%.3f pi/a, ' ...
        'max growth %.3f MHz, PWE/TMM error %.3g Omega.\n'], ...
        item.model.modulation.fmHz/1e6,item.kGap.minimum, ...
        item.kGap.maximum,max(item.growthMaximumHz)/1e6, ...
        item.check.maximumCircularComplexErrorOverOmega);
end
item = results.figure3.bands.dynamic;
fprintf(['  Fig. 3, f_m=700 MHz: unstable fraction %.3f, ' ...
    'max growth %.3f MHz, PWE/TMM error %.3g Omega.\n'], ...
    item.kGap.fraction,max(item.growthMaximumHz)/1e6, ...
    item.check.maximumCircularComplexErrorOverOmega);
fprintf('  Fig. 3 finite-Cblock static k=0 shift: %.3f MHz.\n', ...
    results.figure3.bands.staticBlockShiftHz/1e6);
if results.options.runFftBands
    for index = 1:numel(results.figure2.fdtdFft.cases)
        item = results.figure2.fdtdFft.cases(index);
        fprintf(['  Fig. 2 finite-chain FFT, f_m=%.0f MHz: median/max ' ...
            'dominant-peak error %.3g/%.3g Omega.\n'], ...
            item.model.modulation.fmHz/1e6, ...
            item.validation.medianDominantErrorOverOmega, ...
            item.validation.maximumDominantErrorOverOmega);
    end
    item = results.figure3.fdtdFft.case;
    fprintf(['  Fig. 3 finite-chain FFT: median/max dominant-peak ' ...
        'error %.3g/%.3g Omega.\n'], ...
        item.validation.medianDominantErrorOverOmega, ...
        item.validation.maximumDominantErrorOverOmega);
end
end

% -------------------------------------------------------------------------
function output = merge_struct(first,second)
output = first;
names = fieldnames(second);
for index = 1:numel(names)
    output.(names{index}) = second.(names{index});
end
end

% -------------------------------------------------------------------------
function distance = circular_distance(first,second,period)
distance = abs(mod(first-second+period/2,period)-period/2);
end

% -------------------------------------------------------------------------
function value = read_logical(input,name,defaultValue)
if isfield(input,name) && ~isempty(input.(name))
    value = input.(name);
else
    value = defaultValue;
end
if ~islogical(value) || ~isscalar(value)
    error('options.%s must be a logical scalar.',name);
end
end

function value = read_text(input,name,defaultValue)
if isfield(input,name) && ~isempty(input.(name))
    value = input.(name);
else
    value = defaultValue;
end
if isstring(value) && isscalar(value)
    value = char(value);
end
if ~ischar(value) || size(value,1) ~= 1
    error('options.%s must be a text scalar.',name);
end
end

function value = read_positive_integer(input,name,defaultValue)
if isfield(input,name) && ~isempty(input.(name))
    value = input.(name);
else
    value = defaultValue;
end
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value <= 0 || value ~= round(value)
    error('options.%s must be a positive integer.',name);
end
end

function value = read_positive(input,name,defaultValue)
if isfield(input,name) && ~isempty(input.(name))
    value = input.(name);
else
    value = defaultValue;
end
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value <= 0
    error('options.%s must be a positive finite real scalar.',name);
end
end

function value = read_fraction(input,name,defaultValue)
value = read_positive(input,name,defaultValue);
if value >= 1
    error('options.%s must be smaller than one.',name);
end
end

function value = read_integer_vector(input,name,defaultValue)
if isfield(input,name) && ~isempty(input.(name))
    value = input.(name);
else
    value = defaultValue;
end
if ~isnumeric(value) || isempty(value) || ~isvector(value) || ...
        ~isreal(value) || any(~isfinite(value)) || ...
        any(value ~= round(value)) || numel(unique(value)) ~= numel(value)
    error('options.%s must contain unique finite integers.',name);
end
value = value(:).';
end
