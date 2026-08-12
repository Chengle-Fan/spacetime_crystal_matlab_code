function summary = run_all_reproductions(varargin)
%RUN_ALL_REPRODUCTIONS Reproduce all numerical phenomena in Moussa et al.
%   SUMMARY = RUN_ALL_REPRODUCTIONS() executes the official-source-data
%   plots, analytic/circuit figures, the Kirchhoff/MNA demo, and the
%   independent numerical validation suite.  Figures and MAT files are
%   written below this file in output/.
%
%   RUN_ALL_REPRODUCTIONS('quick') skips the three large oscilloscope-data
%   readers and the paper-scale MNA movie-like field map.  It still runs all
%   analytic figures and validation.
%
%   Name/value options:
%     'IncludeSourceData'       logical (default true)
%     'IncludeCircuit'          logical (default true)
%     'RunValidation'           logical (default true)
%     'CloseGeneratedFigures'  logical (default true)
%     'ContinueOnError'         logical (default false)

paths = otr_setup();
options = local_options(varargin{:});

jobs = local_jobs(options);
entries = repmat(struct('name','','phenomenon','','status','pending', ...
    'elapsedSeconds',NaN,'message',''),numel(jobs),1);

fprintf('\nObservation of temporal reflection -- full reproduction\n');
fprintf('Output: %s\n',paths.output);
fprintf('Jobs: %d (official data: %s, MNA: %s, validation: %s)\n\n', ...
    numel(jobs),local_onoff(options.IncludeSourceData), ...
    local_onoff(options.IncludeCircuit), ...
    local_onoff(options.RunValidation));

runClock = tic;
for jobIndex = 1:numel(jobs)
    job = jobs(jobIndex);
    entries(jobIndex).name = job.name;
    entries(jobIndex).phenomenon = job.phenomenon;
oldFigures = findall(groot,'Type','figure');
    fprintf('[%02d/%02d] %-34s ',jobIndex,numel(jobs),job.name);
    jobClock = tic;
    try
        job.call();
        entries(jobIndex).status = 'passed';
        entries(jobIndex).elapsedSeconds = toc(jobClock);
        fprintf('PASS  %7.2f s\n',entries(jobIndex).elapsedSeconds);
    catch exception
        entries(jobIndex).status = 'failed';
        entries(jobIndex).elapsedSeconds = toc(jobClock);
        entries(jobIndex).message = exception.message;
        fprintf('FAIL  %7.2f s\n',entries(jobIndex).elapsedSeconds);
        fprintf('         %s\n',exception.message);
        local_save_manifest(paths,options,entries,toc(runClock));
        if ~options.ContinueOnError
            rethrow(exception);
        end
    end

    if options.CloseGeneratedFigures
        newFigures = setdiff(findall(groot,'Type','figure'),oldFigures);
        if ~isempty(newFigures)
            close(newFigures);
        end
    end
end

summary = local_save_manifest(paths,options,entries,toc(runClock));
failed = strcmp({entries.status},'failed');
fprintf('\nCompleted %d/%d jobs in %.1f s; failures: %d.\n', ...
    sum(~failed),numel(jobs),summary.elapsedSeconds,sum(failed));
fprintf('Manifest: %s\n\n',fullfile(paths.output,'run_manifest.mat'));

if any(failed)
    error('otr:ReproductionFailed','%d reproduction job(s) failed.', ...
        sum(failed));
end
end

function jobs = local_jobs(options)
jobs = struct('name',{},'phenomenon',{},'call',{});

if options.IncludeSourceData
    jobs(end+1) = local_job('fig01_source_data', ...
        'Measured phase, group velocity, and time-reflection trace', ...
        @fig01_source_data);
    jobs(end+1) = local_job('fig02_source_data', ...
        'Measured redshift/blueshift and frequency sweeps', ...
        @fig02_source_data);
    jobs(end+1) = local_job('fig03_source_data', ...
        'Raw temporal-slab oscilloscope records', ...
        @fig03_source_data);
end

jobs(end+1) = local_job('fig04_boundary_conditions', ...
    'Charge-continuous ON and voltage-continuous OFF boundaries', ...
    @fig04_boundary_conditions);
jobs(end+1) = local_job('fig05_ideal_time_reflection', ...
    'Time reversal, negative polarity, and conserved wavenumber', ...
    @fig05_ideal_time_reflection);
jobs(end+1) = local_job('fig06_temporal_slab_theory', ...
    'Four causal paths and temporal Fabry-Perot interference', ...
    @fig06_temporal_slab_theory);
jobs(end+1) = local_job('fig07_finite_switching_and_homogeneity', ...
    'Finite-rise and switching-synchronization penalties', ...
    @fig07_finite_switching_and_homogeneity);
jobs(end+1) = local_job('fig08_inverted_slab_and_leakage', ...
    'Inverted slab loss and control-line leakage', ...
    @fig08_inverted_slab_and_leakage);
jobs(end+1) = local_job('fig09_tlm_design', ...
    'Loaded transmission-line unit-cell design',@fig09_tlm_design);
jobs(end+1) = local_job('fig10_scattering_spectra', ...
    'Single-interface scattering spectra',@fig10_scattering_spectra);
jobs(end+1) = local_job('fig11_slab_wavepackets', ...
    'Temporal-slab wave packets and duration sweep', ...
    @fig11_slab_wavepackets);

if options.IncludeCircuit
    jobs(end+1) = local_job('demo_circuit_time_interface', ...
        'Kirchhoff/MNA time interface and exact state remapping', ...
        @demo_circuit_time_interface);
end
if options.RunValidation
    jobs(end+1) = local_job('otr_validate', ...
        'Independent analytic, data, and MNA acceptance tests', ...
        @otr_validate);
end
end

function job = local_job(name,phenomenon,call)
job = struct('name',name,'phenomenon',phenomenon,'call',call);
end

function summary = local_save_manifest(paths,options,entries,elapsedSeconds)
summary = struct();
summary.generatedAt = char(datetime('now','Format','yyyyMMdd''T''HHmmss'));
summary.matlabVersion = version;
summary.root = paths.reproduction;
summary.options = options;
summary.entries = entries;
summary.elapsedSeconds = elapsedSeconds;
finished = ~strcmp({entries.status},'pending');
summary.allPassed = all(finished) && ...
    all(strcmp({entries.status},'passed'));
save(fullfile(paths.output,'run_manifest.mat'),'summary');
end

function options = local_options(varargin)
options = struct('IncludeSourceData',true,'IncludeCircuit',true, ...
    'RunValidation',true,'CloseGeneratedFigures',true, ...
    'ContinueOnError',false);

if isscalar(varargin) && (ischar(varargin{1}) || ...
        (isstring(varargin{1}) && isscalar(varargin{1})))
    mode = lower(char(varargin{1}));
    if strcmp(mode,'quick')
        options.IncludeSourceData = false;
        options.IncludeCircuit = false;
        return;
    elseif strcmp(mode,'all')
        return;
    end
end

if mod(numel(varargin),2) ~= 0
    error('otr:InvalidOptions','Options must be name/value pairs or quick.');
end
names = fieldnames(options);
for optionIndex = 1:2:numel(varargin)
    name = char(varargin{optionIndex});
    match = find(strcmpi(name,names),1);
    if isempty(match)
        error('otr:InvalidOption','Unknown option "%s".',name);
    end
    value = varargin{optionIndex+1};
    if ~(islogical(value) && isscalar(value))
        error('otr:InvalidOption','%s must be a logical scalar.',names{match});
    end
    options.(names{match}) = value;
end
end

function label = local_onoff(value)
if value
    label = 'on';
else
    label = 'off';
end
end
