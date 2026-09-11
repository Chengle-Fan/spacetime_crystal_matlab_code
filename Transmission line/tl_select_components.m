function report = tl_select_components(model,catalog,selectCfg)
%TL_SELECT_COMPONENTS Map effective TL targets to candidate components.
%
%   report = tl_select_components(model,catalog,selectCfg)
%
% catalog may contain struct arrays (or tables) named capacitors,
% inductors, and varactors.  Each candidate uses this schema:
%   partNumber, value, toleranceFraction, esrOhm, srfHz, package, source
% value is in F for capacitors/varactors and H for inductors.  A missing
% catalog is allowed and produces target-only rows, never fabricated parts.
%
% C0 is the nominal modulated shunt contribution and is not automatically
% a purchasable capacitor.  Cpar, when nonzero, is an additive extracted
% parasitic.  Selection results must be de-embedded before they become a BOM.

if nargin ~= 3 || ~isstruct(model) || ~isscalar(model) || ...
        ~isstruct(catalog) || ~isscalar(catalog) || ...
        ~isstruct(selectCfg) || ~isscalar(selectCfg)
    error('Use tl_select_components(model,catalog,selectCfg).');
end
retainedSidebandOrder = tl_option('integer', ...
    selectCfg,'retainedSidebandOrder',3,0);
defaultMaximum = model.derived.omegaMaximumEstimate/(2*pi)+ ...
    retainedSidebandOrder*model.modulation.fmHz;
maximumFrequencyHz = tl_option('positive', ...
    selectCfg,'maximumFrequencyHz',defaultMaximum);
srfSafetyFactor = tl_option('positive',selectCfg,'srfSafetyFactor',1.5);
minimumSrfHz = maximumFrequencyHz*srfSafetyFactor;
nSamples = tl_option('integer',selectCfg,'monteCarloSamples',1000,1000);
randomSeed = tl_option('integer',selectCfg,'randomSeed',260417408,0);

capacitors = normalize_catalog(catalog,'capacitors','F');
inductors = normalize_catalog(catalog,'inductors','H');
varactors = normalize_catalog(catalog,'varactors','F');
targets = target_rows(model);
selections = repmat(empty_selection(),numel(targets),1);
for index = 1:numel(targets)
    target = targets(index);
    switch target.catalogType
        case 'capacitors', candidates = capacitors;
        case 'inductors', candidates = inductors;
        case 'varactors', candidates = varactors;
        otherwise, candidates = repmat(empty_candidate(),0,1);
    end
    selections(index) = choose_candidate( ...
        target,candidates,minimumSrfHz,maximumFrequencyHz);
end

previousRandomState = rng;
randomCleanup = onCleanup(@() rng(previousRandomState));
rng(randomSeed,'twister');
monteCarlo = tolerance_monte_carlo( ...
    model,selections,nSamples,selectCfg);
clear randomCleanup;

checks.capacitancePositive = ...
    model.modulation.capacitanceMinimum > 0;
checks.inductanceMatrixPositive = ...
    model.series.Ls > 2*abs(model.series.mutualS);
checks.maximumModeledFrequencyHz = maximumFrequencyHz;
checks.minimumRequiredSrfHz = minimumSrfHz;
checks.allSelectedSrfCompliant = all([selections.srfCompliant] | ...
    strcmp({selections.status},'target-only'));
checks.peakReverseVoltageKnown = isfield(selectCfg,'peakReverseVoltage');
if checks.peakReverseVoltageKnown
    peakReverseVoltage = tl_option('nonnegative',selectCfg,'peakReverseVoltage',0);
else
    peakReverseVoltage = NaN;
end
if isfield(selectCfg,'maximumAllowedReverseVoltage')
    maximumAllowedReverseVoltage = tl_option('positive', ...
        selectCfg,'maximumAllowedReverseVoltage',[]);
    checks.reverseVoltageCompliant = checks.peakReverseVoltageKnown && ...
        peakReverseVoltage <= maximumAllowedReverseVoltage;
else
    maximumAllowedReverseVoltage = NaN;
    checks.reverseVoltageCompliant = false;
end
checks.peakReverseVoltage = peakReverseVoltage;
checks.maximumAllowedReverseVoltage = maximumAllowedReverseVoltage;

report.targets = targets;
report.selections = selections;
report.checks = checks;
report.monteCarlo = monteCarlo;
report.modelSnapshot = model.snapshot;
report.selectionConfig = selectCfg;
report.catalogSummary = struct('capacitorCount',numel(capacitors), ...
    'inductorCount',numel(inductors),'varactorCount',numel(varactors));
report.status = 'simulation-seed component screening';
report.warning = ['Nearest nominal values are not experimental component ' ...
    'values.  Verify RF C(V), ESR/Q, SRF, package and PCB parasitics with ' ...
    'datasheets and measured impedance/S-parameters.'];
end

% -------------------------------------------------------------------------
function targets = target_rows(model)
targets = repmat(struct('name','','value',NaN,'unit','', ...
    'catalogType','','role',''),0,1);
targets(end+1) = make_target('Ls',model.series.Ls,'H','inductors', ...
    'effective series inductance per cell');
targets(end+1) = make_target('C0',model.shunt.C0,'F','capacitors', ...
    ['nominal modulated shunt capacitance excluding separately declared ' ...
    'Cpar; de-embedding required']);
if model.shunt.Cpar > 0
    targets(end+1) = make_target('Cpar',model.shunt.Cpar,'F','none', ...
        'extracted additive parasitic capacitance, not a catalog part');
end
targets(end+1) = make_target('CvarAtBias',model.paper.CvarAtBias,'F', ...
    'varactors','bias-point varactor capacitance reference');
targets(end+1) = make_target('deltaC',model.shunt.deltaC,'F','none', ...
    'required dynamic capacitance swing, not a fixed capacitor');
if strcmp(model.kind,'crow')
    targets(end+1) = make_target('L0',model.resonator.L0,'H', ...
        'inductors','parallel resonator branch inductance');
    if ~isinf(model.resonator.Cblock)
        targets(end+1) = make_target('Cblock',model.resonator.Cblock,'F', ...
            'capacitors','resonator-branch DC-block capacitor');
    end
end
end

function target = make_target(name,value,unit,catalogType,role)
target = struct('name',name,'value',value,'unit',unit, ...
    'catalogType',catalogType,'role',role);
end

% -------------------------------------------------------------------------
function candidates = normalize_catalog(catalog,fieldName,unit)
if ~isfield(catalog,fieldName) || isempty(catalog.(fieldName))
    candidates = repmat(empty_candidate(),0,1);
    return;
end
input = catalog.(fieldName);
if istable(input)
    input = table2struct(input);
end
if ~isstruct(input)
    error('catalog.%s must be a struct array or table.',fieldName);
end
required = {'partNumber','value','toleranceFraction','esrOhm','srfHz', ...
    'package','source'};
missing = required(~isfield(input,required));
if ~isempty(missing)
    error('catalog.%s is missing field(s): %s.', ...
        fieldName,strjoin(missing,', '));
end
candidates = repmat(empty_candidate(),numel(input),1);
for index = 1:numel(input)
    item = input(index);
    candidates(index).partNumber = text_value(item.partNumber,'partNumber');
    candidates(index).value = positive_value(item.value,'value');
    candidates(index).toleranceFraction = ...
        nonnegative_value(item.toleranceFraction,'toleranceFraction');
    if candidates(index).toleranceFraction > 1
        error('Candidate toleranceFraction must be expressed as a fraction.');
    end
    candidates(index).esrOhm = nonnegative_value(item.esrOhm,'esrOhm');
    candidates(index).srfHz = positive_value(item.srfHz,'srfHz');
    candidates(index).package = text_value(item.package,'package');
    candidates(index).source = text_value(item.source,'source');
    candidates(index).unit = unit;
end
end

% -------------------------------------------------------------------------
function selection = choose_candidate( ...
        target,candidates,minimumSrfHz,maximumFrequencyHz)
selection = empty_selection();
selection.name = target.name;
selection.target = target.value;
selection.unit = target.unit;
selection.role = target.role;
if strcmp(target.catalogType,'none') || isempty(candidates)
    selection.status = 'target-only';
    return;
end
compliant = [candidates.srfHz] >= minimumSrfHz;
if ~any(compliant)
    selection.status = 'no-SRF-compliant-candidate';
    return;
end
eligible = find(compliant);
[~,localIndex] = min(abs(log([candidates(eligible).value]/target.value)));
candidate = candidates(eligible(localIndex));
selection.selected = candidate.value;
selection.relativeNominalError = candidate.value/target.value-1;
selection.toleranceFraction = candidate.toleranceFraction;
selection.esrOhm = candidate.esrOhm;
if candidate.esrOhm == 0
    selection.qualityFactorAtMaximumFrequency = Inf;
elseif strcmp(candidate.unit,'H')
    selection.qualityFactorAtMaximumFrequency = ...
        2*pi*maximumFrequencyHz*candidate.value/candidate.esrOhm;
else
    selection.qualityFactorAtMaximumFrequency = ...
        1/(2*pi*maximumFrequencyHz*candidate.value*candidate.esrOhm);
end
selection.srfHz = candidate.srfHz;
selection.srfCompliant = true;
selection.partNumber = candidate.partNumber;
selection.package = candidate.package;
selection.source = candidate.source;
selection.status = 'candidate-selected';
end

% -------------------------------------------------------------------------
function output = tolerance_monte_carlo(model,selections,count,cfg)
[Ls,LsTolerance] = selected_value(selections,'Ls',model.series.Ls,0);
[C0,C0Tolerance] = selected_value(selections,'C0',model.shunt.C0,0);
if strcmp(model.kind,'crow')
    [L0,L0Tolerance] = selected_value( ...
        selections,'L0',model.resonator.L0,0);
    [Cblock,CblockTolerance] = selected_value( ...
        selections,'Cblock',model.resonator.Cblock,0);
else
    L0 = Inf;
    L0Tolerance = 0;
    Cblock = Inf;
    CblockTolerance = 0;
end
mutualTolerance = tl_option('nonnegative',cfg,'mutualToleranceFraction',0);
LsSamples = Ls*(1+LsTolerance*(2*rand(count,1)-1));
C0Samples = C0*(1+C0Tolerance*(2*rand(count,1)-1));
CtotalSamples = C0Samples+model.shunt.Cpar;
SSamples = model.series.mutualS* ...
    (1+mutualTolerance*(2*rand(count,1)-1));
valid = LsSamples > 2*abs(SSamples) & C0Samples > 0 & CtotalSamples > 0;
if strcmp(model.kind,'crow')
    L0Samples = L0*(1+L0Tolerance*(2*rand(count,1)-1));
    CblockSamples = Cblock*(1+CblockTolerance*(2*rand(count,1)-1));
    valid = valid & L0Samples > 0 & CblockSamples > 0;
else
    L0Samples = Inf(count,1);
    CblockSamples = Inf(count,1);
end
boundary = nan(count,1);
impedance = nan(count,1);
boundary(valid) = tl_static_omega(pi,LsSamples(valid),SSamples(valid), ...
    CtotalSamples(valid),L0Samples(valid),CblockSamples(valid))/(2*pi);
impedance(valid) = sqrt(LsSamples(valid)./CtotalSamples(valid));
output.sampleCount = count;
output.randomSeed = tl_option('integer',cfg,'randomSeed',260417408,0);
output.validFraction = mean(valid);
output.boundaryFrequencyHz = percentiles(boundary(valid),[5 50 95]);
output.characteristicImpedanceOhm = percentiles( ...
    impedance(valid),[5 50 95]);
if strcmp(model.kind,'crow')
    fcol = nan(count,1);
    fcol(valid) = tl_static_omega(0,LsSamples(valid),SSamples(valid), ...
        CtotalSamples(valid),L0Samples(valid),CblockSamples(valid))/(2*pi);
    output.fcolHz = percentiles(fcol(valid),[5 50 95]);
end
output.percentileLevels = [5 50 95];
output.maximumDeviceVoltage = [NaN NaN NaN];
output.maximumDeviceCurrent = [NaN NaN NaN];
output.deviceStressStatus = ...
    'unavailable without a specified source amplitude and probe definition';
end

function [value,tolerance] = selected_value(selections,name,defaultValue,defaultTolerance)
index = find(strcmp({selections.name},name),1);
if isempty(index) || ~strcmp(selections(index).status,'candidate-selected')
    value = defaultValue;
    tolerance = defaultTolerance;
else
    value = selections(index).selected;
    tolerance = selections(index).toleranceFraction;
end
end

function values = percentiles(samples,levels)
if isempty(samples)
    values = nan(size(levels));
    return;
end
samples = sort(samples(:));
coordinate = 1+(numel(samples)-1)*levels/100;
lower = floor(coordinate);
upper = ceil(coordinate);
weight = coordinate-lower;
values = samples(lower).*(1-weight(:))+samples(upper).*weight(:);
values = values(:).';
end

% -------------------------------------------------------------------------
function candidate = empty_candidate()
candidate = struct('partNumber','','value',NaN,'toleranceFraction',NaN, ...
    'esrOhm',NaN,'srfHz',NaN,'package','','source','','unit','');
end

function selection = empty_selection()
selection = struct('name','','target',NaN,'selected',NaN,'unit','', ...
    'role','','relativeNominalError',NaN,'toleranceFraction',NaN, ...
    'esrOhm',NaN,'qualityFactorAtMaximumFrequency',NaN, ...
    'srfHz',NaN,'srfCompliant',false,'partNumber','', ...
    'package','','source','','status','');
end

function value = text_value(value,label)
if isstring(value) && isscalar(value), value = char(value); end
if ~ischar(value) || size(value,1) ~= 1
    error('Candidate %s must be a text scalar.',label);
end
end

function value = positive_value(value,label)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value <= 0
    error('Candidate %s must be a positive finite real scalar.',label);
end
end

function value = nonnegative_value(value,label)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value < 0
    error('Candidate %s must be a nonnegative finite real scalar.',label);
end
end
