function model = tl_build_model(physicalCfg)
%TL_BUILD_MODEL Build one traceable SSPP or CROW transmission-line model.
%
%   model = tl_build_model(physicalCfg)
%
% The returned model is the single physical-parameter source used by PWE,
% temporal TMM, finite-chain FDTD, and both FFT paths.  Internal units are
% SI.  The seed values come only from the experimental platform described
% in arXiv:2604.17408v1; parameters absent from that PDF are explicit
% assumptions and require physicalCfg.allowAssumptions=true.
%
% Required field:
%   topology            'sspp' or 'crow'
%
% Common optional overrides:
%   allowAssumptions    logical scalar (default false)
%   a, Ls, mutualS, Rs, C0, deltaC, Gp
%   fmHz, modulationPhase, modulationType, dutyCycle
%   Zsource, Zload, cellCount, phaseByCell, amplitudeScaleByCell
%
% CROW-only overrides:
%   L0, R0, Cblock
%
% The prescribed capacitance is linear and independent of signal voltage.
% Device C-V nonlinearity must be calibrated separately before replacing
% model.modulation.bulkFcn/finiteFcn.

if nargin ~= 1 || ~isstruct(physicalCfg) || ~isscalar(physicalCfg)
    error('Use tl_build_model with one scalar physicalCfg struct.');
end
topology = read_text(physicalCfg,'topology','',{'sspp','crow'});
allowAssumptions = read_logical(physicalCfg,'allowAssumptions',false);

seed = paper_seed(topology);
assumedFields = {'a','Ls','mutualS','Rs','Gp','Zsource','Zload'};
if strcmp(topology,'crow')
    assumedFields = [assumedFields,{'R0'}];
end
assumptionPresent = false(size(assumedFields));
for assumptionIndex = 1:numel(assumedFields)
    assumptionPresent(assumptionIndex) = ...
        isfield(physicalCfg,assumedFields{assumptionIndex});
end
missingAssumptions = assumedFields(~assumptionPresent);
if ~allowAssumptions && ~isempty(missingAssumptions)
    error(['The paper does not identify the following hardware values: %s. ' ...
        'Provide them explicitly or set allowAssumptions=true for a ' ...
        'simulation seed.'],strjoin(missingAssumptions,', '));
end

a = read_positive(physicalCfg,'a',seed.a);
Ls = read_positive(physicalCfg,'Ls',seed.Ls);
mutualS = read_finite_real(physicalCfg,'mutualS',seed.mutualS);
Rs = read_nonnegative(physicalCfg,'Rs',seed.Rs);
C0 = read_positive(physicalCfg,'C0',seed.C0);
deltaC = read_nonnegative(physicalCfg,'deltaC',seed.deltaC);
Gp = read_nonnegative(physicalCfg,'Gp',seed.Gp);
Cpar = read_nonnegative(physicalCfg,'Cpar',0);
fmHz = read_positive(physicalCfg,'fmHz',seed.fmHz);
phase = read_finite_real(physicalCfg,'modulationPhase',0);
modulationType = read_text(physicalCfg,'modulationType','sinusoidal', ...
    {'sinusoidal','square'});
dutyCycle = read_finite_real(physicalCfg,'dutyCycle',0.5);
if dutyCycle <= 0 || dutyCycle >= 1
    error('physicalCfg.dutyCycle must lie strictly between zero and one.');
end
Zsource = read_positive(physicalCfg,'Zsource',50);
Zload = read_positive(physicalCfg,'Zload',50);
cellCount = read_integer(physicalCfg,'cellCount',64,3);

if strcmp(topology,'crow')
    L0 = read_positive(physicalCfg,'L0',seed.L0);
    R0 = read_nonnegative(physicalCfg,'R0',seed.R0);
    Cblock = read_positive_or_inf(physicalCfg,'Cblock',seed.Cblock);
else
    L0 = NaN;
    R0 = NaN;
    Cblock = NaN;
end

if Ls <= 2*abs(mutualS)
    error('Require Ls > 2*abs(mutualS) so the periodic inductance is positive.');
end
if deltaC >= C0
    error('Require deltaC < C0 so the prescribed capacitance remains positive.');
end

phaseByCell = read_cell_array(physicalCfg,'phaseByCell',0,cellCount,true);
amplitudeScaleByCell = read_cell_array(physicalCfg, ...
    'amplitudeScaleByCell',1,cellCount,false);
if any(amplitudeScaleByCell < 0)
    error('physicalCfg.amplitudeScaleByCell must be nonnegative.');
end
if max(amplitudeScaleByCell)*deltaC >= C0
    error('Per-cell modulation would make at least one capacitance nonpositive.');
end

Omega = 2*pi*fmHz;
T = 1/fmHz;
capMin = C0-deltaC;
capMax = C0+deltaC;
if strcmp(modulationType,'square')
    capMin = C0-deltaC;
    capMax = C0+deltaC;
end

if strcmp(topology,'sspp')
    stateDimension = 2;
    stateLabels = {'Q','PhiSeries'};
    observableLabels = {'V','ISeries'};
elseif isinf(Cblock)
    stateDimension = 3;
    stateLabels = {'Q','PhiSeries','Phi0'};
    observableLabels = {'V','ISeries','I0'};
else
    stateDimension = 4;
    stateLabels = {'Q','PhiSeries','Phi0','Qblock'};
    observableLabels = {'V','ISeries','I0','Vblock'};
end

model.kind = topology;
model.version = '0.1.0';
model.units = 'SI';
model.phaseConvention = 'exp(i*k*x-i*omega*t)';
model.cell = struct('a',a,'reportedStripSpacing',seed.reportedW);
model.series = struct('Ls',Ls,'Rs',Rs,'mutualS',mutualS);
model.shunt = struct('C0',C0,'deltaC',deltaC,'Gp',Gp,'Cpar',Cpar);
model.resonator = struct('L0',L0,'R0',R0,'Cblock',Cblock);
model.ports = struct('Zsource',Zsource,'Zload',Zload);
model.bulk = struct('stateDimension',stateDimension, ...
    'stateLabels',{stateLabels},'observableLabels',{observableLabels});
model.finite = struct('cellCount',cellCount,'phaseByCell',phaseByCell, ...
    'amplitudeScaleByCell',amplitudeScaleByCell);
model.modulation = struct('type',modulationType,'fmHz',fmHz, ...
    'OmegaRadPerSec',Omega,'period',T,'phase',phase, ...
    'dutyCycle',dutyCycle,'capacitanceMinimum',capMin, ...
    'capacitanceMaximum',capMax);
model.paper = struct('Vdc',seed.Vdc,'CvarAtBias',seed.Cvar, ...
    'reportedDeltaC',seed.deltaC,'referenceFrequencyHz',seed.referenceHz, ...
    'source','arXiv:2604.17408v1');

model.provenance = make_provenance(seed,physicalCfg,topology);

model.functions.capacitanceBulk = @(t) bulk_capacitance( ...
    t,C0,deltaC,Omega,phase,modulationType,dutyCycle);
model.functions.capacitanceFinite = @(cellIndex,t) finite_capacitance( ...
    cellIndex,t,C0,deltaC,Omega,phase,modulationType,dutyCycle, ...
    phaseByCell,amplitudeScaleByCell);
model.functions.bulkMatrices = @(k) bulk_matrices( ...
    k,a,Ls,mutualS,Rs,Gp,topology,L0,R0,Cblock);
model.functions.bulkStateMatrix = @(k,t) bulk_state_matrix( ...
    k,t,a,Ls,mutualS,Rs,Gp,topology,L0,R0,Cblock, ...
    C0,deltaC,Omega,phase,modulationType,dutyCycle);
model.functions.observableMatrix = @(k,t) observable_matrix( ...
    k,t,a,Ls,mutualS,topology,L0,Cblock, ...
    C0,deltaC,Omega,phase,modulationType,dutyCycle);

omegaMax = estimate_omega_max(model,capMin);
boundaryAcoustic = 2/sqrt((Ls-2*mutualS)*C0);
if strcmp(topology,'crow')
    omegaCol = 1/sqrt(L0*C0);
    boundaryProvisional = sqrt(omegaCol^2+boundaryAcoustic^2);
else
    omegaCol = 0;
    boundaryProvisional = boundaryAcoustic;
end
model.derived = struct();
model.derived.omegaMaximumEstimate = omegaMax;
model.derived.maximumLeapfrogDt = 1.8/omegaMax;
model.derived.fcolHz = omegaCol/(2*pi);
model.derived.provisionalBoundaryFrequencyHz = ...
    boundaryProvisional/(2*pi);
model.derived.Zreference = sqrt(Ls/C0);
model.derived.longWaveVelocity = a/sqrt((Ls+2*mutualS)*C0);
model.derived.effectiveBoundaryInductance = Ls-2*mutualS;

snapshot = model;
snapshot = rmfield(snapshot,'functions');
model.snapshot = snapshot;
end

% -------------------------------------------------------------------------
function seed = paper_seed(topology)
if strcmp(topology,'sspp')
    seed.reportedW = 4e-3;
    seed.a = 4e-3;
    seed.Vdc = 14;
    seed.Cvar = 12e-12;
    seed.C0 = 19.8e-12;
    seed.deltaC = 2.38e-12;
    seed.fmHz = 675e6;
    seed.referenceHz = 371e6;
    seed.Ls = 4/((2*pi*seed.referenceHz)^2*seed.C0);
    seed.mutualS = 0;
    seed.Rs = 0;
    seed.Gp = 0;
    seed.L0 = NaN;
    seed.R0 = NaN;
    seed.Cblock = NaN;
else
    seed.reportedW = 12e-3;
    seed.a = 12e-3;
    seed.Vdc = 10.8;
    seed.Cvar = 17e-12;
    seed.C0 = 27.8e-12;
    seed.deltaC = 5e-12;
    seed.fmHz = 700e6;
    seed.referenceHz = 315e6;
    seed.Ls = 165e-9;
    seed.mutualS = 0;
    seed.Rs = 0;
    seed.Gp = 0;
    seed.L0 = 9.2e-9;
    seed.R0 = 0;
    seed.Cblock = 200e-12;
end
end

% -------------------------------------------------------------------------
function provenance = make_provenance(seed,cfg,topology)
provenance.a = provenance_item(seed.a,'paper w used as a', ...
    'PDF pp. 5/7','assumption','unknown/not-reported');
provenance.Ls = provenance_item(seed.Ls,'paper-derived/reference', ...
    'PDF p. 5 or Fig. 3e','paper-derived','unknown/not-reported');
provenance.mutualS = provenance_item(seed.mutualS,'not reported', ...
    'none','assumption','unknown/not-reported');
provenance.Rs = provenance_item(seed.Rs,'not reported', ...
    'none','assumption','unknown/not-reported');
provenance.C0 = provenance_item(seed.C0,'effective capacitance', ...
    'PDF pp. 5/7','paper-given','unknown/not-reported');
provenance.deltaC = provenance_item(seed.deltaC,'capacitance modulation', ...
    'PDF pp. 5/8','paper-given','unknown/not-reported');
provenance.Gp = provenance_item(seed.Gp,'not reported', ...
    'none','assumption','unknown/not-reported');
provenance.fmHz = provenance_item(seed.fmHz,'modulation frequency', ...
    'PDF pp. 5-8','paper-given','unknown/not-reported');
provenance.Zsource = provenance_item(50,'not reported; 50 ohm seed', ...
    'none','assumption','unknown/not-reported');
provenance.Zload = provenance_item(50,'not reported; 50 ohm seed', ...
    'none','assumption','unknown/not-reported');
if strcmp(topology,'crow')
    provenance.L0 = provenance_item(seed.L0,'parallel resonator inductor', ...
        'PDF p. 7','paper-given','unknown/not-reported');
    provenance.R0 = provenance_item(seed.R0,'not reported', ...
        'none','assumption','unknown/not-reported');
    provenance.Cblock = provenance_item(seed.Cblock,'DC-block capacitor', ...
        'PDF p. 7','paper-given','unknown/not-reported');
end

names = fieldnames(provenance);
for iname = 1:numel(names)
    name = names{iname};
    if isfield(cfg,name)
        item = provenance.(name);
        item.value = cfg.(name);
        item.source = ['physicalCfg.' name];
        item.location = 'user input';
        item.status = 'assumption';
        if isfield(cfg,'provenance') && isstruct(cfg.provenance) && ...
                isfield(cfg.provenance,name)
            supplied = cfg.provenance.(name);
            if isstruct(supplied)
                fields = fieldnames(supplied);
                for jf = 1:numel(fields)
                    item.(fields{jf}) = supplied.(fields{jf});
                end
            end
        end
        provenance.(name) = item;
    end
end
end

% -------------------------------------------------------------------------
function item = provenance_item(value,source,location,status,uncertainty)
item = struct('value',value,'source',source,'location',location, ...
    'status',status,'uncertainty',uncertainty);
end

% -------------------------------------------------------------------------
function C = bulk_capacitance(t,C0,deltaC,Omega,phase,type,duty)
if ~isnumeric(t) || ~isreal(t) || any(~isfinite(t),'all')
    error('Capacitance time input must contain finite real values.');
end
switch type
    case 'sinusoidal'
        C = C0+deltaC*cos(Omega*t+phase);
    case 'square'
        high = mod(Omega*t+phase,2*pi) < 2*pi*duty;
        C = C0+deltaC*(2*double(high)-1);
end
if any(~isfinite(C),'all') || any(C <= 0,'all')
    error('The prescribed bulk capacitance became nonpositive.');
end
end

% -------------------------------------------------------------------------
function C = finite_capacitance(cellIndex,t,C0,deltaC,Omega,phase,type,duty, ...
        phaseByCell,amplitudeScale)
if ~isnumeric(cellIndex) || isempty(cellIndex) || ~isvector(cellIndex) || ...
        any(~isfinite(cellIndex)) || any(cellIndex ~= round(cellIndex)) || ...
        any(cellIndex < 1)
    error('cellIndex must contain positive finite integers.');
end
if ~isnumeric(t) || ~isreal(t) || any(~isfinite(t),'all')
    error('Finite capacitance time input must contain finite real values.');
end
indices = cellIndex(:);
phaseValues = expand_cell_parameter(phaseByCell,indices,'phaseByCell');
amplitudeValues = expand_cell_parameter(amplitudeScale,indices, ...
    'amplitudeScaleByCell');
timeRow = t(:).';
argument = Omega*timeRow+phase+phaseValues;
switch type
    case 'sinusoidal'
        C = C0+deltaC*amplitudeValues.*cos(argument);
    case 'square'
        high = mod(argument,2*pi) < 2*pi*duty;
        C = C0+deltaC*amplitudeValues.*(2*double(high)-1);
end
if isscalar(t)
    C = C(:,1);
elseif isscalar(cellIndex)
    C = reshape(C,size(t));
end
if any(~isfinite(C),'all') || any(C <= 0,'all')
    error('The prescribed finite-chain capacitance became nonpositive.');
end
end

% -------------------------------------------------------------------------
function values = expand_cell_parameter(parameter,indices,label)
if isscalar(parameter)
    values = repmat(parameter,numel(indices),1);
elseif max(indices) <= numel(parameter)
    values = parameter(indices);
    values = values(:);
else
    error('%s does not cover all requested cell indices.',label);
end
end

% -------------------------------------------------------------------------
function [Aconstant,AinverseC,Lk] = bulk_matrices( ...
        k,a,Ls,S,Rs,Gp,topology,L0,R0,Cblock)
if ~isnumeric(k) || ~isscalar(k) || ~isreal(k) || ~isfinite(k)
    error('Bloch k must be a finite real scalar.');
end
Lk = Ls+2*S*cos(k*a);
if Lk <= 0
    error('The Bloch series inductance is nonpositive.');
end
d = 1-exp(-1i*k*a);
dc = conj(d);
if strcmp(topology,'sspp')
    Aconstant = [0,-d/Lk;0,-Rs/Lk];
    AinverseC = [-Gp,0;dc,0];
elseif isinf(Cblock)
    Aconstant = [0,-d/Lk,-1/L0; ...
                 0,-Rs/Lk,0; ...
                 0,0,-R0/L0];
    AinverseC = [-Gp,0,0;dc,0,0;1,0,0];
else
    Aconstant = [0,-d/Lk,-1/L0,0; ...
                 0,-Rs/Lk,0,0; ...
                 0,0,-R0/L0,-1/Cblock; ...
                 0,0,1/L0,0];
    AinverseC = [-Gp,0,0,0;dc,0,0,0;1,0,0,0;0,0,0,0];
end
end

% -------------------------------------------------------------------------
function A = bulk_state_matrix(k,t,a,Ls,S,Rs,Gp,topology,L0,R0,Cblock, ...
        C0,deltaC,Omega,phase,type,duty)
if ~isscalar(t) || ~isreal(t) || ~isfinite(t)
    error('bulkStateMatrix requires a finite real scalar time.');
end
[Aconstant,AinverseC] = bulk_matrices( ...
    k,a,Ls,S,Rs,Gp,topology,L0,R0,Cblock);
C = bulk_capacitance(t,C0,deltaC,Omega,phase,type,duty);
A = Aconstant+AinverseC/C;
end

% -------------------------------------------------------------------------
function O = observable_matrix(k,t,a,Ls,S,topology,L0,Cblock, ...
        C0,deltaC,Omega,phase,type,duty)
C = bulk_capacitance(t,C0,deltaC,Omega,phase,type,duty);
Lk = Ls+2*S*cos(k*a);
if strcmp(topology,'sspp')
    O = [1/C,0;0,1/Lk];
elseif isinf(Cblock)
    O = [1/C,0,0;0,1/Lk,0;0,0,1/L0];
else
    O = [1/C,0,0,0;0,1/Lk,0,0; ...
         0,0,1/L0,0;0,0,0,1/Cblock];
end
end

% -------------------------------------------------------------------------
function omegaMax = estimate_omega_max(model,capMin)
kValues = linspace(-pi/model.cell.a,pi/model.cell.a,129);
omegaMax = 0;
for k = kValues
    [Aconstant,AinverseC] = model.functions.bulkMatrices(k);
    roots = eig(Aconstant+AinverseC/capMin);
    omegaMax = max(omegaMax,max(abs(roots)));
end
if ~isfinite(omegaMax) || omegaMax <= 0
    error('Unable to determine a positive finite circuit frequency bound.');
end
end

% -------------------------------------------------------------------------
function value = read_text(cfg,name,defaultValue,allowed)
if isfield(cfg,name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
if isstring(value) && isscalar(value)
    value = char(value);
elseif ~ischar(value) || size(value,1) ~= 1
    error('physicalCfg.%s must be a text scalar.',name);
end
match = strcmpi(value,allowed);
if ~any(match)
    error('physicalCfg.%s must be one of: %s.',name,strjoin(allowed,', '));
end
value = allowed{find(match,1)};
end

% -------------------------------------------------------------------------
function value = read_logical(cfg,name,defaultValue)
if isfield(cfg,name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
if ~islogical(value) || ~isscalar(value)
    error('physicalCfg.%s must be a logical scalar.',name);
end
end

% -------------------------------------------------------------------------
function value = read_positive(cfg,name,defaultValue)
value = read_finite_real(cfg,name,defaultValue);
if value <= 0
    error('physicalCfg.%s must be positive.',name);
end
end

% -------------------------------------------------------------------------
function value = read_positive_or_inf(cfg,name,defaultValue)
if isfield(cfg,name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        isnan(value) || value <= 0
    error('physicalCfg.%s must be positive and may be Inf.',name);
end
end

% -------------------------------------------------------------------------
function value = read_nonnegative(cfg,name,defaultValue)
value = read_finite_real(cfg,name,defaultValue);
if value < 0
    error('physicalCfg.%s must be nonnegative.',name);
end
end

% -------------------------------------------------------------------------
function value = read_finite_real(cfg,name,defaultValue)
if isfield(cfg,name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value)
    error('physicalCfg.%s must be a finite real scalar.',name);
end
end

% -------------------------------------------------------------------------
function value = read_integer(cfg,name,defaultValue,minimumValue)
value = read_finite_real(cfg,name,defaultValue);
if value ~= round(value) || value < minimumValue
    error('physicalCfg.%s must be an integer not smaller than %d.', ...
        name,minimumValue);
end
end

% -------------------------------------------------------------------------
function value = read_cell_array(cfg,name,defaultValue,cellCount,allowSigned)
if isfield(cfg,name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
if ~isnumeric(value) || ~isvector(value) || ~isreal(value) || ...
        any(~isfinite(value)) || ~(isscalar(value) || numel(value) == cellCount)
    error('physicalCfg.%s must be a scalar or contain cellCount values.',name);
end
value = value(:);
if ~allowSigned && any(value < 0)
    error('physicalCfg.%s must be nonnegative.',name);
end
end
