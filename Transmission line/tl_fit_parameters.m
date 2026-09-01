function fitResult = tl_fit_parameters(basePhysicalCfg,fitData,fitCfg)
%TL_FIT_PARAMETERS Fit identifiable TL parameters with Base MATLAB.
%
%   fitResult = tl_fit_parameters(basePhysicalCfg,fitData,fitCfg)
%
% fitCfg.parameterNames is a cell array selected from a, Ls, mutualS, Rs,
% C0, deltaC, Gp, Cpar, fmHz, Zsource, Zload, L0, R0, Cblock, and
% modulationPhase.  Supported data blocks are:
%
%   fitData.staticDispersion.kRadPerM, frequencyHz [,bandIndex,uncertaintyHz]
%   fitData.imagOmega.kRadPerM, imagOmegaRadPerSec
%       [,bandIndex,uncertaintyRadPerSec]
%   fitData.impedance.valueOhm [,uncertaintyOhm]
%   fitData.anchors.boundaryFrequencyHz [,boundaryUncertaintyHz]
%   fitData.anchors.fcolHz [,fcolUncertaintyHz]       (CROW only)
%
% Positive parameters use logarithmic coordinates.  deltaC/C0 and
% mutualS/Ls use bounded transforms.  The routine reports a local numerical
% Jacobian; rank deficiency means the supplied data do not identify all
% requested parameters and should not be hidden by a best-fit number.

if nargin ~= 3 || ~isstruct(basePhysicalCfg) || ...
        ~isscalar(basePhysicalCfg) || ~isstruct(fitData) || ...
        ~isscalar(fitData) || ~isstruct(fitCfg) || ~isscalar(fitCfg)
    error(['Use tl_fit_parameters(basePhysicalCfg,fitData,fitCfg) with ' ...
        'scalar structs.']);
end
if ~isfield(fitCfg,'parameterNames') || isempty(fitCfg.parameterNames)
    error('fitCfg.parameterNames is required.');
end
parameterNames = fitCfg.parameterNames;
if isstring(parameterNames), parameterNames = cellstr(parameterNames); end
if ~iscellstr(parameterNames) || numel(unique(parameterNames)) ~= ...
        numel(parameterNames)
    error('fitCfg.parameterNames must contain unique parameter names.');
end
allowed = {'a','Ls','mutualS','Rs','C0','deltaC','Gp','Cpar', ...
    'fmHz','Zsource','Zload','L0','R0','Cblock','modulationPhase'};
if any(~ismember(parameterNames,allowed))
    error('fitCfg.parameterNames contains an unsupported parameter.');
end

baseModel = tl_build_model(basePhysicalCfg);
if strcmp(baseModel.kind,'sspp') && ...
        any(ismember(parameterNames,{'L0','R0','Cblock'}))
    error('SSPP data cannot fit CROW resonator parameters.');
end
data = prepare_data(fitData,fitCfg,baseModel);
theta0 = encode_parameters(baseModel,parameterNames);
maxIterations = read_integer(fitCfg,'maxIterations',400,1);
maxEvaluations = read_integer(fitCfg,'maxEvaluations',2000,1);
displayMode = read_text(fitCfg,'display','off',{'off','iter','final','notify'});
options = optimset('Display',displayMode,'MaxIter',maxIterations, ...
    'MaxFunEvals',maxEvaluations,'TolX',1e-9,'TolFun',1e-9);
objective = @(theta) objective_value(theta,basePhysicalCfg, ...
    parameterNames,data);
[theta,objectiveValue,exitFlag,output] = fminsearch( ...
    objective,theta0,options);
[bestCfg,parameterValues] = decode_parameters( ...
    theta,basePhysicalCfg,baseModel,parameterNames);
residual = residual_vector(bestCfg,data);

[jacobian,singularValues,rankValue,covarianceTheta,standardTheta] = ...
    local_uncertainty(theta,basePhysicalCfg,baseModel,parameterNames,data);
confidence95 = struct();
for index = 1:numel(parameterNames)
    lowerTheta = theta;
    upperTheta = theta;
    lowerTheta(index) = lowerTheta(index)-1.96*standardTheta(index);
    upperTheta(index) = upperTheta(index)+1.96*standardTheta(index);
    [~,lowerValues] = decode_parameters(lowerTheta,basePhysicalCfg, ...
        baseModel,parameterNames);
    [~,upperValues] = decode_parameters(upperTheta,basePhysicalCfg, ...
        baseModel,parameterNames);
    interval = sort([lowerValues.(parameterNames{index}), ...
        upperValues.(parameterNames{index})]);
    confidence95.(parameterNames{index}) = interval;
    if ~isfield(bestCfg,'provenance') || ~isstruct(bestCfg.provenance)
        bestCfg.provenance = struct();
    end
    bestCfg.provenance.(parameterNames{index}) = struct( ...
        'status','fit','source','tl_fit_parameters', ...
        'location','fitData','uncertainty',interval);
end
bestModel = tl_build_model(bestCfg);

fitResult.parameterNames = parameterNames;
fitResult.parameterValues = parameterValues;
fitResult.confidence95 = confidence95;
fitResult.model = bestModel;
fitResult.physicalConfig = bestCfg;
fitResult.residual = residual;
fitResult.weightedRms = sqrt(mean(residual.^2));
fitResult.objectiveValue = objectiveValue;
fitResult.exitFlag = exitFlag;
fitResult.output = output;
fitResult.theta = theta;
fitResult.jacobian = jacobian;
fitResult.covarianceTheta = covarianceTheta;
fitResult.jacobianSingularValues = singularValues;
fitResult.identifiableRank = rankValue;
fitResult.parameterCount = numel(parameterNames);
fitResult.isLocallyIdentifiable = rankValue == numel(parameterNames);
fitResult.dataSummary = data.summary;
fitResult.warning = ['Confidence intervals are local linearized estimates. ' ...
    'Use held-out measurements and profile/Monte-Carlo analysis before ' ...
    'calling any value an experimental component parameter.'];
end

% -------------------------------------------------------------------------
function data = prepare_data(input,cfg,model)
data.blocks = cell(1,5);
blockCount = 0;
summary = struct('staticDispersionCount',0,'imagOmegaCount',0, ...
    'impedanceCount',0,'anchorCount',0);
if isfield(input,'staticDispersion') && ~isempty(input.staticDispersion)
    block = input.staticDispersion;
    block.k = read_vector_field(block,'kRadPerM');
    block.value = read_vector_field(block,'frequencyHz');
    assert_same_size(block.k,block.value,'static dispersion');
    block.band = read_band(block,numel(block.k));
    defaultSigma = read_positive(cfg,'frequencyUncertaintyHz', ...
        max(median(abs(block.value))*0.01,1));
    block.sigma = read_uncertainty(block,'uncertaintyHz', ...
        numel(block.k),defaultSigma);
    block.type = 'frequency';
    blockCount = blockCount+1;
    data.blocks{blockCount} = block;
    summary.staticDispersionCount = numel(block.k);
end
if isfield(input,'imagOmega') && ~isempty(input.imagOmega)
    block = input.imagOmega;
    block.k = read_vector_field(block,'kRadPerM');
    block.value = read_vector_field(block,'imagOmegaRadPerSec');
    assert_same_size(block.k,block.value,'imaginary frequency');
    block.band = read_band(block,numel(block.k));
    defaultSigma = read_positive(cfg,'imagOmegaUncertaintyRadPerSec', ...
        max(median(abs(block.value))*0.05,1));
    block.sigma = read_uncertainty(block,'uncertaintyRadPerSec', ...
        numel(block.k),defaultSigma);
    block.type = 'imagOmega';
    blockCount = blockCount+1;
    data.blocks{blockCount} = block;
    summary.imagOmegaCount = numel(block.k);
end
if isfield(input,'impedance') && ~isempty(input.impedance)
    block = input.impedance;
    block.value = read_vector_field(block,'valueOhm');
    defaultSigma = read_positive(cfg,'impedanceUncertaintyOhm', ...
        max(median(abs(block.value))*0.02,1e-3));
    block.sigma = read_uncertainty(block,'uncertaintyOhm', ...
        numel(block.value),defaultSigma);
    block.type = 'impedance';
    blockCount = blockCount+1;
    data.blocks{blockCount} = block;
    summary.impedanceCount = numel(block.value);
end
if isfield(input,'anchors') && ~isempty(input.anchors)
    anchors = input.anchors;
    if isfield(anchors,'boundaryFrequencyHz')
        block.type = 'boundary';
        block.value = validate_scalar(anchors.boundaryFrequencyHz, ...
            'boundaryFrequencyHz');
        if isfield(anchors,'boundaryUncertaintyHz')
            block.sigma = validate_positive_scalar( ...
                anchors.boundaryUncertaintyHz,'boundaryUncertaintyHz');
        else
            block.sigma = max(block.value*0.01,1);
        end
        blockCount = blockCount+1;
        data.blocks{blockCount} = block;
        summary.anchorCount = summary.anchorCount+1;
    end
    if isfield(anchors,'fcolHz')
        if strcmp(model.kind,'sspp')
            error('fitData.anchors.fcolHz is only meaningful for CROW.');
        end
        block.type = 'fcol';
        block.value = validate_scalar(anchors.fcolHz,'fcolHz');
        if isfield(anchors,'fcolUncertaintyHz')
            block.sigma = validate_positive_scalar( ...
                anchors.fcolUncertaintyHz,'fcolUncertaintyHz');
        else
            block.sigma = max(block.value*0.01,1);
        end
        blockCount = blockCount+1;
        data.blocks{blockCount} = block;
        summary.anchorCount = summary.anchorCount+1;
    end
end
if blockCount == 0
    error('fitData contains no supported measurement block.');
end
data.blocks = data.blocks(1:blockCount);
data.summary = summary;
end

% -------------------------------------------------------------------------
function value = objective_value(theta,baseCfg,names,data)
try
    seedModel = tl_build_model(baseCfg);
    [cfg,~] = decode_parameters(theta,baseCfg,seedModel,names);
    residual = residual_vector(cfg,data);
    if any(~isfinite(residual))
        value = 1e100;
    else
        value = sum(residual.^2);
    end
catch
    value = 1e100;
end
end

% -------------------------------------------------------------------------
function residual = residual_vector(cfg,data)
model = tl_build_model(cfg);
residual = zeros(0,1);
for index = 1:numel(data.blocks)
    block = data.blocks{index};
    switch block.type
        case {'frequency','imagOmega'}
            prediction = zeros(size(block.k));
            for sample = 1:numel(block.k)
                omega = static_mode(model,block.k(sample),block.band(sample));
                if strcmp(block.type,'frequency')
                    prediction(sample) = real(omega)/(2*pi);
                else
                    prediction(sample) = imag(omega);
                end
            end
        case 'impedance'
            prediction = model.derived.Zreference*ones(size(block.value));
        case 'boundary'
            prediction = model.derived.provisionalBoundaryFrequencyHz;
        case 'fcol'
            prediction = model.derived.fcolHz;
    end
    residual = [residual;(prediction(:)-block.value(:))./block.sigma(:)]; %#ok<AGROW>
end
end

% -------------------------------------------------------------------------
function omega = static_mode(model,k,band)
[Aconstant,AinverseC] = model.functions.bulkMatrices(k);
values = 1i*eig(Aconstant+AinverseC/model.shunt.totalC0);
scale = max(abs(values));
positive = values(real(values) > -1e-10*max(scale,1));
[~,order] = sort(real(positive),'ascend');
positive = positive(order);
if numel(positive) >= 2 && abs(real(positive(1))) < 1e-10*max(scale,1)
    positive = positive(2:end);
end
if band > numel(positive)
    error('Requested static bandIndex %d does not exist at k=%.16g.',band,k);
end
omega = positive(band);
end

% -------------------------------------------------------------------------
function theta = encode_parameters(model,names)
theta = zeros(numel(names),1);
for index = 1:numel(names)
    name = names{index};
    value = model_value(model,name);
    switch name
        case 'mutualS'
            ratio = min(max(value/(0.499*model.series.Ls),-0.999),0.999);
            theta(index) = atanh(ratio);
        case 'deltaC'
            ratio = min(max(value/model.shunt.C0,1e-9),1-1e-9);
            theta(index) = log(ratio/(1-ratio));
        case 'modulationPhase'
            theta(index) = value;
        case {'Rs','Gp','Cpar','R0'}
            theta(index) = log(max(value,parameter_floor(name)));
        otherwise
            if isinf(value)
                error('An infinite %s cannot be used as a fit initial value.',name);
            end
            theta(index) = log(value);
    end
end
end

% -------------------------------------------------------------------------
function [cfg,values] = decode_parameters(theta,baseCfg,baseModel,names)
cfg = baseCfg;
values = struct();
for index = 1:numel(names)
    name = names{index};
    if ismember(name,{'mutualS','deltaC'})
        continue;
    elseif strcmp(name,'modulationPhase')
        value = theta(index);
    else
        value = exp(theta(index));
    end
    cfg.(name) = value;
    values.(name) = value;
end
if isfield(cfg,'Ls'), Ls = cfg.Ls; else, Ls = baseModel.series.Ls; end
if isfield(cfg,'C0'), C0 = cfg.C0; else, C0 = baseModel.shunt.C0; end
for index = 1:numel(names)
    name = names{index};
    if strcmp(name,'mutualS')
        value = 0.499*Ls*tanh(theta(index));
        cfg.mutualS = value;
        values.mutualS = value;
    elseif strcmp(name,'deltaC')
        ratio = 1/(1+exp(-theta(index)));
        value = C0*ratio;
        cfg.deltaC = value;
        values.deltaC = value;
    end
end
end

% -------------------------------------------------------------------------
function [J,singularValues,rankValue,covariance,standard] = ...
        local_uncertainty(theta,baseCfg,baseModel,names,data)
baseResidual = residual_vector(decode_cfg(theta,baseCfg,baseModel,names),data);
J = zeros(numel(baseResidual),numel(theta));
for index = 1:numel(theta)
    step = 1e-5*(1+abs(theta(index)));
    upper = theta;
    lower = theta;
    upper(index) = upper(index)+step;
    lower(index) = lower(index)-step;
    rUpper = residual_vector(decode_cfg(upper,baseCfg,baseModel,names),data);
    rLower = residual_vector(decode_cfg(lower,baseCfg,baseModel,names),data);
    J(:,index) = (rUpper-rLower)/(2*step);
end
singularValues = svd(J);
if isempty(singularValues)
    rankValue = 0;
else
    rankValue = sum(singularValues > max(size(J))*eps(max(singularValues)));
end
dof = max(numel(baseResidual)-numel(theta),1);
variance = sum(baseResidual.^2)/dof;
covariance = pinv(J'*J)*variance;
standard = sqrt(max(real(diag(covariance)),0));
if rankValue < numel(theta)
    standard(:) = Inf;
end
end

function cfg = decode_cfg(theta,baseCfg,baseModel,names)
[cfg,~] = decode_parameters(theta,baseCfg,baseModel,names);
end

% -------------------------------------------------------------------------
function value = model_value(model,name)
switch name
    case 'a', value = model.cell.a;
    case {'Ls','Rs','mutualS'}, value = model.series.(name);
    case {'C0','deltaC','Gp','Cpar'}, value = model.shunt.(name);
    case {'L0','R0','Cblock'}, value = model.resonator.(name);
    case {'Zsource','Zload'}, value = model.ports.(name);
    case 'fmHz', value = model.modulation.fmHz;
    case 'modulationPhase', value = model.modulation.phase;
end
end

function value = parameter_floor(name)
switch name
    case {'Rs','R0'}, value = 1e-6;
    case 'Gp', value = 1e-12;
    case 'Cpar', value = 1e-16;
end
end

% -------------------------------------------------------------------------
function vector = read_vector_field(block,name)
if ~isstruct(block) || ~isscalar(block) || ~isfield(block,name)
    error('Measurement block requires field %s.',name);
end
vector = block.(name);
if ~isnumeric(vector) || ~isvector(vector) || isempty(vector) || ...
        ~isreal(vector) || any(~isfinite(vector))
    error('%s must be a nonempty finite real vector.',name);
end
vector = vector(:);
end

function band = read_band(block,count)
if isfield(block,'bandIndex') && ~isempty(block.bandIndex)
    band = block.bandIndex;
else
    band = ones(count,1);
end
if ~isnumeric(band) || ~isvector(band) || ...
        ~(isscalar(band) || numel(band) == count) || ...
        any(~isfinite(band)) || any(band ~= round(band)) || any(band < 1)
    error('bandIndex must contain positive integers.');
end
if isscalar(band), band = repmat(band,count,1); else, band = band(:); end
end

function sigma = read_uncertainty(block,name,count,defaultValue)
if isfield(block,name) && ~isempty(block.(name))
    sigma = block.(name);
else
    sigma = defaultValue;
end
if ~isnumeric(sigma) || ~isvector(sigma) || ...
        ~(isscalar(sigma) || numel(sigma) == count) || ...
        ~isreal(sigma) || any(~isfinite(sigma)) || any(sigma <= 0)
    error('%s must contain positive finite values.',name);
end
if isscalar(sigma), sigma = repmat(sigma,count,1); else, sigma = sigma(:); end
end

function assert_same_size(a,b,label)
if numel(a) ~= numel(b), error('%s vectors must have equal lengths.',label); end
end

function value = validate_scalar(value,label)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value)
    error('%s must be a finite real scalar.',label);
end
end

function value = validate_positive_scalar(value,label)
value = validate_scalar(value,label);
if value <= 0, error('%s must be positive.',label); end
end

function value = read_positive(cfg,name,defaultValue)
if isfield(cfg,name) && ~isempty(cfg.(name)), value = cfg.(name); else, value = defaultValue; end
if isempty(value) || ~isnumeric(value) || ~isscalar(value) || ...
        ~isreal(value) || ~isfinite(value) || value <= 0
    error('fitCfg.%s must be a positive finite real scalar.',name);
end
end

function value = read_integer(cfg,name,defaultValue,minimum)
if isfield(cfg,name) && ~isempty(cfg.(name)), value = cfg.(name); else, value = defaultValue; end
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value ~= round(value) || value < minimum
    error('fitCfg.%s must be an integer not smaller than %d.',name,minimum);
end
end

function value = read_text(cfg,name,defaultValue,allowed)
if isfield(cfg,name) && ~isempty(cfg.(name)), value = cfg.(name); else, value = defaultValue; end
if isstring(value) && isscalar(value), value = char(value); end
if ~ischar(value) || size(value,1) ~= 1 || ~any(strcmp(value,allowed))
    error('fitCfg.%s has an unsupported text value.',name);
end
end
