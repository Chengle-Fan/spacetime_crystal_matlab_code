function result = tl_tmm_bands(model,kScan,tmmCfg)
%TL_TMM_BANDS Temporal monodromy bands of an SSPP or CROW unit cell.
%
% This is a time-domain Floquet transfer matrix at fixed Bloch k, not a
% spatial ABCD cascade.  The conservative Q/Phi state is continuous across
% piecewise-constant capacitance interfaces.

if nargin ~= 3 || ~isstruct(model) || ~isscalar(model) || ...
        ~isstruct(tmmCfg) || ~isscalar(tmmCfg)
    error('Use tl_tmm_bands(model,kScan,tmmCfg).');
end
if ~isnumeric(kScan) || isempty(kScan) || ~isvector(kScan) || ...
        ~isreal(kScan) || any(~isfinite(kScan))
    error('kScan must be a nonempty finite real vector.');
end
kScan = kScan(:).';
nState = model.bulk.stateDimension;
nK = numel(kScan);
T = model.modulation.period;
Omega = model.modulation.OmegaRadPerSec;

if strcmp(model.modulation.type,'square')
    [segmentMidpoints,segmentDurations] = square_segments(model);
    integrationName = 'exact-square-layers';
else
    temporalSlices = read_integer(tmmCfg,'temporalSlices',512,4);
    dt = T/temporalSlices;
    segmentMidpoints = ((1:temporalSlices)-0.5)*dt;
    segmentDurations = dt*ones(1,temporalSlices);
    integrationName = 'midpoint-expm';
end

Uall = complex(zeros(nState,nState,nK));
lambda = complex(zeros(nState,nK));
omegaRaw = complex(zeros(nState,nK));
omegaFolded = complex(zeros(nState,nK));
rightAll = complex(zeros(nState,nState,nK));
leftAll = complex(zeros(nState,nState,nK));
conditionNumber = zeros(nState,nK);
selectionConfidence = ones(nState,nK);
selectedIndices = zeros(nState,nK);
previousLeft = [];

for ik = 1:nK
    U = eye(nState);
    for segment = 1:numel(segmentDurations)
        A = model.functions.bulkStateMatrix( ...
            kScan(ik),segmentMidpoints(segment));
        U = expm(A*segmentDurations(segment))*U;
    end
    if any(~isfinite(U),'all')
        error('The monodromy matrix became non-finite at k=%.16g.',kScan(ik));
    end
    [rightVectors,lambdaMatrix,leftVectors] = eig(U);
    values = diag(lambdaMatrix);
    if any(values == 0) || any(~isfinite(values))
        error('The Floquet multiplier is zero or non-finite at k=%.16g.', ...
            kScan(ik));
    end
    [rightVectors,leftVectors,condition] = ...
        normalize_eigenvectors(rightVectors,leftVectors);
    raw = 1i*log(values)/T;
    folded = fold_frequency(raw,Omega);

    if isempty(previousLeft)
        [~,order] = sortrows([real(folded),imag(folded)],[1 2]);
    else
        overlap = abs(previousLeft'*rightVectors);
        order = best_permutation(overlap);
        matched = overlap(sub2ind([nState nState],(1:nState).',order(:)));
        selectionConfidence(:,ik) = matched;
    end
    rightVectors = rightVectors(:,order);
    leftVectors = leftVectors(:,order);
    values = values(order);
    raw = raw(order);
    folded = folded(order);
    condition = condition(order);

    Uall(:,:,ik) = U;
    lambda(:,ik) = values;
    omegaRaw(:,ik) = raw;
    omegaFolded(:,ik) = folded;
    rightAll(:,:,ik) = rightVectors;
    leftAll(:,:,ik) = leftVectors;
    conditionNumber(:,ik) = condition;
    selectedIndices(:,ik) = order;
    previousLeft = leftVectors;
end

result.k = kScan;
result.kaOverPi = kScan*model.cell.a/pi;
result.U = Uall;
result.lambda = lambda;
result.omegaRaw = omegaRaw;
result.omegaFolded = omegaFolded;
result.omegaSelected = omegaFolded;
result.selectedRawIndices = selectedIndices;
result.replicaOrder = zeros(size(omegaFolded));
result.selectionConfidence = selectionConfidence;
result.conditionNumber = conditionNumber;
result.rightEigenvectors = rightAll;
result.leftEigenvectors = leftAll;
result.integrationName = integrationName;
result.segmentMidpoints = segmentMidpoints;
result.segmentDurations = segmentDurations;
result.Omega = Omega;
result.modelSnapshot = model.snapshot;
result.tmmConfig = tmmCfg;
end

% -------------------------------------------------------------------------
function [midpoints,durations] = square_segments(model)
T = model.modulation.period;
Omega = model.modulation.OmegaRadPerSec;
phase = model.modulation.phase;
duty = model.modulation.dutyCycle;
phaseBoundaries = [0,2*pi*duty];
times = [0,T];
for boundary = phaseBoundaries
    for order = -3:3
        value = (boundary-phase+2*pi*order)/Omega;
        if value > 64*eps(T) && value < T-64*eps(T)
            times(end+1) = value; %#ok<AGROW>
        end
    end
end
times = unique(sort(times));
durations = diff(times);
midpoints = times(1:end-1)+durations/2;
if any(durations <= 0) || abs(sum(durations)-T) > 128*eps(T)
    error('Unable to construct exact square-wave time layers.');
end
end

% -------------------------------------------------------------------------
function [V,W,conditionNumber] = normalize_eigenvectors(V,W)
nMode = size(V,2);
conditionNumber = zeros(nMode,1);
for mode = 1:nMode
    nv = norm(V(:,mode));
    nw = norm(W(:,mode));
    if nv == 0 || nw == 0 || ~isfinite(nv) || ~isfinite(nw)
        error('TMM returned an invalid left or right eigenvector.');
    end
    V(:,mode) = V(:,mode)/nv;
    W(:,mode) = W(:,mode)/nw;
    conditionNumber(mode) = ...
        1/max(abs(W(:,mode)'*V(:,mode)),realmin);
end
end

% -------------------------------------------------------------------------
function order = best_permutation(overlap)
n = size(overlap,1);
candidate = perms(1:n);
scores = zeros(size(candidate,1),1);
rows = (1:n).';
for index = 1:size(candidate,1)
    cols = candidate(index,:).';
    scores(index) = sum(overlap(sub2ind([n n],rows,cols)));
end
[~,best] = max(scores);
order = candidate(best,:);
end

% -------------------------------------------------------------------------
function folded = fold_frequency(omega,Omega)
folded = mod(real(omega)+Omega/2,Omega)-Omega/2+1i*imag(omega);
end

% -------------------------------------------------------------------------
function value = read_integer(cfg,name,defaultValue,minimumValue)
if isfield(cfg,name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value ~= round(value) || value < minimumValue
    error('tmmCfg.%s must be an integer not smaller than %d.', ...
        name,minimumValue);
end
end
