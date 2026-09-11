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

temporalSlices = 512;
if ~strcmp(model.modulation.type,'square')
    temporalSlices = tl_option('integer',tmmCfg,'temporalSlices',512,4);
end
[segmentMidpoints,segmentDurations,integrationName] = ...
    tl_time_layers(model,temporalSlices);

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
        tl_normalize_eigenvectors(rightVectors,leftVectors);
    raw = 1i*log(values)/T;
    folded = fold_frequency(raw,Omega);

    if isempty(previousLeft)
        [~,order] = sortrows([real(folded),imag(folded)],[1 2]);
    else
        overlap = abs(previousLeft'*rightVectors);
        order = tl_match_modes(overlap);
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
function folded = fold_frequency(omega,Omega)
folded = mod(real(omega)+Omega/2,Omega)-Omega/2+1i*imag(omega);
end
