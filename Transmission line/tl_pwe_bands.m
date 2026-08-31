function result = tl_pwe_bands(fourier,model,kScan,pweCfg)
%TL_PWE_BANDS Complex Floquet bands of a periodic transmission-line cell.
%
% The state convention is exp(i*k*x-i*omega*t).  All temporal replicas are
% retained in omegaRaw.  A physical representative set with one mode per
% conservative circuit state is selected by central-harmonic participation,
% then tracked across k using left/right biorthogonal overlap.

if nargin ~= 4 || ~isstruct(fourier) || ~isscalar(fourier) || ...
        ~isstruct(model) || ~isscalar(model) || ...
        ~isstruct(pweCfg) || ~isscalar(pweCfg)
    error('Use tl_pwe_bands(fourier,model,kScan,pweCfg).');
end
if ~isnumeric(kScan) || isempty(kScan) || ~isvector(kScan) || ...
        ~isreal(kScan) || any(~isfinite(kScan))
    error('kScan must be a nonempty finite real vector.');
end
kScan = kScan(:).';
required = {'Omega','Mtime','mList','Cinverse','stateDimension','modelKind'};
missing = required(~isfield(fourier,required));
if ~isempty(missing)
    error('fourier is missing field(s): %s.',strjoin(missing,', '));
end
if fourier.stateDimension ~= model.bulk.stateDimension || ...
        ~strcmp(fourier.modelKind,model.kind)
    error('fourier and model describe different circuit state spaces.');
end

Mtime = fourier.Mtime;
mList = fourier.mList;
nHarmonic = numel(mList);
nState = model.bulk.stateDimension;
nFull = nState*nHarmonic;
nK = numel(kScan);
Omega = fourier.Omega;
Ih = eye(nHarmonic);
Is = eye(nState);
W = diag(mList*Omega);
centerHarmonic = Mtime+1;
centerRows = (centerHarmonic-1)*nState+(1:nState);

omegaRaw = complex(zeros(nFull,nK));
omegaFoldedRaw = complex(zeros(nFull,nK));
centralWeight = zeros(nFull,nK);
observableWeight = zeros(nFull,nK);
omegaSelected = complex(zeros(nState,nK));
omegaSelectedRaw = complex(zeros(nState,nK));
selectedRawIndices = zeros(nState,nK);
replicaOrder = zeros(nState,nK);
selectionConfidence = zeros(nState,nK);
conditionNumberSelected = zeros(nState,nK);
rightSelected = complex(zeros(nFull,nState,nK));
leftSelected = complex(zeros(nFull,nState,nK));

previousLeft = [];
for ik = 1:nK
    [Aconstant,AinverseC] = model.functions.bulkMatrices(kScan(ik));
    H = 1i*kron(Ih,Aconstant)+1i*kron(fourier.Cinverse,AinverseC) ...
        -kron(W,Is);
    [rightVectors,eigenvalueMatrix,leftVectors] = eig(H);
    raw = diag(eigenvalueMatrix);
    [rightVectors,leftVectors,conditionNumber] = ...
        normalize_eigenvectors(rightVectors,leftVectors);
    foldedRaw = fold_frequency(raw,Omega);

    rightNorm = sum(abs(rightVectors).^2,1).';
    central = sum(abs(rightVectors(centerRows,:)).^2,1).'./rightNorm;
    O0 = model.functions.observableMatrix(kScan(ik),0);
    observable = sum(abs(O0*rightVectors(centerRows,:)).^2,1).';
    observable = observable./max(max(observable),realmin);
    score = central.*(0.25+0.75*observable);
    [~,ranked] = sort(score,'descend');
    chosen = ranked(1:nState);

    candidateRight = rightVectors(:,chosen);
    candidateFolded = foldedRaw(chosen);
    if isempty(previousLeft)
        [~,order] = sortrows([real(candidateFolded),imag(candidateFolded)],[1 2]);
    else
        overlap = abs(previousLeft'*candidateRight);
        order = best_permutation(overlap);
    end
    chosen = chosen(order);
    candidateRight = rightVectors(:,chosen);
    candidateLeft = leftVectors(:,chosen);
    candidateFolded = foldedRaw(chosen);
    candidateRaw = raw(chosen);

    omegaRaw(:,ik) = raw;
    omegaFoldedRaw(:,ik) = foldedRaw;
    centralWeight(:,ik) = central;
    observableWeight(:,ik) = observable;
    omegaSelected(:,ik) = candidateFolded;
    omegaSelectedRaw(:,ik) = candidateRaw;
    selectedRawIndices(:,ik) = chosen;
    replicaOrder(:,ik) = round((real(candidateRaw)-real(candidateFolded))/Omega);
    selectionConfidence(:,ik) = score(chosen);
    conditionNumberSelected(:,ik) = conditionNumber(chosen);
    rightSelected(:,:,ik) = candidateRight;
    leftSelected(:,:,ik) = candidateLeft;
    previousLeft = candidateLeft;
end

result.k = kScan;
result.kaOverPi = kScan*model.cell.a/pi;
result.omegaRaw = omegaRaw;
result.omegaFoldedRaw = omegaFoldedRaw;
result.omegaSelectedRaw = omegaSelectedRaw;
result.omegaSelected = omegaSelected;
result.omegaFolded = omegaSelected;
result.selectedRawIndices = selectedRawIndices;
result.replicaOrder = replicaOrder;
result.selectionConfidence = selectionConfidence;
result.centralHarmonicWeight = centralWeight;
result.observableWeight = observableWeight;
result.conditionNumber = conditionNumberSelected;
result.rightEigenvectors = rightSelected;
result.leftEigenvectors = leftSelected;
result.stateDimension = nState;
result.harmonicOrders = mList;
result.Omega = Omega;
result.modelSnapshot = model.snapshot;
result.pweConfig = pweCfg;
end

% -------------------------------------------------------------------------
function [V,W,conditionNumber] = normalize_eigenvectors(V,W)
nMode = size(V,2);
conditionNumber = zeros(nMode,1);
for mode = 1:nMode
    nv = norm(V(:,mode));
    nw = norm(W(:,mode));
    if nv == 0 || nw == 0 || ~isfinite(nv) || ~isfinite(nw)
        error('PWE returned an invalid left or right eigenvector.');
    end
    V(:,mode) = V(:,mode)/nv;
    W(:,mode) = W(:,mode)/nw;
    product = W(:,mode)'*V(:,mode);
    conditionNumber(mode) = 1/max(abs(product),realmin);
end
end

% -------------------------------------------------------------------------
function order = best_permutation(overlap)
n = size(overlap,1);
if size(overlap,2) ~= n
    error('Mode-overlap matrix must be square.');
end
if n > 8
    [~,order] = max(overlap,[],2);
    if numel(unique(order)) ~= n
        [~,order] = sort(max(overlap,[],1),'descend');
    end
    order = order(:).';
    return;
end
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
