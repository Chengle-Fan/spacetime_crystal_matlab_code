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
required = {'Omega','Mtime','mList','Cinverse','stateDimension','modelKind', ...
    'modelSnapshot'};
missing = required(~isfield(fourier,required));
if ~isempty(missing)
    error('fourier is missing field(s): %s.',strjoin(missing,', '));
end
if fourier.stateDimension ~= model.bulk.stateDimension || ...
        ~strcmp(fourier.modelKind,model.kind)
    error('fourier and model describe different circuit state spaces.');
end
if ~isequaln(fourier.modelSnapshot,model.snapshot)
    error(['fourier and model do not share the same physical snapshot; ' ...
        'rerun tl_pwe_fourier after changing any model parameter.']);
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
        tl_normalize_eigenvectors(rightVectors,leftVectors);
    foldedRaw = fold_frequency(raw,Omega);

    rightNorm = sum(abs(rightVectors).^2,1).';
    central = sum(abs(rightVectors(centerRows,:)).^2,1).'./rightNorm;
    O0 = model.functions.observableMatrix(kScan(ik),0);
    observable = sum(abs(O0*rightVectors(centerRows,:)).^2,1).';
    observable = observable./max(max(observable),realmin);
    score = central.*(0.25+0.75*observable);
    [~,ranked] = sort(score,'descend');
    chosen = select_representatives(ranked,raw,rightVectors,nState,Omega);

    candidateRight = rightVectors(:,chosen);
    candidateFolded = foldedRaw(chosen);
    if isempty(previousLeft)
        [~,order] = sortrows([real(candidateFolded),imag(candidateFolded)],[1 2]);
    else
        overlap = abs(previousLeft'*candidateRight);
        order = tl_match_modes(overlap);
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
function folded = fold_frequency(omega,Omega)
folded = mod(real(omega)+Omega/2,Omega)-Omega/2+1i*imag(omega);
end

function chosen = select_representatives(ranked,omega,vectors,nState,Omega)
% Replicas have frequencies separated by j*Omega and the same eigenvector
% after shifting its harmonic index by j. Keep independent degenerate modes,
% but do not let two replicas of a strong resonance displace a weak branch.
nHarmonic = size(vectors,1)/nState;
chosen = zeros(0,1);
for index = ranked(:).'
    duplicate = false;
    candidate = reshape(vectors(:,index),nState,nHarmonic);
    for previous = chosen(:).'
        shift = round(real(omega(index)-omega(previous))/Omega);
        if shift == 0 || abs(shift) >= nHarmonic || ...
                abs(omega(index)-omega(previous)-shift*Omega) > 1e-5*Omega
            continue;
        end
        rows = max(1,1-shift):min(nHarmonic,nHarmonic-shift);
        reference = reshape(vectors(:,previous),nState,nHarmonic);
        a = candidate(:,rows);
        b = reference(:,rows+shift);
        overlap = abs(a(:)'*b(:))/max(norm(a(:))*norm(b(:)),realmin);
        if overlap > 1-1e-5
            duplicate = true;
            break;
        end
    end
    if ~duplicate, chosen(end+1,1) = index; end %#ok<AGROW>
    if numel(chosen) == nState, return; end
end
error('Unable to select independent Floquet representatives; increase Mtime/Nt.');
end
