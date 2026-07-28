function mode = temporal_domain_wall_mode(kValues, ...
    epsLeft, muLeft, durationsLeft, ...
    epsRight, muRight, durationsRight, nLeft, nRight)
%TEMPORAL_DOMAIN_WALL_MODE Match growing/decaying Floquet states in time.
%
% The left temporal crystal occupies NLEFT periods before the temporal
% domain wall and the right crystal occupies NRIGHT periods after it.
% At a common momentum gap, a state localized at the wall grows toward the
% wall under the left monodromy and decays away from it under the right
% monodromy. The routine scans conserved k and minimizes the determinant
% mismatch between those two one-dimensional Floquet eigenspaces.
%
% This establishes a bulk-interface mode of the specified transfer
% matrices. It does not by itself calculate or prove a quantized temporal
% Zak phase; that additionally requires a symmetry and unit-cell analysis.

if abs(sum(durationsLeft)-sum(durationsRight)) > ...
        1e-12*max(1,sum(durationsLeft))
    error('The two temporal crystals must have the same period.');
end
if any([nLeft,nRight] < 1) || ...
        any([nLeft,nRight] ~= round([nLeft,nRight]))
    error('nLeft and nRight must be positive integers.');
end

kValues = kValues(:).';
Nk = numel(kValues);
mismatch = inf(1,Nk);
growthMagnitude = nan(1,Nk);
decayMagnitude = nan(1,Nk);
leftHalfTrace = complex(nan(1,Nk));
rightHalfTrace = complex(nan(1,Nk));
gapTolerance = 1e-7;

for ik = 1:Nk
    UL = temporal_crystal_monodromy(kValues(ik), ...
        epsLeft,muLeft,durationsLeft);
    UR = temporal_crystal_monodromy(kValues(ik), ...
        epsRight,muRight,durationsRight);
    [VL,lambdaL] = eig(UL,'vector');
    [VR,lambdaR] = eig(UR,'vector');
    [growthMagnitude(ik),idGrowth] = max(abs(lambdaL));
    [decayMagnitude(ik),idDecay] = min(abs(lambdaR));
    leftHalfTrace(ik) = trace(UL)/2;
    rightHalfTrace(ik) = trace(UR)/2;

    if growthMagnitude(ik) <= 1+gapTolerance || ...
            decayMagnitude(ik) >= 1-gapTolerance
        continue;
    end
    vGrowth = VL(:,idGrowth);
    vDecay = VR(:,idDecay);
    vGrowth = vGrowth/norm(vGrowth);
    vDecay = vDecay/norm(vDecay);
    mismatch(ik) = abs(det([vGrowth,vDecay]));
end

[bestMismatch,bestId] = min(mismatch);
if ~isfinite(bestMismatch)
    error(['No common momentum gap was found. Expand kValues or check ', ...
        'the temporal-crystal parameters.']);
end

kBest = kValues(bestId);
if bestId > 1 && bestId < Nk && ...
        isfinite(mismatch(bestId-1)) && isfinite(mismatch(bestId+1))
    tolerance = 1e-11*max(1,abs(kBest));
    options = optimset('TolX',tolerance,'Display','off');
    [kRefined,mismatchRefined] = fminbnd( ...
        @(kq) wall_mismatch(kq, ...
        epsLeft,muLeft,durationsLeft, ...
        epsRight,muRight,durationsRight,gapTolerance), ...
        kValues(bestId-1),kValues(bestId+1),options);
    if mismatchRefined < bestMismatch
        kBest = kRefined;
        bestMismatch = mismatchRefined;
    end
end

UL = temporal_crystal_monodromy(kBest, ...
    epsLeft,muLeft,durationsLeft);
UR = temporal_crystal_monodromy(kBest, ...
    epsRight,muRight,durationsRight);
[VL,lambdaL] = eig(UL,'vector');
[VR,lambdaR] = eig(UR,'vector');
[~,idGrowth] = max(abs(lambdaL));
[~,idDecay] = min(abs(lambdaR));
vGrowth = VL(:,idGrowth)/norm(VL(:,idGrowth));
vDecay = VR(:,idDecay)/norm(VR(:,idDecay));

cellIndex = -nLeft:nRight;
states = complex(zeros(2,numel(cellIndex)));
states(:,1) = vGrowth;
for q = 1:nLeft
    states(:,q+1) = UL*states(:,q);
end
interfaceId = nLeft+1;
for q = 1:nRight
    states(:,interfaceId+q) = UR*states(:,interfaceId+q-1);
end
interfaceScale = norm(states(:,interfaceId));
states = states/interfaceScale;

mode.kGrid = kValues;
mode.mismatch = mismatch;
mode.growthMagnitude = growthMagnitude;
mode.decayMagnitude = decayMagnitude;
mode.leftHalfTrace = leftHalfTrace;
mode.rightHalfTrace = rightHalfTrace;
mode.k = kBest;
mode.bestId = bestId;
mode.bestMismatch = bestMismatch;
mode.leftMultiplier = lambdaL(idGrowth);
mode.rightMultiplier = lambdaR(idDecay);
mode.leftEigenvector = vGrowth;
mode.rightEigenvector = vDecay;
mode.cellIndex = cellIndex;
mode.states = states;
mode.stateNorm = sqrt(sum(abs(states).^2,1));
mode.interfaceId = interfaceId;
mode.leftMonodromy = UL;
mode.rightMonodromy = UR;
mode.T = sum(durationsLeft);
end

function mismatch = wall_mismatch(k, ...
    epsLeft,muLeft,durationsLeft, ...
    epsRight,muRight,durationsRight,gapTolerance)
UL = temporal_crystal_monodromy(k,epsLeft,muLeft,durationsLeft);
UR = temporal_crystal_monodromy(k,epsRight,muRight,durationsRight);
[VL,lambdaL] = eig(UL,'vector');
[VR,lambdaR] = eig(UR,'vector');
[growthMagnitude,idGrowth] = max(abs(lambdaL));
[decayMagnitude,idDecay] = min(abs(lambdaR));
if growthMagnitude <= 1+gapTolerance || ...
        decayMagnitude >= 1-gapTolerance
    mismatch = inf;
    return;
end
vGrowth = VL(:,idGrowth)/norm(VL(:,idGrowth));
vDecay = VR(:,idDecay)/norm(VR(:,idDecay));
mismatch = abs(det([vGrowth,vDecay]));
end
