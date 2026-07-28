function response = temporal_finite_crystal_response( ...
    kValues, epsSequence, muSequence, durations, nPeriods, ...
    epsBackground, muBackground, inputDirectional)
%TEMPORAL_FINITE_CRYSTAL_RESPONSE Exact response of a finite time crystal.
%
% A spatially uniform crystal is active for NPERIODS and is embedded
% between identical static temporal backgrounds. The input and output are
% the directional amplitudes [E_plus;E_minus] in that background. D and B
% are continuous when the modulation is turned on and off.

if nargin < 8 || isempty(inputDirectional)
    inputDirectional = [1;0];
end
if size(inputDirectional,1) ~= 2 || size(inputDirectional,2) ~= 1
    error('inputDirectional must be a 2-by-1 column vector.');
end
if ~isscalar(nPeriods) || nPeriods < 0 || nPeriods ~= round(nPeriods)
    error('nPeriods must be a nonnegative integer.');
end

kValues = kValues(:).';
Nk = numel(kValues);
outputDirectional = complex(zeros(2,Nk));
finalState = complex(zeros(2,Nk));
periodMatrix = complex(zeros(2,2,Nk));
multipliers = complex(zeros(2,Nk));
initialState = temporal_directional_to_db( ...
    inputDirectional,epsBackground,muBackground);

for ik = 1:Nk
    U = temporal_crystal_monodromy(kValues(ik),epsSequence, ...
        muSequence,durations);
    state = (U^nPeriods)*initialState;
    outputDirectional(:,ik) = temporal_db_to_directional( ...
        state,epsBackground,muBackground);
    finalState(:,ik) = state;
    periodMatrix(:,:,ik) = U;
    multipliers(:,ik) = eig(U);
end

response.k = kValues;
response.inputDirectional = inputDirectional;
response.outputDirectional = outputDirectional;
response.forward = outputDirectional(1,:);
response.backward = outputDirectional(2,:);
response.relativePhase = angle(response.backward./response.forward);
response.finalState = finalState;
response.periodMatrix = periodMatrix;
response.multipliers = multipliers;
response.nPeriods = nPeriods;
response.T = sum(durations);
end
