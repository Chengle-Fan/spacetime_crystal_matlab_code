function [TM, details] = temporal_multilayer_tmm(omegaInitial, ...
    epsInitial, muInitial, epsSlabs, muSlabs, durations, ...
    epsFinal, muFinal)
%TEMPORAL_MULTILAYER_TMM Transfer matrix of an arbitrary temporal stack.
%
% Fields are column vectors. Therefore the chronological map is
%   TM = MM_final*DM_M*...*MM_2*DM_1*MM_1.
%
% The printed product in Eq. (11) of Ramaccia et al. should be interpreted
% with its propagation ordering; the explicit order above removes any
% ambiguity caused by left- versus right-acting conventions.

epsSlabs = epsSlabs(:).';
durations = durations(:).';
M = numel(epsSlabs);

if M == 0
    [TM,tau,rho] = temporal_interface_matrix( ...
        epsInitial,muInitial,epsFinal,muFinal);
    details.tau = tau;
    details.rho = rho;
    details.nInitial = sqrt(epsInitial*muInitial);
    details.nSlabs = [];
    details.omegaSlabs = [];
    details.matchingMatrices = {TM};
    details.delayMatrices = {};
    return;
end

if isscalar(muSlabs)
    muSlabs = repmat(muSlabs, 1, M);
else
    muSlabs = muSlabs(:).';
end
if numel(muSlabs) ~= M || numel(durations) ~= M
    error('epsSlabs, muSlabs, and durations must have equal lengths.');
end

nInitial = sqrt(epsInitial*muInitial);
nSlabs = sqrt(epsSlabs.*muSlabs);

[MMfirst, tau(1), rho(1)] = temporal_interface_matrix( ...
    epsInitial, muInitial, epsSlabs(1), muSlabs(1));
TM = MMfirst;
matchingMatrices = cell(1,M+1);
delayMatrices = cell(1,M);
matchingMatrices{1} = MMfirst;

for m = 1:M
    DM = temporal_delay_matrix(omegaInitial, nInitial, ...
        nSlabs(m), durations(m));
    delayMatrices{m} = DM;

    if m < M
        epsNext = epsSlabs(m+1);
        muNext = muSlabs(m+1);
    else
        epsNext = epsFinal;
        muNext = muFinal;
    end

    [MMnext, tau(m+1), rho(m+1)] = temporal_interface_matrix( ...
        epsSlabs(m), muSlabs(m), epsNext, muNext);
    matchingMatrices{m+1} = MMnext;
    TM = MMnext*DM*TM;
end

details.tau = tau;
details.rho = rho;
details.nInitial = nInitial;
details.nSlabs = nSlabs;
details.omegaSlabs = (nInitial./nSlabs)*omegaInitial;
details.matchingMatrices = matchingMatrices;
details.delayMatrices = delayMatrices;
end
