function U = temporal_crystal_monodromy(k, epsSequence, ...
    muSequence, durations)
%TEMPORAL_CRYSTAL_MONODROMY One-period D/B evolution at conserved k.
%
% Normalized units eps0=mu0=c0=1 are used. For a plane wave exp(i*k*x),
%   d/dt [D;B] = -i*k*[0 1/mu; 1/eps 0]*[D;B].
%
% D and B are continuous at every time interface, so no separate interface
% matrix is required in this state representation.

epsSequence = epsSequence(:).';
durations = durations(:).';
M = numel(epsSequence);
if isscalar(muSequence)
    muSequence = repmat(muSequence,1,M);
else
    muSequence = muSequence(:).';
end
if numel(muSequence) ~= M || numel(durations) ~= M
    error('epsSequence, muSequence, and durations must have equal lengths.');
end

U = eye(2);
for m = 1:M
    generator = -1i*k*[0, 1/muSequence(m); ...
        1/epsSequence(m), 0];
    U = expm(generator*durations(m))*U;
end
end
