function result = tmm_bands(kValues, epsLayers, muLayers, durations)
%TMM_BANDS Floquet quasifrequencies of a piecewise-constant temporal crystal.
%
% Continuous [D;B] state, single-period monodromy, Floquet quasifrequencies.
%
% For a plane wave exp(i*k*x) propagating through layers that are constant in
% space but switch in time, the D/B state evolves by
%     d/dt [D;B] = -i*k*[0 1/mu; 1/eps 0]*[D;B]
% within each layer. D and B are continuous at every time interface, so the
% one-period evolution is a product of layer exponentials (no separate
% interface matrix is needed). If U is the one-period monodromy with
% eigenvalues lambda, the Floquet quasifrequencies are
%     omega = i*log(lambda)/T,
% whose principal logarithm places Re(omega) in the first temporal Brillouin
% zone. The gap condition is |trace(U)/2| > 1, so that omega acquires a
% nonzero imaginary part.
%
% INPUTS
%   kValues   : numeric vector of conserved wavenumbers -> treated as a row.
%   epsLayers : row/column vector, length M >= 1 (permittivity per layer).
%   muLayers  : scalar (broadcast to M layers) or length-M vector.
%   durations : row/column vector, length M (duration of each layer).
%
% OUTPUT (struct)
%   result.k          : 1 x nK row of wavenumbers.
%   result.omega      : 2 x nK complex quasifrequencies, rows sorted by
%                       [real(omega), imag(omega)] (omega(1,:) <= omega(2,:)).
%   result.lambda     : 2 x nK monodromy eigenvalues (same row order).
%   result.halfTrace  : 1 x nK complex trace(U)/2.
%   result.gapMask    : 1 x nK logical, true where |real(halfTrace)| > 1.
%   result.T          : total period T = sum(durations).
%   result.Omega      : temporal frequency 2*pi/T.
%   result.epsLayers  : 1 x M permittivity layers.
%   result.muLayers   : 1 x M permeability layers.
%   result.durations  : 1 x M layer durations.

kValues = kValues(:).';

epsLayers = epsLayers(:).';
durations = durations(:).';
M = numel(epsLayers);

if isscalar(muLayers)
    muLayers = repmat(muLayers, 1, M);
else
    muLayers = muLayers(:).';
end
if numel(muLayers) ~= M || numel(durations) ~= M
    error('epsLayers, muLayers, and durations must have equal lengths.');
end

T = sum(durations);
Omega = 2*pi/T;

nK = numel(kValues);
omega = complex(zeros(2, nK));
lambda = complex(zeros(2, nK));
halfTrace = complex(zeros(1, nK));
gapMask = false(1, nK);

for ik = 1:nK
    % Single-period D/B monodromy (layer 1 applied first). Never raise U to
    % high powers: only the one-period U and its eig are used, avoiding
    % overflow in strong-gain regimes.
    U = eye(2);
    for m = 1:M
        generator = -1i*kValues(ik)*[0, 1/muLayers(m); 1/epsLayers(m), 0];
        U = expm(generator*durations(m))*U;
    end

    lam = eig(U);
    w = 1i*log(lam)/T;
    [~, order] = sortrows([real(w), imag(w)], [1 2]);
    omega(:, ik) = w(order);
    lambda(:, ik) = lam(order);

    halfTrace(ik) = trace(U)/2;
    gapMask(ik) = abs(real(halfTrace(ik))) > 1 + 1e-12;
end

result.k = kValues;
result.omega = omega;
result.lambda = lambda;
result.halfTrace = halfTrace;
result.gapMask = gapMask;
result.T = T;
result.Omega = Omega;
result.epsLayers = epsLayers;
result.muLayers = muLayers;
result.durations = durations;
end
