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
% Scope: scalar, non-dispersive, lossless media in normalized units.
%   kValues are real (complex k, e.g. for damped continua, is not supported by
%   the real half-trace gap criterion). eps/mu must be finite and nonzero
%   real scalars/layers. This method applies only to space-uniform,
%   piecewise-constant, instantaneously switched eps(t) with continuous D/B
%   at each time interface.
%
% INPUTS
%   kValues   : nonempty vector of conserved real wavenumbers (row output).
%   epsLayers : row/column vector, length M >= 1 (permittivity per layer).
%   muLayers  : scalar (broadcast to M layers) or length-M vector.
%   durations : row/column vector, length M >= 1 (duration of each layer;
%               each entry finite nonnegative real, total T = sum > 0).
%
% OUTPUT (struct)
%   result.k          : 1 x nK row of wavenumbers.
%   result.U          : 2 x 2 x nK one-period D/B monodromy per wavenumber
%                       (layer 1 applied first).
%   result.omega      : 2 x nK complex quasifrequencies, rows sorted by
%                       [real(omega), imag(omega)] (omega(1,:) <= omega(2,:)).
%   result.lambda     : 2 x nK monodromy eigenvalues (same row order).
%   result.halfTrace  : 1 x nK complex trace(U)/2.
%   result.gapMask    : 1 x nK logical, true where |real(halfTrace)| > 1
%                       (with a 1e-12 safety margin, see below).
%   result.T          : total period T = sum(durations).
%   result.Omega      : temporal frequency 2*pi/T.
%   result.epsLayers  : 1 x M permittivity layers.
%   result.muLayers   : 1 x M permeability layers.
%   result.durations  : 1 x M layer durations.

%==========================================================================
% Input validation (no NaN/Inf or silently-empty pseudo-results)
%==========================================================================
if nargin < 4
    error('tmm_bands requires kValues, epsLayers, muLayers, and durations.');
end

kValues = kValues(:).';
if isempty(kValues)
    error('kValues must be a nonempty vector of real wavenumbers.');
end
if any(~isreal(kValues)) || any(~isfinite(kValues))
    error('kValues must contain only finite real wavenumbers.');
end

epsLayers = epsLayers(:).';
muLayers  = muLayers(:).';
durations = durations(:).';

M = numel(epsLayers);
if M < 1 || numel(muLayers) ~= 1 && numel(muLayers) ~= M || numel(durations) ~= M
    error(['epsLayers and durations must have equal, nonempty lengths, and ', ...
        'muLayers must be a scalar or match that length.']);
end

% eps/mu: finite real nonzero (lossless, non-singular constitutive response).
if any(~isreal(epsLayers)) || any(~isfinite(epsLayers)) || any(epsLayers == 0)
    error('epsLayers must be finite real nonzero values.');
end
if any(~isreal(muLayers)) || any(~isfinite(muLayers)) || any(muLayers == 0)
    error('muLayers must be finite real nonzero values.');
end

% Durations: finite nonnegative real, total period strictly positive.
if any(~isreal(durations)) || any(~isfinite(durations)) || any(durations < 0)
    error('durations must be finite nonnegative real values.');
end
T = sum(durations);
if ~(T > 0)
    error('The total modulation period T = sum(durations) must be positive.');
end

% Broadcast a scalar permeability to the M layers (validated above: either
% scalar or already length M), so the layer loop can index muLayers(m).
if isscalar(muLayers)
    muLayers = muLayers * ones(1, M);
end
Omega = 2*pi/T;

%==========================================================================
% Scan: single-period D/B monodromy and its Floquet eigenvalues
%==========================================================================
nK = numel(kValues);
Uall = complex(zeros(2, 2, nK));
omega = complex(zeros(2, nK));
lambda = complex(zeros(2, nK));
halfTrace = complex(zeros(1, nK));
gapMask = false(1, nK);

for ik = 1:nK
    % Single-period D/B monodromy (layer 1 applied first). Never raise U to
    % high powers: only the one-period U and its roots are used, avoiding
    % overflow in strong-gain regimes.
    U = eye(2);
    for m = 1:M
        generator = -1i*kValues(ik)*[0, 1/muLayers(m); 1/epsLayers(m), 0];
        U = expm(generator*durations(m))*U;
    end
    Uall(:, :, ik) = U;

    % Robust 2x2 eigenvalue pair: closed-form quadratic plus the reciprocal
    % relation det(U) = lambda1*lambda2. Direct eig(U) can round the small
    % root of a strongly gainful unitary monodromy to ~0; the product
    % relation recovers it from the other root with full relative precision.
    % The two names lam1/lam2 are neutral: which one carries the larger
    % magnitude flips at tr < -2 (both roots real negative), and the rows
    % are sorted below anyway.
    tr  = U(1,1) + U(2,2);
    det = U(1,1)*U(2,2) - U(1,2)*U(2,1);
    if ~isfinite(tr) || ~isfinite(det)
        error(['Monodromy at k = %.16g is not finite (unstable/overflow). ', ...
            'Reduce the scan range or layer durations.'], kValues(ik));
    end
    s = sqrt(tr*tr - 4*det);
    lam1 = 0.5*(tr + s);
    lam2 = det/lam1;              % keeps lam1*lam2 == det exactly

    w1 = 1i*log(lam1)/T;
    w2 = 1i*log(lam2)/T;
    [~, order] = sortrows([real(w1), imag(w1);
                           real(w2), imag(w2)], [1 2]);
    wPair = [w1; w2];
    lPair = [lam1; lam2];
    omega(:, ik) = wPair(order);
    lambda(:, ik) = lPair(order);

    halfTrace(ik) = tr/2;
    % Gap criterion |real(tr/2)| > 1 with a 1e-12 safety margin: a monodromy
    % computed in floating point can land epsilon above 1 right AT a passband
    % edge (where |real(halfTrace)| = 1 exactly), which would otherwise flip
    % that single column to an in-gap (spurious-gap) classification.
    gapMask(ik) = abs(real(halfTrace(ik))) > 1 + 1e-12;
end

result.k = kValues;
result.U = Uall;
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
