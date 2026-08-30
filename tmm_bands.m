function result = tmm_bands(kScan, epsLayers, muLayers, durations)
%TMM_BANDS Complex quasifrequency bands of a square-wave time crystal.
%
% The two temporal layers are spatially uniform, lossless and nondispersive.
% At each instantaneous time interface, D and B are continuous. Within one
% layer the state follows
%
%   d/dt [D;B] = -i*k*[0 1/mu; 1/eps 0]*[D;B].
%
% The eigenvalues lambda of the one-period 2-by-2 evolution matrix give
% omega=i*log(lambda)/T. The principal logarithm places Re(omega) in the
% first temporal Floquet zone.
%
% INPUTS
%   kScan     : nonempty vector of finite real wavenumbers
%   epsLayers : [eps1 eps2], finite and strictly positive
%   muLayers  : positive scalar or [mu1 mu2]
%   durations : [dt1 dt2], finite and strictly positive
%
% OUTPUT
%   result.k     : 1 x nK scan vector
%   result.omega : 2 x nK complex quasifrequencies

if nargin ~= 4
    error('Use tmm_bands(kScan,epsLayers,muLayers,durations).');
end
if ~isnumeric(kScan) || isempty(kScan) || ~isvector(kScan) || ...
        ~isreal(kScan) || any(~isfinite(kScan))
    error('kScan must be a nonempty vector of finite real wavenumbers.');
end
if ~isnumeric(epsLayers) || ~isvector(epsLayers) || numel(epsLayers) ~= 2 || ...
        ~isreal(epsLayers) || any(~isfinite(epsLayers)) || any(epsLayers <= 0)
    error('epsLayers must contain two finite, real, strictly positive values.');
end
if ~isnumeric(muLayers) || ~isvector(muLayers) || ...
        ~(isscalar(muLayers) || numel(muLayers) == 2) || ...
        ~isreal(muLayers) || any(~isfinite(muLayers)) || any(muLayers <= 0)
    error('muLayers must be one or two finite, real, strictly positive values.');
end
if ~isnumeric(durations) || ~isvector(durations) || numel(durations) ~= 2 || ...
        ~isreal(durations) || any(~isfinite(durations)) || any(durations <= 0)
    error('durations must contain two finite, real, strictly positive values.');
end

kScan = kScan(:).';
epsLayers = epsLayers(:).';
durations = durations(:).';
if isscalar(muLayers)
    muLayers = [muLayers muLayers];
else
    muLayers = muLayers(:).';
end
T = sum(durations);
omega = complex(zeros(2,numel(kScan)));

for ik = 1:numel(kScan)
    U = eye(2);
    for layer = 1:2
        epsr = epsLayers(layer);
        mur = muLayers(layer);
        theta = kScan(ik)*durations(layer)/sqrt(epsr*mur);
        c = cos(theta);
        s = sin(theta);

        % Exact layer propagator for the continuous [D;B] state.
        P = [c, -1i*s*sqrt(epsr/mur); ...
             -1i*s*sqrt(mur/epsr), c];
        U = P*U;
    end

    halfTrace = (U(1,1)+U(2,2))/2;
    determinant = U(1,1)*U(2,2)-U(1,2)*U(2,1);
    root = sqrt(halfTrace^2-determinant);
    candidates = [halfTrace+root, halfTrace-root];
    [~,largeIndex] = max(abs(candidates));
    lambdaLarge = candidates(largeIndex);
    if ~isfinite(lambdaLarge) || lambdaLarge == 0 || ...
            ~isfinite(determinant) || determinant == 0
        error('The one-period evolution is singular or non-finite at k=%.16g.', ...
            kScan(ik));
    end
    lambdaSmall = determinant/lambdaLarge;

    pair = 1i*log([lambdaLarge;lambdaSmall])/T;
    [~,order] = sortrows([real(pair),imag(pair)],[1 2]);
    omega(:,ik) = pair(order);
end

result.k = kScan;
result.omega = omega;
end
