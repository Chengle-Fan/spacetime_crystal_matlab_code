function fourier = pwe_fourier(epsFun, muFun, pweCfg)
%PWE_FOURIER Fourier coefficients of a 1D spacetime material sampled on the
%endpoint-free periodic cell, covering the full convolution-difference range.
%
% The material is evaluated on the endpoint-free grid
%     x = (0:Nx-1)*Lambda/Nx,   t = (0:Nt-1)*T/Nt,
% and expanded in the convention
%     eps(x,t) = sum_{m,n} eps_mn exp(i*n*g*x - i*m*Omega*t),
%     g = 2*pi/Lambda,  Omega = 2*pi/T.
% The returned tables span the temporal orders m = -2*Mtime:2*Mtime and the
% spatial orders n = -2*Nspace:2*Nspace, exactly the difference range needed
% by the convolution matrices assembled in pwe_bands. A lookup outside this
% range raises an error instead of silently zero-padding.
%
% Anti-aliasing contract: the sample counts must satisfy
%     Nx >= 4*Nspace + 1   and   Nt >= 4*Mtime + 1,
% so that the requested difference orders (-2N..2N, -2M..2M) never alias onto
% each other on the discrete grid. This lower bound does NOT stop true
% Nyquist-external harmonics of the material (e.g. the sharp edges of a
% binary profile) from aliasing in: keep the default oversampling and confirm
% convergence by doubling Nx/Nt. As a regression, a constant medium must give
% a single nonzero coefficient at (m,n) = (0,0) and (numerically) zero
% coefficients at every other requested order.
%
% INPUTS
%   epsFun : function handle epsFun(X,TIME) -> eps_r on an Nt-by-Nx array
%            (must be finite everywhere).
%   muFun  : optional function handle muFun(X,TIME); default
%            @(x,t) ones(size(x)) (vacuum mu_r = 1).
%   pweCfg : struct with required fields
%              Lambda   (positive finite real scalar: spatial period)
%              T        (positive finite real scalar: temporal period)
%              Nspace   (nonnegative finite integer: spatial truncation)
%              Mtime    (nonnegative finite integer: temporal truncation)
%            and optional fields
%              Nx (default max(128, 8*(2*Nspace+1))), must satisfy
%                 Nx >= 4*Nspace+1
%              Nt (default max(128, 8*(2*Mtime+1))), must satisfy
%                 Nt >= 4*Mtime+1.
%
% OUTPUT (struct)
%   fourier.Lambda/.T/.g/.Omega : geometry and reciprocal frequencies.
%   fourier.Nspace/.Mtime/.Nx/.Nt : truncation and sample counts.
%   fourier.mOrders : 1 x (4*Mtime+1) row of temporal orders.
%   fourier.nOrders : 1 x (4*Nspace+1) row of spatial orders.
%   fourier.epsilon : (4*Mtime+1) x (4*Nspace+1) eps coefficients; entry
%                     (im,in) <-> m = mOrders(im), n = nOrders(in).
%   fourier.mu      : same layout for the mu coefficients.
%   fourier.epsCoeff : @(m,n) scalar lookup (m,n finite integers within the
%                     sampled difference range; hard validation).
%   fourier.muCoeff  : @(m,n) scalar lookup (same validation).

%==========================================================================
% Validate configuration (periods, truncation, sample counts)
%==========================================================================
if nargin < 1 || isempty(epsFun)
    error('pwe_fourier requires an epsFun function handle.');
end
if ~isa(epsFun, 'function_handle')
    error('epsFun must be a function handle.');
end
if nargin < 2 || isempty(muFun)
    muFun = @(x,t) ones(size(x));
end
if ~isa(muFun, 'function_handle')
    error('muFun must be a function handle.');
end

validate_pwe_periods(pweCfg, 'Lambda');
validate_pwe_periods(pweCfg, 'T');
Lambda = pweCfg.Lambda;
T       = pweCfg.T;
Nspace  = validate_pwe_truncation(pweCfg, 'Nspace');
Mtime   = validate_pwe_truncation(pweCfg, 'Mtime');

% Sample counts (oversample 8x the kept order by default). The lower bounds
% 4*Nspace+1 / 4*Mtime+1 prevent the requested difference orders from
% colliding on the discrete grid (P1-01).
minNx = 4*Nspace + 1;
minNt = 4*Mtime + 1;
if isfield(pweCfg, 'Nx') && ~isempty(pweCfg.Nx)
    Nx = pweCfg.Nx;
    if ~isscalar(Nx) || ~isreal(Nx) || ~isfinite(Nx) || Nx < 1 || ...
            Nx ~= round(Nx)
        error('pweCfg.Nx must be a positive finite integer scalar.');
    end
    if Nx < minNx
        error(['pweCfg.Nx = %d is too small for Nspace = %d: need ', ...
            'Nx >= 4*Nspace+1 = %d to avoid aliasing between the requested ', ...
            'difference orders.'], Nx, Nspace, minNx);
    end
else
    Nx = max(128, 8*(2*Nspace+1));
end
if isfield(pweCfg, 'Nt') && ~isempty(pweCfg.Nt)
    Nt = pweCfg.Nt;
    if ~isscalar(Nt) || ~isreal(Nt) || ~isfinite(Nt) || Nt < 1 || ...
            Nt ~= round(Nt)
        error('pweCfg.Nt must be a positive finite integer scalar.');
    end
    if Nt < minNt
        error(['pweCfg.Nt = %d is too small for Mtime = %d: need ', ...
            'Nt >= 4*Mtime+1 = %d to avoid aliasing between the requested ', ...
            'difference orders.'], Nt, Mtime, minNt);
    end
else
    Nt = max(128, 8*(2*Mtime+1));
end

g     = 2*pi/Lambda;
Omega = 2*pi/T;

% Endpoint-free periodic grid (no repeated endpoint at x=Lambda, t=T).
x = (0:Nx-1)*Lambda/Nx;
t = (0:Nt-1)*T/Nt;
[X, Time] = meshgrid(x, t);            % Nt-by-Nx

% Difference-order ranges covered by the convolution matrices.
mOrders = -2*Mtime:2*Mtime;            % row
nOrders = -2*Nspace:2*Nspace;          % row

% Permittivity.
values = epsFun(X, Time);
if ~isequal(size(values), size(X))
    error('epsFun must return an Nt-by-Nx array matching the sample grid.');
end
if ~all(isfinite(values), 'all')
    error('epsFun returned non-finite values on the sample grid.');
end
epsilon = sample_table(values, mOrders, nOrders, g, Omega, X, Time);

% Permeability.
mvalues = muFun(X, Time);
if ~isequal(size(mvalues), size(X))
    error('muFun must return an Nt-by-Nx array matching the sample grid.');
end
if ~all(isfinite(mvalues), 'all')
    error('muFun returned non-finite values on the sample grid.');
end
mu = sample_table(mvalues, mOrders, nOrders, g, Omega, X, Time);

fourier.Lambda   = Lambda;
fourier.T        = T;
fourier.g        = g;
fourier.Omega    = Omega;
fourier.Nspace   = Nspace;
fourier.Mtime    = Mtime;
fourier.Nx       = Nx;
fourier.Nt       = Nt;
fourier.mOrders  = mOrders;
fourier.nOrders  = nOrders;
fourier.epsilon  = epsilon;
fourier.mu       = mu;
fourier.epsCoeff = @(m,n) coeff_lookup(epsilon, Nspace, Mtime, m, n);
fourier.muCoeff  = @(m,n) coeff_lookup(mu,      Nspace, Mtime, m, n);
end

% -------------------------------------------------------------------------
function validate_pwe_periods(pweCfg, fieldName)
%VALIDATE_PWE_PERIODS A period (Lambda or T) must be a positive finite real.
if ~isfield(pweCfg, fieldName)
    error('pweCfg.%s is required.', fieldName);
end
val = pweCfg.(fieldName);
if ~isscalar(val) || ~isreal(val) || ~isfinite(val) || val <= 0
    error('pweCfg.%s must be a finite positive real scalar.', fieldName);
end
end

% -------------------------------------------------------------------------
function ord = validate_pwe_truncation(pweCfg, fieldName)
%VALIDATE_PWE_TRUNCATION A truncation (Nspace/Mtime) must be a nonneg integer.
if ~isfield(pweCfg, fieldName)
    error('pweCfg.%s is required.', fieldName);
end
ord = pweCfg.(fieldName);
if ~isscalar(ord) || ~isreal(ord) || ~isfinite(ord) || ord < 0 || ...
        ord ~= round(ord)
    error('pweCfg.%s must be a nonnegative finite integer scalar.', fieldName);
end
end

% -------------------------------------------------------------------------
function coeff = sample_table(values, mOrders, nOrders, g, Omega, X, Time)
%SAMPLE_TABLE Evaluate p_{mn} = <p, exp(i n g x - i m Omega t)> on the
%endpoint-free grid for every (m,n) in the order lists.
coeff = complex(zeros(numel(mOrders), numel(nOrders)));
for im = 1:numel(mOrders)
    m = mOrders(im);
    for in = 1:numel(nOrders)
        n = nOrders(in);
        coeff(im,in) = mean(values .* exp(-1i*n*g*X + 1i*m*Omega*Time), 'all');
    end
end
end

% -------------------------------------------------------------------------
function c = coeff_lookup(coeffMat, Nspace, Mtime, m, n)
%COEFF_LOOKUP Hard-range coefficient lookup. Never silently zero-pads.
if ~isscalar(m) || ~isscalar(n) || ~isreal(m) || ~isreal(n) || ...
        ~isfinite(m) || ~isfinite(n) || m ~= round(m) || n ~= round(n)
    error('Fourier coefficient lookup requires finite integer scalar (m,n).');
end
if abs(m) > 2*Mtime || abs(n) > 2*Nspace
    error('Fourier coefficient (m=%d,n=%d) outside sampled difference range [-2M,2M]x[-2N,2N]', m, n);
end
im = m + 2*Mtime + 1;
in = n + 2*Nspace + 1;
c = coeffMat(im, in);
end
