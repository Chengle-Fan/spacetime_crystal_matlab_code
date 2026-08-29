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
% INPUTS
%   epsFun : function handle epsFun(X,TIME) -> eps_r on an Nt-by-Nx array.
%   muFun  : optional function handle muFun(X,TIME); default
%            @(x,t) ones(size(x)) (vacuum mu_r = 1).
%   pweCfg : struct with required fields
%              Lambda, T, Nspace, Mtime
%            and optional fields
%              Nx (default max(128, 8*(2*Nspace+1)))
%              Nt (default max(128, 8*(2*Mtime+1))).
%
% OUTPUT (struct)
%   fourier.Lambda/.T/.g/.Omega : geometry and reciprocal frequencies.
%   fourier.Nspace/.Mtime/.Nx/.Nt : truncation and sample counts.
%   fourier.mOrders : 1 x (4*Mtime+1) row of temporal orders.
%   fourier.nOrders : 1 x (4*Nspace+1) row of spatial orders.
%   fourier.epsilon : (4*Mtime+1) x (4*Nspace+1) eps coefficients; entry
%                     (im,in) <-> m = mOrders(im), n = nOrders(in).
%   fourier.mu      : same layout for the mu coefficients.
%   fourier.epsCoeff : @(m,n) scalar lookup with hard range validation.
%   fourier.muCoeff  : @(m,n) scalar lookup with hard range validation.

% Default permeability (vacuum mu_r = 1).
if nargin < 2 || isempty(muFun)
    muFun = @(x,t) ones(size(x));
end

Lambda = pweCfg.Lambda;
T       = pweCfg.T;
Nspace  = pweCfg.Nspace;
Mtime   = pweCfg.Mtime;

% Sample counts (oversample 8x the kept order by default).
if isfield(pweCfg, 'Nx') && ~isempty(pweCfg.Nx)
    Nx = pweCfg.Nx;
else
    Nx = max(128, 8*(2*Nspace+1));
end
if isfield(pweCfg, 'Nt') && ~isempty(pweCfg.Nt)
    Nt = pweCfg.Nt;
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
epsilon = sample_table(values, mOrders, nOrders, g, Omega, X, Time);

% Permeability.
mvalues = muFun(X, Time);
if ~isequal(size(mvalues), size(X))
    error('muFun must return an Nt-by-Nx array matching the sample grid.');
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
if ~isscalar(m) || ~isscalar(n)
    error('Fourier coefficient lookup requires scalar (m,n).');
end
if abs(m) > 2*Mtime || abs(n) > 2*Nspace
    error('Fourier coefficient (m=%d,n=%d) outside sampled difference range [-2M,2M]x[-2N,2N]', m, n);
end
im = m + 2*Mtime + 1;
in = n + 2*Nspace + 1;
c = coeffMat(im, in);
end
