function fourier = pwe_fourier(epsFun, pweCfg)
%PWE_FOURIER Fourier convolution matrix of a photonic time crystal.
%
% The material is spatially uniform, periodic in time, lossless and
% nondispersive, with mu_r = 1. The convention is
%
%   epsilon(t) = sum_m epsilon_m exp(-i*m*Omega*t),  Omega = 2*pi/T.
%
% epsFun is sampled on the endpoint-free grid t=(0:Nt-1)*T/Nt. To prevent
% the Fourier orders required by the convolution matrix from colliding, Nt
% must satisfy Nt >= 4*Mtime+1. A square wave should use a much larger Nt and
% be checked for sampling convergence.
%
% INPUTS
%   epsFun : function handle epsFun(t), returning an array the same size as t.
%            Values must be finite, real and strictly positive.
%   pweCfg : scalar struct with
%              T      : positive temporal period
%              Mtime  : nonnegative temporal-harmonic truncation
%              Nt     : positive integer sample count, Nt >= 4*Mtime+1
%
% OUTPUT
%   fourier.T, fourier.Omega, fourier.Mtime
%   fourier.mList    : harmonic orders -Mtime:Mtime (column)
%   fourier.Cepsilon : permittivity convolution matrix

if nargin ~= 2 || ~isa(epsFun,'function_handle') || ...
        ~isstruct(pweCfg) || ~isscalar(pweCfg)
    error('Use pwe_fourier(epsFun,pweCfg) with a scalar configuration struct.');
end

T = positive_scalar(pweCfg,'T');
Mtime = nonnegative_integer(pweCfg,'Mtime');
Nt = positive_integer(pweCfg,'Nt');
minNt = 4*Mtime + 1;
if Nt < minNt
    error('pweCfg.Nt must satisfy Nt >= 4*Mtime+1 = %d.',minNt);
end

Omega = 2*pi/T;
t = (0:Nt-1)*T/Nt;
epsilon = epsFun(t);
if ~isnumeric(epsilon) || ~isequal(size(epsilon),size(t)) || ...
        ~isreal(epsilon) || any(~isfinite(epsilon)) || any(epsilon <= 0)
    error(['epsFun(t) must return finite, real, strictly positive values ', ...
        'with exactly the same size as t.']);
end

% The eigenproblem retains m=-Mtime:Mtime, so convolution differences require
% coefficients over -2*Mtime:2*Mtime.
mDifference = -2*Mtime:2*Mtime;
epsilonCoeff = complex(zeros(size(mDifference)));
for im = 1:numel(mDifference)
    m = mDifference(im);
    epsilonCoeff(im) = mean(epsilon .* exp(1i*m*Omega*t));
end

mList = (-Mtime:Mtime).';
dm = mList - mList.';
Cepsilon = epsilonCoeff(dm + 2*Mtime + 1);

% A positive real material produces a Hermitian positive-definite convolution
% matrix. These checks catch Fourier-sign and indexing regressions early.
scale = max(1,norm(Cepsilon,'fro'));
if norm(Cepsilon-Cepsilon','fro')/scale > 1e-11
    error('The permittivity convolution matrix is not Hermitian.');
end
Cepsilon = (Cepsilon+Cepsilon')/2;
[~,cholFlag] = chol(Cepsilon);
if cholFlag ~= 0
    error(['The permittivity convolution matrix is not positive definite; ', ...
        'increase Nt or inspect epsilon(t).']);
end

fourier.T = T;
fourier.Omega = Omega;
fourier.Mtime = Mtime;
fourier.mList = mList;
fourier.Cepsilon = Cepsilon;
end

% -------------------------------------------------------------------------
function value = positive_scalar(cfg,name)
if ~isfield(cfg,name)
    error('pweCfg.%s is required.',name);
end
value = cfg.(name);
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value <= 0
    error('pweCfg.%s must be a positive finite real scalar.',name);
end
end

% -------------------------------------------------------------------------
function value = nonnegative_integer(cfg,name)
if ~isfield(cfg,name)
    error('pweCfg.%s is required.',name);
end
value = cfg.(name);
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value < 0 || value ~= round(value)
    error('pweCfg.%s must be a nonnegative finite integer.',name);
end
end

% -------------------------------------------------------------------------
function value = positive_integer(cfg,name)
value = nonnegative_integer(cfg,name);
if value < 1
    error('pweCfg.%s must be a positive integer.',name);
end
end
