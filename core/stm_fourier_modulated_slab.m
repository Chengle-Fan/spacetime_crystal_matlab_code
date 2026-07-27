function coeff = stm_fourier_modulated_slab(m, n, p)
%STM_FOURIER_MODULATED_SLAB Analytic spacetime Fourier coefficient for a modulated slab.
%
%   coeff = stm_fourier_modulated_slab(m, n, p)
%
% Convention:
%   eps(x,t) = sum_{m,n} eps_{m,n} exp(i*n*g*x - i*m*Omega*t)
%
% Inputs:
%   m  -- temporal harmonic index
%   n  -- spatial harmonic index
%   p  -- parameter struct from stm_preset_modulated_slab()
%
% Because the modulation is purely sinusoidal, only m = 0, +/-1 are
% non-zero. The integral splits into static and modulated spatial regions.
%
% See also STM_PRESET_MODULATED_SLAB, STM_PERMITTIVITY_MODULATED_SLAB, STM_INTERVAL_FOURIER

if nargin < 3 || isempty(p)
    p = stm_preset_modulated_slab();
end

staticPart = stm_interval_fourier(n, 0, p.xModStart, p.Lambda);
modPart = stm_interval_fourier(n, p.xModStart, p.xModEnd, p.Lambda);

switch m
    case 0                             % DC temporal component
        coeff = p.eps1*staticPart + p.epsc*modPart;
    case 1                             % fundamental (e^{-i*Omega*t})
        coeff = 1i*(p.epsc*p.modDepth/2)*modPart;
    case -1                            % fundamental (e^{+i*Omega*t})
        coeff = -1i*(p.epsc*p.modDepth/2)*modPart;
    otherwise
        coeff = 0;                     % sinusoidal modulation: only m = 0, +/-1
end
end
