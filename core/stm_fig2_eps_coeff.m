function coeff = stm_fig2_eps_coeff(m, n, p)
%STM_FIG2_EPS_COEFF Analytic space-time Fourier coefficient of epsilon.
%
% Convention:
%   eps(x,t) = sum_{m,n} eps_{m,n} exp(i*n*g*x - i*m*Omega*t).

if nargin < 3 || isempty(p)
    p = stm_fig2_parameters();
end

staticPart = stm_interval_fourier(n, 0, p.xModStart, p.Lambda);
modPart = stm_interval_fourier(n, p.xModStart, p.xModEnd, p.Lambda);

switch m
    case 0
        coeff = p.eps1*staticPart + p.epsc*modPart;
    case 1
        coeff = 1i*(p.epsc*p.modDepth/2)*modPart;
    case -1
        coeff = -1i*(p.epsc*p.modDepth/2)*modPart;
    otherwise
        coeff = 0;
end
end
