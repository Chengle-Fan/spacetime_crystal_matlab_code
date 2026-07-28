function p = stm_fig2_parameters()
%STM_FIG2_PARAMETERS Normalized parameters used in Park and Min, Fig. 2.
%
% Normalization: Lambda = c0 = 1. Hence
%   kBar = k*Lambda/(2*pi)
%   fBar = omega*Lambda/(2*pi*c0)
% and the modulation period is T = 1/OmegaBar.

p.Lambda = 1;
p.c0 = 1;
p.g = 2*pi/p.Lambda;

p.eps1 = 2;
p.epsc = 6;
p.modDepth = 0.6;
p.xModStart = 3*p.Lambda/4;
p.xModEnd = p.Lambda;

p.OmegaBar = 0.20;
p.Omega = 2*pi*p.c0/p.Lambda*p.OmegaBar;
p.T = 2*pi/p.Omega;
end
