function p = stm_preset_modulated_slab()
%STM_PRESET_MODULATED_SLAB Parameters for a unit cell with a partially modulated slab.
%
%   p = stm_preset_modulated_slab();
%
% Returns a parameter struct for a 1D spacetime crystal where part of the
% unit cell is static and the remaining part is sinusoidally modulated in
% time. This is a common reference configuration for ST-PWE studies.
%
% Unit-cell layout (normalized Lambda = c0 = 1):
%   [0, xModStart)        — static dielectric, eps = eps1
%   [xModStart, Lambda)   — temporally modulated, eps = epsc*(1+modDepth*sin(Omega*t))
%
% Normalization:
%   kBar = k*Lambda/(2*pi)
%   fBar = omega*Lambda/(2*pi*c0)
%
% See also STM_PERMITTIVITY_MODULATED_SLAB, STM_FOURIER_MODULATED_SLAB, STM_INTERVAL_FOURIER

p.Lambda = 1;           % normalized spatial period
p.c0 = 1;               % normalized speed of light
p.g = 2*pi/p.Lambda;    % spatial reciprocal lattice vector

p.eps1 = 2;             % static region permittivity
p.epsc = 6;             % modulated region background permittivity
p.modDepth = 0.6;       % modulation depth (fraction of epsc)
p.xModStart = 3*p.Lambda/4;   % start of modulated region
p.xModEnd = p.Lambda;         % end of modulated region

p.OmegaBar = 0.20;      % normalized modulation frequency
p.Omega = 2*pi*p.c0/p.Lambda * p.OmegaBar;  % physical modulation angular frequency
p.T = 2*pi/p.Omega;     % modulation period
end
