function epsr = stm_permittivity_modulated_slab(x, t, p)
%STM_PERMITTIVITY_MODULATED_SLAB Real-space permittivity of a partially modulated unit cell.
%
%   epsr = stm_permittivity_modulated_slab(x, t, p)
%
% Inputs:
%   x  -- spatial coordinate (scalar or compatible array)
%   t  -- time coordinate (scalar or compatible array)
%   p  -- parameter struct from stm_preset_modulated_slab()
%
% Output:
%   epsr -- relative permittivity at each (x,t), same size as x and t
%
% The unit cell is piecewise:
%   [0, xModStart)       — eps1 (static)
%   [xModStart, Lambda)  — epsc * (1 + modDepth * sin(Omega*t))
%
% See also STM_PRESET_MODULATED_SLAB, STM_FOURIER_MODULATED_SLAB

if nargin < 3 || isempty(p)
    p = stm_preset_modulated_slab();
end

xGrid = x + zeros(size(t));
tGrid = t + zeros(size(x));
xCell = mod(xGrid, p.Lambda);
isModulated = (xCell >= p.xModStart) & (xCell < p.xModEnd);

epsr = p.eps1*ones(size(xGrid));
modulatedValue = p.epsc*(1 + p.modDepth*sin(p.Omega*tGrid));
epsr(isModulated) = modulatedValue(isModulated);
end
