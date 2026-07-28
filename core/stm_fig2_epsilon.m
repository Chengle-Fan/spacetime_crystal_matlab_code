function epsr = stm_fig2_epsilon(x, t, p)
%STM_FIG2_EPSILON Relative permittivity of the Fig. 2 Floquet crystal.
%
%   epsr = stm_fig2_epsilon(x,t,p)
%
% x and t may be scalars or compatible arrays. The modulated quarter-cell
% obeys eps2(t) = epsc*(1 + modDepth*sin(Omega*t)).

if nargin < 3 || isempty(p)
    p = stm_fig2_parameters();
end

xGrid = x + zeros(size(t));
tGrid = t + zeros(size(x));
xCell = mod(xGrid, p.Lambda);
isModulated = (xCell >= p.xModStart) & (xCell < p.xModEnd);

epsr = p.eps1*ones(size(xGrid));
modulatedValue = p.epsc*(1 + p.modDepth*sin(p.Omega*tGrid));
epsr(isModulated) = modulatedValue(isModulated);
end
