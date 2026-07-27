function field = stpwe_reconstruct_field(mode, x, t, includeGrowth)
%STPWE_RECONSTRUCT_FIELD Reconstruct the real electric field of one mode.
%
% The full field is
% exp(i*k*x-i*omega*t) sum E_mn exp(i*n*g*x-i*m*Omega*t).

if nargin < 4
    includeGrowth = false;
end

[XX, TT] = meshgrid(x, t);
periodicPart = complex(zeros(size(XX)));
for q = 1:numel(mode.Ecoef)
    periodicPart = periodicPart + mode.Ecoef(q).*exp(1i*( ...
        mode.nList(q)*mode.g*XX - mode.mList(q)*mode.Omega*TT));
end

if includeGrowth
    omegaForPlot = mode.omega;
else
    omegaForPlot = real(mode.omega);
end
carrier = exp(1i*(mode.k*XX - omegaForPlot*TT));
field = real(carrier.*periodicPart);

peak = max(abs(field(:)));
if peak > 0
    field = field/peak;
end
end
