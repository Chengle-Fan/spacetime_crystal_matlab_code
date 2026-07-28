function state = temporal_directional_to_db(amplitudes, epsr, mur)
%TEMPORAL_DIRECTIONAL_TO_DB Convert [E_plus;E_minus] to [D;B].
%
% The convention is exp(i*k*x-i*omega*t), k>0, in normalized units
% eps0=mu0=c0=1. The two temporal-frequency branches obey
%
%   D = epsr*(E_plus+E_minus),
%   B = n*(E_plus-E_minus),  n=sqrt(epsr*mur).

if size(amplitudes,1) ~= 2
    error('amplitudes must have size 2-by-N.');
end
if ~isscalar(epsr) || ~isscalar(mur) || epsr == 0 || mur == 0
    error('epsr and mur must be nonzero scalars.');
end

n = sqrt(epsr*mur);
state = [epsr*(amplitudes(1,:)+amplitudes(2,:)); ...
    n*(amplitudes(1,:)-amplitudes(2,:))];
end
