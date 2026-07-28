function amplitudes = temporal_db_to_directional(state, epsr, mur)
%TEMPORAL_DB_TO_DIRECTIONAL Convert [D;B] to [E_plus;E_minus].
%
% This is the inverse of TEMPORAL_DIRECTIONAL_TO_DB.

if size(state,1) ~= 2
    error('state must have size 2-by-N.');
end
if ~isscalar(epsr) || ~isscalar(mur) || epsr == 0 || mur == 0
    error('epsr and mur must be nonzero scalars.');
end

n = sqrt(epsr*mur);
amplitudes = 0.5*[state(1,:)/epsr + state(2,:)/n; ...
    state(1,:)/epsr - state(2,:)/n];
end
