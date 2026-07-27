function [MM, tau, rho] = temporal_interface_matrix( ...
    epsBefore, muBefore, epsAfter, muAfter)
%TEMPORAL_INTERFACE_MATRIX Morgenthaler matrix at an abrupt time boundary.
%
% The column vector is [E_plus; E_minus]. D and B continuity gives
%   after = MM*before,  MM=[tau rho; rho tau].
%
% This is Eq. (5) of Ramaccia et al., APL 118, 101901 (2021).

ratioD = epsBefore/epsAfter;
ratioDB = sqrt((epsBefore*muBefore)/(epsAfter*muAfter));
tau = 0.5*(ratioD + ratioDB);
rho = 0.5*(ratioD - ratioDB);
MM = [tau, rho; rho, tau];
end
