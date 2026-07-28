function [MM, tau, rho] = temporal_interface_matrix( ...
    epsBefore, muBefore, epsAfter, muAfter)
%TEMPORAL_INTERFACE_MATRIX Morgenthaler matrix at an abrupt time boundary.
%
% The column vector is [E_plus; E_minus]. D and B continuity gives
%   after = MM*before,  MM=[tau rho; rho tau].
%
% This is Eq. (5) of Ramaccia et al., APL 118, 101901 (2021).
%
% D/B continuity is a constitutive-switching model, not a universal rule
% for every experimental time boundary. Use
% TEMPORAL_INTERFACE_MATRIX_JUMP when the switching circuit or microscopic
% model preserves E, injects charge, or otherwise changes the jump laws.

[MM, tau, rho] = temporal_interface_matrix_jump( ...
    epsBefore, muBefore, epsAfter, muAfter, 'DB');
end
