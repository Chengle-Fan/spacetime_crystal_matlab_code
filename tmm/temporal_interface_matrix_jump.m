function [MM, tau, rho, details] = temporal_interface_matrix_jump( ...
    epsBefore, muBefore, epsAfter, muAfter, model)
%TEMPORAL_INTERFACE_MATRIX_JUMP Time-interface matrix with explicit jumps.
%
%   after = MM*before
%
% acts on [E_plus; E_minus]. The jump factors are defined by
%
%   D_after = jumpD*D_before,
%   B_after = jumpB*B_before.
%
% MODEL may be:
%   'DB'  D and B continuous (standard bulk time-boundary model)
%   'EB'  E and B continuous
%   'DH'  D and H continuous
%   'EH'  E and H continuous
%   struct('jumpD',qD,'jumpB',qB) for a user-supplied microscopic model
%
% The named E/B/D/H choices are idealized jump laws. A switched circuit
% containing sources, disconnected charge, spatially localized switching,
% or temporal dispersion may require a larger state and cannot in general
% be represented by two scalar jump factors.

if nargin < 5 || isempty(model)
    model = 'DB';
end
if epsBefore == 0 || muBefore == 0 || epsAfter == 0 || muAfter == 0
    error('Constitutive parameters must be nonzero at the interface.');
end

if isstruct(model)
    if ~isfield(model,'jumpD') || ~isfield(model,'jumpB')
        error('A custom model must contain jumpD and jumpB.');
    end
    jumpD = model.jumpD;
    jumpB = model.jumpB;
    label = 'custom';
elseif ischar(model) || (isstring(model) && isscalar(model))
    label = upper(regexprep(char(model),'[^A-Za-z]',''));
    switch label
        case 'DB'
            jumpD = 1;
            jumpB = 1;
        case 'EB'
            jumpD = epsAfter/epsBefore;
            jumpB = 1;
        case 'DH'
            jumpD = 1;
            jumpB = muAfter/muBefore;
        case 'EH'
            jumpD = epsAfter/epsBefore;
            jumpB = muAfter/muBefore;
        otherwise
            error('Unknown temporal-interface jump model ''%s''.',label);
    end
else
    error('model must be a named jump law or a jump-factor struct.');
end

ratioElectric = jumpD*epsBefore/epsAfter;
ratioMagnetic = jumpB*sqrt( ...
    (epsBefore*muBefore)/(epsAfter*muAfter));
tau = 0.5*(ratioElectric + ratioMagnetic);
rho = 0.5*(ratioElectric - ratioMagnetic);
MM = [tau, rho; rho, tau];

details.model = label;
details.jumpD = jumpD;
details.jumpB = jumpB;
details.electricRatio = ratioElectric;
details.magneticRatio = ratioMagnetic;
end
