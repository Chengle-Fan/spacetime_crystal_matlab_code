function mode = stpwe_select_mode(sys, kTarget, omegaTarget, omegaWindow)
%STPWE_SELECT_MODE Select the visible Floquet mode nearest a target point.

if nargin < 4 || isempty(omegaWindow)
    omegaWindow = [-inf, inf];
end

sol = stpwe_solve_omega(sys, kTarget);
valid = isfinite(sol.omega) ...
    & real(sol.omega) >= omegaWindow(1) ...
    & real(sol.omega) <= omegaWindow(2);
ids = find(valid);
if isempty(ids)
    error('No finite ST-PWE mode exists in the requested window.');
end

scale = max(abs(omegaTarget), sys.Omega);
score = abs(real(sol.omega(ids)) - omegaTarget)/scale ...
    + 0.15*abs(imag(sol.omega(ids)))/scale ...
    + 0.02*(1 - sol.m0Weight(ids));
[~, localId] = min(score);
id = ids(localId);

vec = sol.R(:,id);
[~, phaseId] = max(abs(vec(1:sys.S)));
vec = vec*exp(-1i*angle(vec(phaseId)));

mode.k = kTarget;
mode.omega = sol.omega(id);
mode.Ecoef = vec(1:sys.S);
mode.Hcoef = vec(sys.S+1:end);
mode.nList = sys.nList;
mode.mList = sys.mList;
mode.g = sys.g;
mode.Omega = sys.Omega;
mode.m0Weight = sol.m0Weight(id);
end
