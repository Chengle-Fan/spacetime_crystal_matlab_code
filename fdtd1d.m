function out = fdtd1d(cfg)
%FDTD1D  1-D D/B-Yee leapfrog FDTD for scalar, non-dispersive media.
%
% Directly advances the electric displacement D and the magnetic induction B
% (not E and H). This is the correct field pair across an ideal temporal
% interface: D and B are continuous at a permittivity/permeability switch,
% so no separate interface jump matrix is needed. E = D/eps(x,t) and
% H = B/mu(x,t) are recovered by the constitutive relations at each step.
%
% Leapfrog updates (periodic boundary; sponge boundary drops the end cells):
%   Faraday: B^{n+1/2} = B^{n-1/2} - (dt/dx) * (E^n_{i+1} - E^n_i)
%   Ampere : D^{n+1}   = D^n     - (dt/dx) * (H^{n+1/2}_i - H^{n+1/2}_{i-1})
% An explicit temporal interface centered on an integer E-node uses the
% central correction E_interface = 0.5*(D/eps(t-) + D/eps(t+)) when advancing
% B across the switch, avoiding a one-sided (first-order in time) error.
%
% Scope: ordinary positive media only (eps and mu finite positive real).
% Complex eps/mu, negative or zero values are rejected — the stability and
% material contract is defined only for lossless non-dispersive media.
%
% CFL: the stability audit is RIGOROUS. The kernel ALWAYS samples eps on the
% E grid at every integer node (and both sides of each temporal interface)
% and mu on the H grid at every half node, and uses
%   maxWaveSpeed = 1/sqrt(min(eps) * min(mu))
% over those exact states. A caller-supplied cfg.maxWaveSpeed is validated
% as a strict upper bound and cross-checked against this estimate: a value
% below the estimated bound is REJECTED with an error (it would under-state
% the Courant number and hide instability). This covers the actual staggered
% grids and every stepping
% time, so a reported sampledCourant < 1 is a strict guarantee for the
% numerical scheme (a coarse 9-point audit cannot guarantee this and is no
% longer used).
%
% INPUTS  cfg (struct)
%   x             : uniform E-grid, >= 5 points, strictly increasing finite
%                   real vector (row or column).
%   dt            : time step, positive real.
%   nSteps        : nonnegative integer number of leapfrog steps.
%   epsFun(x,t)   : relative permittivity; returns a finite positive real
%                   array the same size as x.
%   muFun(x,t)    : optional, default @(x,t) ones(size(x)) (mu_r = 1).
%   boundary      : 'sponge' (default) or 'periodic'.
%   spongeCells   : absorbing-layer width in cells, nonnegative integer,
%                   default min(80,floor(Nx/8)); the two sponges must not
%                   overlap (2*spongeCells < Nx).
%   spongeStrength: damping strength, finite nonnegative real, default 0.12.
%   E0            : E at t = 0 on the E grid (finite, complex allowed).
%   Hhalf0        : H at t = -dt/2 on the H grid (finite, complex allowed).
%   recordEvery   : record a snapshot every N steps (positive finite integer,
%                   default 1). Snapshots are taken at t = 0 and at every
%                   multiple of recordEvery up to the last one <= nSteps; a
%                   final state not landing on a multiple is omitted.
%   temporalInterfaces : strictly increasing switch times that must align with
%                   integer E-grid nodes n*dt (default: none).
%   temporalInterfaceTolerance : alignment tolerance (default ~dt/100).
%   sourceD(x,t,n): optional soft source: array same size as x, added to D
%                   after the Ampere update at step n (finite, complex
%                   allowed). It is an INCREMENT of D added per step: the
%                   physical current density it represents is sourceD/dt.
%   maxWaveSpeed  : optional validated upper bound on 1/sqrt(eps*mu) over the
%                   whole run (finite, real, strictly positive). The sampled
%                   estimate is always computed and cross-checked against it:
%                   a value below the estimate is rejected with an error.
%
% OUTPUTS (struct)
%   out.x, out.xH  : E grid and H grid (staggered by dx/2).
%   out.t          : nRecords x 1 record times (t = 0 first).
%   out.D, out.E, out.H : nRecords x Nx histories (H interpolated to E grid);
%                   D is ALWAYS stored (no storeD switch).
%   out.dx, out.dt, out.boundary.
%   out.sampledMaxWaveSpeed, out.maxWaveSpeedSource ('user'/'estimated'),
%   out.sampledCourant : CFL audit result.
%   out.temporalInterfaces, out.processedTemporalInterfaceCount.

%==========================================================================
% Grid and step validation (strictly increasing finite real x, dx > 0)
%==========================================================================
if ~isstruct(cfg)
    error('cfg must be a struct.');
end
if ~isfield(cfg, 'x') || isempty(cfg.x)
    error('cfg.x is required.');
end
if ~isvector(cfg.x)
    error('cfg.x must be a vector.');
end
x = cfg.x(:).';
Nx = numel(x);
if Nx < 5
    error('The E grid must contain at least five points.');
end
if ~isreal(x) || ~all(isfinite(x))
    error('cfg.x must contain only finite real values.');
end
if any(diff(x) <= 0)
    error('cfg.x must be strictly increasing (dx > 0); descending grids are not supported.');
end
dxVector = diff(x);
dx = mean(dxVector);
if max(abs(dxVector - dx)) > 1e-10*max(1, abs(dx))
    error('cfg.x must be uniformly spaced.');
end
dt = cfg.dt;
nSteps = cfg.nSteps;
if ~isscalar(dt) || ~isreal(dt) || ~isfinite(dt) || dt <= 0
    error('cfg.dt must be a finite positive real scalar.');
end
if ~isscalar(nSteps) || ~isreal(nSteps) || ~isfinite(nSteps) || ...
        nSteps < 0 || nSteps ~= round(nSteps)
    error('cfg.nSteps must be a nonnegative integer scalar.');
end

%==========================================================================
% Optional-parameter defaults
%==========================================================================
if isfield(cfg, 'muFun') && ~isempty(cfg.muFun)
    muFun = cfg.muFun;
else
    muFun = @(xq,tq) ones(size(xq));
end

if isfield(cfg, 'boundary') && ~isempty(cfg.boundary)
    boundary = lower(cfg.boundary);
else
    boundary = 'sponge';
end
if ~ismember(boundary, {'sponge','periodic'})
    error('cfg.boundary must be ''sponge'' or ''periodic''.');
end

% recordEvery must be a positive finite integer (P1-04): validate before any
% allocation. Snapshots at t=0 and every recordEvery-th step <= nSteps.
if isfield(cfg, 'recordEvery') && ~isempty(cfg.recordEvery)
    recordEvery = cfg.recordEvery;
else
    recordEvery = 1;
end
if ~isscalar(recordEvery) || ~isreal(recordEvery) || ...
        ~isfinite(recordEvery) || recordEvery < 1 || ...
        recordEvery ~= round(recordEvery)
    error('cfg.recordEvery must be a positive finite integer scalar.');
end

%==========================================================================
% Temporal interfaces: align switch times with integer E-grid nodes
%==========================================================================
timeScale = max(1, abs(nSteps*dt));
defaultInterfaceTolerance = min(max(128*eps(timeScale), 1e-12*dt), dt/100);
if isfield(cfg, 'temporalInterfaceTolerance') && ...
        ~isempty(cfg.temporalInterfaceTolerance)
    temporalInterfaceTolerance = cfg.temporalInterfaceTolerance;
    if ~isscalar(temporalInterfaceTolerance) || ...
            ~isreal(temporalInterfaceTolerance) || ...
            ~isfinite(temporalInterfaceTolerance) || ...
            temporalInterfaceTolerance <= 0 || ...
            temporalInterfaceTolerance >= dt/4
        error(['cfg.temporalInterfaceTolerance must be a finite positive ', ...
            'real scalar smaller than cfg.dt/4.']);
    end
else
    temporalInterfaceTolerance = defaultInterfaceTolerance;
end

if isfield(cfg, 'temporalInterfaces') && ~isempty(cfg.temporalInterfaces)
    temporalInterfacesInput = double(cfg.temporalInterfaces(:).');
    if any(~isreal(temporalInterfacesInput)) || any(~isfinite(temporalInterfacesInput)) || ...
            any(diff(temporalInterfacesInput) <= 0)
        error('cfg.temporalInterfaces must be strictly increasing finite reals.');
    end
    temporalInterfaceNodes = round(temporalInterfacesInput/dt);
    temporalInterfaces = temporalInterfaceNodes*dt;
    if any(abs(temporalInterfacesInput - temporalInterfaces) > ...
            temporalInterfaceTolerance)
        error(['Each cfg.temporalInterfaces value must align with an E-grid ', ...
            'time n*cfg.dt within cfg.temporalInterfaceTolerance.']);
    end
    if any(temporalInterfaceNodes < 0) || ...
            any(temporalInterfaceNodes >= nSteps)
        error(['cfg.temporalInterfaces must lie in [0, cfg.nSteps*cfg.dt). ', ...
            'An interface at the final time has no following B half-step.']);
    end
    if any(diff(temporalInterfaceNodes) <= 0)
        error('Distinct cfg.temporalInterfaces must map to distinct E-grid nodes.');
    end
else
    temporalInterfaces = zeros(1, 0);
    temporalInterfaceNodes = zeros(1, 0);
end

temporalInterfaceIdByStep = zeros(1, nSteps, 'uint32');
if ~isempty(temporalInterfaceNodes)
    temporalInterfaceIdByStep(temporalInterfaceNodes + 1) = ...
        uint32(1:numel(temporalInterfaceNodes));
end
temporalInterfaceSideOffset = dt/4;

%==========================================================================
% Yee grids and initial conditions (size + finiteness validated)
%==========================================================================
switch boundary
    case 'periodic'
        xH = x + dx/2;                    % Nx points (wrapped curl)
    case 'sponge'
        xH = x(1:end-1) + dx/2;           % Nx-1 points (open curls)
end

E = eval_initial_field(cfg, 'E0', x, 0);              % complex allowed
H = eval_initial_field(cfg, 'Hhalf0', xH, -dt/2);

% Constitutive init: D = eps*E at t=0, B = mu*H at t=-dt/2.
epsNow = eval_material(cfg.epsFun, x, 0, 'epsFun');
muHalf = eval_material(muFun, xH, -dt/2, 'muFun');
D = epsNow.*E;
B = muHalf.*H;

%==========================================================================
% Rigorous CFL audit (P1-02/P1-03): covers the actual E/H grids and every
% stepping time. Caller-provided maxWaveSpeed is validated as a strict
% positive finite real upper bound; otherwise it is estimated from the min
% eps on the E grid (all integer nodes + both sides of each interface) and
% the min mu on the H grid (all half nodes): 1/sqrt(min_eps * min_mu).
%==========================================================================
% The estimated bound is ALWAYS computed from the actual E/H grids so a
% caller-supplied bound can be cross-checked against it (a user value below
% the true 1/sqrt(min_eps*min_mu) would under-state the Courant number and
% hide instability). The estimate covers every stepping time on the E grid
% (all integer nodes + both sides of each temporal interface) and every B
% half-step on the H grid (all half nodes).
epsMin = inf;
auditTimesE = (0:nSteps)*dt;
if ~isempty(temporalInterfaceNodes)
    auditTimesE = [auditTimesE, temporalInterfaces - dt/4, ...
        temporalInterfaces + dt/4];
end
for tn = auditTimesE
    v = eval_material(cfg.epsFun, x, tn, 'epsFun');
    epsMin = min(epsMin, min(v));
end
muMin = inf;
auditTimesH = ((1:nSteps) - 0.5)*dt;
for th = auditTimesH
    v = eval_material(muFun, xH, th, 'muFun');
    muMin = min(muMin, min(v));
end
maxWaveSpeedEst = 1/sqrt(epsMin*muMin);

if isfield(cfg, 'maxWaveSpeed') && ~isempty(cfg.maxWaveSpeed)
    maxWaveSpeed = cfg.maxWaveSpeed;
    if ~isscalar(maxWaveSpeed) || ~isreal(maxWaveSpeed) || ...
            ~isfinite(maxWaveSpeed) || maxWaveSpeed <= 0
        error(['cfg.maxWaveSpeed (caller-provided upper bound on ', ...
            '1/sqrt(eps*mu)) must be a finite, real, strictly positive ', ...
            'scalar.']);
    end
    if maxWaveSpeed < maxWaveSpeedEst*(1 - 1e-9)
        error(['cfg.maxWaveSpeed = %.6g is BELOW the estimated grid-bound ', ...
            '1/sqrt(min_eps*min_mu) = %.6g; the effective CFL number would ', ...
            'be under-stated and the leapfrog scheme may be unstable. A ', ...
            'caller-supplied cfg.maxWaveSpeed must be a valid UPPER bound ', ...
            'on the sampled estimate.'], maxWaveSpeed, maxWaveSpeedEst);
    end
    maxWaveSpeedSource = 'user';
else
    maxWaveSpeed = maxWaveSpeedEst;
    maxWaveSpeedSource = 'estimated';
end
sampledCourant = maxWaveSpeed*dt/dx;
if ~isfinite(sampledCourant)
    error(['CFL audit produced a non-finite Courant number (%.6g). ', ...
        'Check the material and grid.'], sampledCourant);
end
if sampledCourant >= 1
    error(['One-dimensional CFL number is %.6g >= 1 using the %s upper ', ...
        'bound maxWaveSpeed = %.6g. Reduce cfg.dt or provide a validated ', ...
        'tighter cfg.maxWaveSpeed.'], sampledCourant, maxWaveSpeedSource, ...
        maxWaveSpeed);
end

%==========================================================================
% Sponge absorbing layer (simplified absorption, NOT a PML)
%==========================================================================
dampE = ones(size(x));
dampH = ones(size(xH));
if strcmp(boundary, 'sponge')
    if isfield(cfg, 'spongeCells') && ~isempty(cfg.spongeCells)
        spongeCells = cfg.spongeCells;
    else
        spongeCells = min(80, floor(Nx/8));
    end
    if ~isscalar(spongeCells) || ~isreal(spongeCells) || ...
            ~isfinite(spongeCells) || spongeCells < 0 || ...
            spongeCells ~= round(spongeCells)
        error('cfg.spongeCells must be a nonnegative finite integer scalar.');
    end
    if 2*spongeCells >= Nx
        error(['cfg.spongeCells = %d makes the two absorbing layers ', ...
            'overlap on a %d-point grid (need 2*spongeCells < Nx).'], ...
            spongeCells, Nx);
    end
    if isfield(cfg, 'spongeStrength') && ~isempty(cfg.spongeStrength)
        spongeStrength = cfg.spongeStrength;
    else
        spongeStrength = 0.12;
    end
    if ~isscalar(spongeStrength) || ~isreal(spongeStrength) || ...
            ~isfinite(spongeStrength) || spongeStrength < 0
        error('cfg.spongeStrength must be a finite nonnegative real scalar.');
    end
    edgeDistanceE = min((x - x(1))/dx, (x(end) - x)/dx);
    edgeDistanceH = min((xH - x(1))/dx, (x(end) - xH)/dx);
    maskE = edgeDistanceE < spongeCells;
    maskH = edgeDistanceH < spongeCells;
    sE = (spongeCells - edgeDistanceE(maskE))/spongeCells;
    sH = (spongeCells - edgeDistanceH(maskH))/spongeCells;
    dampE(maskE) = exp(-spongeStrength*sE.^3);
    dampH(maskH) = exp(-spongeStrength*sH.^3);
end

%==========================================================================
% Record arrays and initial snapshot (exact record index bookkeeping)
%==========================================================================
recordSteps = recordEvery:recordEvery:nSteps;
nRecords = 1 + numel(recordSteps);
DHistory = complex(zeros(nRecords, Nx));
EHistory = complex(zeros(nRecords, Nx));
HHistory = complex(zeros(nRecords, Nx));
tHistory = zeros(nRecords, 1);

recordId = 1;
HOnE = yee_h_to_e_grid(H, boundary, Nx);
DHistory(recordId, :) = D;
EHistory(recordId, :) = E;
HHistory(recordId, :) = HOnE;
tHistory(recordId) = 0;
recordStepIndex = 1;                 % next target step = recordSteps(idx)

%==========================================================================
% Leapfrog time-stepping
%==========================================================================
processedTemporalInterfaceCount = 0;
for step = 1:nSteps
    % ---- Faraday: B^{n+1/2} from E^n (with central interface correction) ----
    temporalInterfaceId = temporalInterfaceIdByStep(step);
    if temporalInterfaceId ~= 0
        processedTemporalInterfaceCount = ...
            processedTemporalInterfaceCount + 1;
        tInterface = temporalInterfaces(double(temporalInterfaceId));
        epsBefore = eval_material(cfg.epsFun, x, ...
            tInterface - temporalInterfaceSideOffset, 'epsFun');
        epsAfter  = eval_material(cfg.epsFun, x, ...
            tInterface + temporalInterfaceSideOffset, 'epsFun');
        EForBUpdate = 0.5*(D./epsBefore + D./epsAfter);
    else
        EForBUpdate = E;
    end

    if strcmp(boundary, 'periodic')
        curlE = circshift(EForBUpdate, -1) - EForBUpdate;
    else
        curlE = diff(EForBUpdate);
    end
    B = B - (dt/dx)*curlE;
    B = B.*dampH;
    tHalf = (step - 0.5)*dt;
    H = B./eval_material(muFun, xH, tHalf, 'muFun');

    % ---- Ampere: D^{n+1} from H^{n+1/2} ----
    if strcmp(boundary, 'periodic')
        curlH = H - circshift(H, 1);
        D = D - (dt/dx)*curlH;
    else
        D(2:end-1) = D(2:end-1) - (dt/dx)*(H(2:end) - H(1:end-1));
    end

    tNow = step*dt;

    if isfield(cfg, 'sourceD') && ~isempty(cfg.sourceD)
        D = D + eval_source(cfg.sourceD, x, tNow, step, 'sourceD');
    end
    D = D.*dampE;
    E = D./eval_material(cfg.epsFun, x, tNow, 'epsFun');

    % ---- record ----
    if recordStepIndex <= numel(recordSteps) && step == recordSteps(recordStepIndex)
        recordId = recordId + 1;
        HOnE = yee_h_to_e_grid(H, boundary, Nx);
        DHistory(recordId, :) = D;
        EHistory(recordId, :) = E;
        HHistory(recordId, :) = HOnE;
        tHistory(recordId) = tNow;
        recordStepIndex = recordStepIndex + 1;
    end
end

if processedTemporalInterfaceCount ~= numel(temporalInterfaces)
    error(['Processed %d temporal interfaces, but %d aligned interfaces ', ...
        'were configured.'], processedTemporalInterfaceCount, ...
        numel(temporalInterfaces));
end

%==========================================================================
% Output
%==========================================================================
out.x = x;
out.xH = xH;
out.t = tHistory;
out.D = DHistory;
out.E = EHistory;
out.H = HHistory;
out.dx = dx;
out.dt = dt;
out.boundary = boundary;
out.sampledMaxWaveSpeed = maxWaveSpeed;
out.maxWaveSpeedSource = maxWaveSpeedSource;
out.sampledCourant = sampledCourant;
out.temporalInterfaces = temporalInterfaces;
out.processedTemporalInterfaceCount = processedTemporalInterfaceCount;
end

%==========================================================================
function v = eval_material(fun, grid, tval, label)
%EVAL_MATERIAL Evaluate eps/mu on a grid and enforce the positive-real,
%grid-size and finiteness contract (P1-05).
if ~isa(fun, 'function_handle')
    error('%s must be a function handle.', label);
end
v = fun(grid, tval);
if ~isequal(size(v), size(grid))
    error(['%s must return an array the same size as its grid argument ', ...
        '(got %s on a %s grid).'], label, mat2str(size(v)), mat2str(size(grid)));
end
v = v(:).';
if ~all(isfinite(v))
    error('%s returned non-finite values at t = %.16g.', label, tval);
end
if any(~isreal(v)) || any(real(v) <= 0)
    error(['%s must return positive real values (ordinary non-dispersive ', ...
        'lossless media) at t = %.16g.'], label, tval);
end
end

%==========================================================================
function v = eval_source(fun, grid, tval, step, label)
%EVAL_SOURCE Evaluate a soft source callback; finite, size-matched, complex
%allowed (P1-05).
if ~isa(fun, 'function_handle')
    error('%s must be a function handle.', label);
end
v = fun(grid, tval, step);
if ~isequal(size(v), size(grid))
    error(['%s must return an array the same size as its grid argument ', ...
        '(got %s on a %s grid).'], label, mat2str(size(v)), mat2str(size(grid)));
end
v = v(:).';
if ~all(isfinite(v))
    error('%s returned non-finite values at step %d (t = %.16g).', ...
        label, step, tval);
end
end

%==========================================================================
function f = eval_initial_field(cfg, fieldName, grid, tval)
%EVAL_INITIAL_FIELD Read a finite complex-allowed initial field of exact size.
if isfield(cfg, fieldName) && ~isempty(cfg.(fieldName))
    v = cfg.(fieldName);
    if ~isvector(v)
        error('cfg.%s must be a vector (no implicit matrix flattening).', fieldName);
    end
    if numel(v) ~= numel(grid)
        error('cfg.%s has %d elements; the selected boundary needs %d.', ...
            fieldName, numel(v), numel(grid));
    end
    f = v(:).';
    if ~all(isfinite(f))
        error('cfg.%s must be finite at t = %.16g.', fieldName, tval);
    end
else
    f = zeros(size(grid));
end
end

%==========================================================================
function HOnE = yee_h_to_e_grid(H, boundary, Nx)
%YEE_H_TO_E_GRID Interpolate H from the staggered H grid onto the E grid.
if strcmp(boundary, 'periodic')
    HOnE = 0.5*(H + circshift(H, 1));
else
    HOnE = zeros(1, Nx, 'like', H);
    HOnE(2:end-1) = 0.5*(H(1:end-1) + H(2:end));
    HOnE(1) = H(1);
    HOnE(end) = H(end);
end
end
