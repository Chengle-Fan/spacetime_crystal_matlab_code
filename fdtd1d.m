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
% INPUTS  cfg (struct)
%   x             : uniform E-grid, >= 5 points (row or column).
%   dt            : time step, positive real.
%   nSteps        : nonnegative integer number of leapfrog steps.
%   epsFun(x,t)   : relative permittivity; returns an array the same size as x.
%   muFun(x,t)    : optional, default @(x,t) ones(size(x)) (mu_r = 1).
%   boundary      : 'sponge' (default) or 'periodic'.
%   spongeCells   : absorbing-layer width in cells, default min(80,floor(Nx/8)).
%   spongeStrength: damping strength, default 0.12.
%   E0            : E at t = 0 on the E grid (default 0).
%   Hhalf0        : H at t = -dt/2 on the H grid (default 0).
%   recordEvery   : record a snapshot every N steps (default 1).
%   temporalInterfaces : strictly increasing switch times that must align with
%                   integer E-grid nodes n*dt (default: none).
%   temporalInterfaceTolerance : alignment tolerance (default ~dt/100).
%   sourceD(x,t,n): optional soft source: array same size as x, added to D
%                   after the Ampere update at step n.
%   maxWaveSpeed  : optional manual CFL bound (skips the sampled audit).
%   storeD        : accepted for compatibility; D history is always stored.
%
% OUTPUTS (struct)
%   out.x, out.xH  : E grid and H grid (staggered by dx/2).
%   out.t          : nRecords x 1 record times.
%   out.D, out.E, out.H : nRecords x Nx histories (H interpolated to E grid).
%   out.dx, out.dt, out.boundary.
%   out.sampledMaxWaveSpeed, out.sampledCourant : CFL audit result.
%   out.temporalInterfaces, out.processedTemporalInterfaceCount.

%==========================================================================
% Grid and step validation
%==========================================================================
x = cfg.x(:).';
Nx = numel(x);
if Nx < 5
    error('The E grid must contain at least five points.');
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

if isfield(cfg, 'recordEvery') && ~isempty(cfg.recordEvery)
    recordEvery = cfg.recordEvery;
else
    recordEvery = 1;
end

if isfield(cfg, 'E0') && ~isempty(cfg.E0)
    E = cfg.E0(:).';
else
    E = zeros(1, Nx);
end
if numel(E) ~= Nx
    error('cfg.E0 must have the same length as cfg.x.');
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
    if any(~isfinite(temporalInterfacesInput)) || ...
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
% Yee grids and initial conditions
%==========================================================================
switch boundary
    case 'periodic'
        xH = x + dx/2;                    % Nx points (wrapped curl)
    case 'sponge'
        xH = x(1:end-1) + dx/2;           % Nx-1 points (open curls)
end

if isfield(cfg, 'Hhalf0') && ~isempty(cfg.Hhalf0)
    H = cfg.Hhalf0(:).';
else
    H = zeros(size(xH));
end
if numel(H) ~= numel(xH)
    error('cfg.Hhalf0 has the wrong length for the selected boundary.');
end

% Constitutive init: D = eps*E at t=0, B = mu*H at t=-dt/2.
epsNow = cfg.epsFun(x, 0);
muHalf = muFun(xH, -dt/2);
if numel(epsNow) ~= Nx || numel(muHalf) ~= numel(xH)
    error('epsFun or muFun returned an array with the wrong grid size.');
end
D = epsNow.*E;
B = muHalf.*H;

%==========================================================================
% Sampled CFL audit (scalar non-dispersive media)
%==========================================================================
if isfield(cfg, 'maxWaveSpeed') && ~isempty(cfg.maxWaveSpeed)
    sampledMaxWaveSpeed = cfg.maxWaveSpeed;
else
    stabilityTimes = linspace(0, nSteps*dt, 9);
    sampledMaxWaveSpeed = 0;
    for tq = stabilityTimes
        epsSample = cfg.epsFun(x, tq);
        muSample = muFun(x, tq);
        if numel(epsSample) ~= Nx || numel(muSample) ~= Nx
            error(['epsFun and muFun must return the same size as their ', ...
                'input grid during the CFL audit.']);
        end
        ordinaryPositive = all(abs(imag(epsSample)) < 1e-12) ...
            && all(abs(imag(muSample)) < 1e-12) ...
            && all(real(epsSample) > 0) && all(real(muSample) > 0);
        if ~ordinaryPositive
            sampledMaxWaveSpeed = nan;
            break;
        end
        sampledMaxWaveSpeed = max(sampledMaxWaveSpeed, ...
            max(1./sqrt(real(epsSample).*real(muSample))));
    end
end
sampledCourant = sampledMaxWaveSpeed*dt/dx;
if isfinite(sampledCourant) && sampledCourant >= 1
    error(['Sampled one-dimensional CFL number is %.6g >= 1. ', ...
        'Reduce cfg.dt or provide a validated constitutive update.'], ...
        sampledCourant);
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
    if isfield(cfg, 'spongeStrength') && ~isempty(cfg.spongeStrength)
        spongeStrength = cfg.spongeStrength;
    else
        spongeStrength = 0.12;
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
% Record arrays and initial snapshot
%==========================================================================
nRecords = floor(nSteps/recordEvery) + 1;
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
        epsBefore = cfg.epsFun(x, tInterface - temporalInterfaceSideOffset);
        epsAfter  = cfg.epsFun(x, tInterface + temporalInterfaceSideOffset);
        if numel(epsBefore) ~= Nx || numel(epsAfter) ~= Nx
            error(['epsFun returned an array with the wrong grid size ', ...
                'around temporal interface %.16g.'], tInterface);
        end
        epsBefore = epsBefore(:).';
        epsAfter = epsAfter(:).';
        if any(abs(epsBefore) == 0) || any(abs(epsAfter) == 0)
            error(['epsFun must return finite nonzero permittivity on both ', ...
                'sides of temporal interface %.16g.'], tInterface);
        end
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
    H = B./muFun(xH, tHalf);

    % ---- Ampere: D^{n+1} from H^{n+1/2} ----
    if strcmp(boundary, 'periodic')
        curlH = H - circshift(H, 1);
        D = D - (dt/dx)*curlH;
    else
        D(2:end-1) = D(2:end-1) - (dt/dx)*(H(2:end) - H(1:end-1));
    end

    tNow = step*dt;

    if isfield(cfg, 'sourceD') && ~isempty(cfg.sourceD)
        D = D + cfg.sourceD(x, tNow, step);
    end
    D = D.*dampE;
    E = D./cfg.epsFun(x, tNow);

    % ---- record ----
    if mod(step, recordEvery) == 0
        recordId = recordId + 1;
        HOnE = yee_h_to_e_grid(H, boundary, Nx);
        DHistory(recordId, :) = D;
        EHistory(recordId, :) = E;
        HHistory(recordId, :) = HOnE;
        tHistory(recordId) = tNow;
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
out.sampledMaxWaveSpeed = sampledMaxWaveSpeed;
out.sampledCourant = sampledCourant;
out.temporalInterfaces = temporalInterfaces;
out.processedTemporalInterfaceCount = processedTemporalInterfaceCount;
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
