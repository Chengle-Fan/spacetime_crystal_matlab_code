function out = fdtd1d_db(cfg)
%FDTD1D_DB One-dimensional Yee FDTD for space-time-varying media.
%
% Required cfg fields:
%   x             E/D grid (uniform row or column vector)
%   dt            time step
%   nSteps        number of updates
%   epsFun(x,t)   relative permittivity
%
% Optional cfg fields:
%   muFun(x,t)       relative permeability; default 1
%   E0               initial E on the E grid; default 0
%   Hhalf0           initial H at t=-dt/2 on the H grid; default 0
%   boundary         'sponge' (default) or 'periodic'
%   spongeCells      default min(80,floor(Nx/8))
%   spongeStrength   per-step damping strength, default 0.12
%   recordEvery      default 1
%   storeFields      store full E/H histories; default true
%   probeIndices     E-grid indices recorded even if storeFields=false
%   stabilityTimes   times sampled for a positive-medium CFL audit
%   maxWaveSpeed     optional certified maximum wave speed
%   sourceD(x,t,n)   optional additive D increment after each update
%
% D and B, rather than E and H, are advanced. Consequently an abrupt
% temporal jump automatically preserves the correct temporal boundary
% conditions D(t+)=D(t-) and B(t+)=B(t-).

x = cfg.x(:).';
Nx = numel(x);
if Nx < 5
    error('The E grid must contain at least five points.');
end
dxVector = diff(x);
dx = mean(dxVector);
if max(abs(dxVector-dx)) > 1e-10*max(1,abs(dx))
    error('cfg.x must be uniformly spaced.');
end
dt = cfg.dt;
nSteps = cfg.nSteps;

if ~isfield(cfg,'muFun') || isempty(cfg.muFun)
    muFun = @(xq,tq) ones(size(xq));
else
    muFun = cfg.muFun;
end
if ~isfield(cfg,'boundary') || isempty(cfg.boundary)
    boundary = 'sponge';
else
    boundary = lower(cfg.boundary);
end
if ~isfield(cfg,'recordEvery') || isempty(cfg.recordEvery)
    recordEvery = 1;
else
    recordEvery = cfg.recordEvery;
end
if ~isfield(cfg,'storeFields') || isempty(cfg.storeFields)
    storeFields = true;
else
    storeFields = logical(cfg.storeFields);
end
if ~isfield(cfg,'probeIndices') || isempty(cfg.probeIndices)
    probeIndices = zeros(1,0);
else
    probeIndices = unique(cfg.probeIndices(:).');
    if any(probeIndices < 1) || any(probeIndices > Nx) || ...
            any(probeIndices ~= round(probeIndices))
        error('cfg.probeIndices must contain valid integer E-grid indices.');
    end
end

switch boundary
    case 'periodic'
        xH = x + dx/2;
    case 'sponge'
        xH = x(1:end-1) + dx/2;
    otherwise
        error('cfg.boundary must be ''periodic'' or ''sponge''.');
end

if isfield(cfg,'E0') && ~isempty(cfg.E0)
    E = cfg.E0(:).';
else
    E = zeros(1,Nx);
end
if numel(E) ~= Nx
    error('cfg.E0 must have the same length as cfg.x.');
end

if isfield(cfg,'Hhalf0') && ~isempty(cfg.Hhalf0)
    H = cfg.Hhalf0(:).';
else
    H = zeros(size(xH));
end
if numel(H) ~= numel(xH)
    error('cfg.Hhalf0 has the wrong length for the selected boundary.');
end

epsNow = cfg.epsFun(x,0);
muHalf = muFun(xH,-dt/2);
if numel(epsNow) ~= Nx || numel(muHalf) ~= numel(xH)
    error('epsFun or muFun returned an array with the wrong grid size.');
end
D = epsNow.*E;
B = muHalf.*H;

% A sampled CFL audit is useful for ordinary positive nondispersive media.
% For dispersive, negative-index, or active constitutive models the scalar
% phase speed is not a sufficient stability criterion.
if isfield(cfg,'maxWaveSpeed') && ~isempty(cfg.maxWaveSpeed)
    sampledMaxWaveSpeed = cfg.maxWaveSpeed;
else
    if isfield(cfg,'stabilityTimes') && ~isempty(cfg.stabilityTimes)
        stabilityTimes = cfg.stabilityTimes(:).';
    else
        stabilityTimes = linspace(0,nSteps*dt,9);
    end
    sampledMaxWaveSpeed = 0;
    for tq = stabilityTimes
        epsSample = cfg.epsFun(x,tq);
        muSample = muFun(x,tq);
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

dampE = ones(size(x));
dampH = ones(size(xH));
if strcmp(boundary,'sponge')
    if ~isfield(cfg,'spongeCells') || isempty(cfg.spongeCells)
        spongeCells = min(80,floor(Nx/8));
    else
        spongeCells = cfg.spongeCells;
    end
    if ~isfield(cfg,'spongeStrength') || isempty(cfg.spongeStrength)
        spongeStrength = 0.12;
    else
        spongeStrength = cfg.spongeStrength;
    end
    edgeDistanceE = min((x-x(1))/dx, (x(end)-x)/dx);
    edgeDistanceH = min((xH-x(1))/dx, (x(end)-xH)/dx);
    maskE = edgeDistanceE < spongeCells;
    maskH = edgeDistanceH < spongeCells;
    sE = (spongeCells-edgeDistanceE(maskE))/spongeCells;
    sH = (spongeCells-edgeDistanceH(maskH))/spongeCells;
    dampE(maskE) = exp(-spongeStrength*sE.^3);
    dampH(maskH) = exp(-spongeStrength*sH.^3);
end

nRecords = floor(nSteps/recordEvery) + 1;
if storeFields
    EHistory = complex(zeros(nRecords,Nx));
    HHistory = complex(zeros(nRecords,Nx));
else
    EHistory = complex(zeros(0,Nx));
    HHistory = complex(zeros(0,Nx));
end
probeE = complex(zeros(nRecords,numel(probeIndices)));
probeH = complex(zeros(nRecords,numel(probeIndices)));
tHistory = zeros(nRecords,1);
energy = zeros(nRecords,1);

recordId = 1;
[HOnE, BOnE] = yee_h_to_e_grid(H,B,boundary,Nx);
if storeFields
    EHistory(recordId,:) = E;
    HHistory(recordId,:) = HOnE;
end
probeE(recordId,:) = E(probeIndices);
probeH(recordId,:) = HOnE(probeIndices);
energy(recordId) = 0.5*dx*sum(real(conj(E).*D + conj(HOnE).*BOnE));

for step = 1:nSteps
    % B at (n+1/2)dt from E at n*dt.
    if strcmp(boundary,'periodic')
        curlE = circshift(E,-1) - E;
    else
        curlE = diff(E);
    end
    B = B - (dt/dx)*curlE;
    B = B.*dampH;

    tHalf = (step-0.5)*dt;
    muHalf = muFun(xH,tHalf);
    H = B./muHalf;

    % D at (n+1)dt from H at (n+1/2)dt.
    if strcmp(boundary,'periodic')
        curlH = H - circshift(H,1);
        D = D - (dt/dx)*curlH;
    else
        D(2:end-1) = D(2:end-1) ...
            - (dt/dx)*(H(2:end)-H(1:end-1));
    end

    tNow = step*dt;
    if isfield(cfg,'sourceD') && ~isempty(cfg.sourceD)
        D = D + cfg.sourceD(x,tNow,step);
    end
    D = D.*dampE;
    epsNow = cfg.epsFun(x,tNow);
    E = D./epsNow;

    if mod(step,recordEvery) == 0
        recordId = recordId + 1;
        [HOnE, BOnE] = yee_h_to_e_grid(H,B,boundary,Nx);
        if storeFields
            EHistory(recordId,:) = E;
            HHistory(recordId,:) = HOnE;
        end
        probeE(recordId,:) = E(probeIndices);
        probeH(recordId,:) = HOnE(probeIndices);
        tHistory(recordId) = tNow;
        energy(recordId) = 0.5*dx*sum(real( ...
            conj(E).*D + conj(HOnE).*BOnE));
    end
end

out.x = x;
out.xH = xH;
out.t = tHistory;
out.E = EHistory;
out.H = HHistory;
out.probeIndices = probeIndices;
out.probeX = x(probeIndices);
out.probeE = probeE;
out.probeH = probeH;
out.energy = energy;
out.finalD = D;
out.finalB = B;
out.finalE = E;
out.finalH = H;
out.dx = dx;
out.dt = dt;
out.boundary = boundary;
out.storeFields = storeFields;
out.sampledMaxWaveSpeed = sampledMaxWaveSpeed;
out.sampledCourant = sampledCourant;
end

function [HOnE, BOnE] = yee_h_to_e_grid(H,B,boundary,Nx)
if strcmp(boundary,'periodic')
    HOnE = 0.5*(H + circshift(H,1));
    BOnE = 0.5*(B + circshift(B,1));
else
    HOnE = zeros(1,Nx,'like',H);
    BOnE = zeros(1,Nx,'like',B);
    HOnE(2:end-1) = 0.5*(H(1:end-1)+H(2:end));
    BOnE(2:end-1) = 0.5*(B(1:end-1)+B(2:end));
    HOnE(1) = H(1);
    HOnE(end) = H(end);
    BOnE(1) = B(1);
    BOnE(end) = B(end);
end
end
