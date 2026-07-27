function band = stpwe_track_band(sys, kGrid, omegaSeed, options)
%STPWE_TRACK_BAND Track one Floquet band by eigenvector continuity.
%
% options.omegaWindow = [minimum maximum]
% options.maxImag      = maximum |Im(omega)|
% options.minM0Weight  = minimum m=0 participation
%
% A low overlap or a closing neighbor gap signals that a single-band Zak
% phase is not reliable; use a composite-band Wilson loop in that case.

if nargin < 4
    options = struct();
end
if ~isfield(options,'omegaWindow')
    options.omegaWindow = [-inf inf];
end
if ~isfield(options,'maxImag')
    options.maxImag = inf;
end
if ~isfield(options,'minM0Weight')
    options.minM0Weight = 0;
end

kGrid = kGrid(:).';
Nk = numel(kGrid);
nState = 2*sys.S;
Rtrack = complex(zeros(nState,Nk));
Ltrack = complex(zeros(nState,Nk));
omegaTrack = complex(zeros(1,Nk));
weightTrack = zeros(1,Nk);
selectedId = zeros(1,Nk);
continuity = ones(1,Nk);

for ik = 1:Nk
    sol = stpwe_solve_omega(sys,kGrid(ik));
    valid = isfinite(sol.omega) ...
        & real(sol.omega) >= options.omegaWindow(1) ...
        & real(sol.omega) <= options.omegaWindow(2) ...
        & abs(imag(sol.omega)) <= options.maxImag ...
        & sol.m0Weight >= options.minM0Weight;
    ids = find(valid);
    if isempty(ids)
        error('No candidate mode remains at k index %d.',ik);
    end

    if ik == 1
        scale = max(abs(omegaSeed),sys.Omega);
        cost = abs(sol.omega(ids)-omegaSeed)/scale ...
            + 0.01*(1-sol.m0Weight(ids));
    else
        prevR = Rtrack(:,ik-1);
        overlap = complex(zeros(size(ids)));
        for q = 1:numel(ids)
            candidate = sol.R(:,ids(q));
            overlap(q) = (prevR'*candidate) ...
                /(norm(prevR)*norm(candidate));
        end
        frequencyStep = abs(sol.omega(ids)-omegaTrack(ik-1)) ...
            /max(sys.Omega,abs(omegaTrack(ik-1)));
        cost = 1-abs(overlap) + 0.15*frequencyStep ...
            + 0.01*(1-sol.m0Weight(ids));
    end

    [~,localId] = min(cost);
    id = ids(localId);
    R = sol.R(:,id);
    L = sol.L(:,id);

    if ik > 1
        link = Ltrack(:,ik-1)'*sys.Bomega*R;
        continuity(ik) = abs(link);
        if abs(link) > 1e-14
            phase = exp(-1i*angle(link));
            R = R*phase;
            L = L*phase;
        end
    end

    Rtrack(:,ik) = R;
    Ltrack(:,ik) = L;
    omegaTrack(ik) = sol.omega(id);
    weightTrack(ik) = sol.m0Weight(id);
    selectedId(ik) = id;
end

band.k = kGrid;
band.omega = omegaTrack;
band.R = Rtrack;
band.L = Ltrack;
band.m0Weight = weightTrack;
band.selectedId = selectedId;
band.continuity = continuity;
end
