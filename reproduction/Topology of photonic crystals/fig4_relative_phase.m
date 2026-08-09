function fig4_relative_phase()
%FIG4_RELATIVE_PHASE Physically constrained reproduction of paper Fig. 4.
%
% Fig. 4 in Lustig et al. contains phases extracted from FDTD Fourier
% components.  The main paper does not specify the probe position, Fourier
% gate/reference time, carrier spectra for gaps 2--6, or the precise phase
% of the switch-off boundary.  Those choices change the phase curve while
% Eq. (6) fixes only its sign in a gap.  Consequently the six published
% curves cannot be recovered uniquely from the stated material parameters.
%
% This script therefore plots an explicit analytic proxy, rather than
% silently fitting the published curves.  At every k it follows the Floquet
% branch which grows inside the corresponding gap, uses the inversion-
% symmetric eps1-centred unit cell, and evaluates
%
%                  phi = arg(-E_minus/E_plus)
%
% in a fixed free-space directional basis at the symmetry plane.  In a
% strict momentum gap this is the asymptotic (long-crystal) finite-TMM
% phase in this convention.  A growing branch is not uniquely defined in
% the portions of the published windows outside a strict gap, so those
% portions are deliberately left unplotted.  The robust result is the gap-sign sequence
% (+,+,-,-,-,+), not the convention-dependent ordinate of every point.

rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(rootDir,'startup_stm.m'));

eps1 = 3;
eps2 = 1;
T = 2;                              % fs
k0 = 2*pi/T;                        % c=1; horizontal axis is k/k0
epsCell = [eps1 eps2 eps1];
muCell = ones(size(epsCell));
durations = [T/4 T/2 T/4];          % inversion centre in eps1
epsReference = 1;
muReference = 1;

% Axes read from the published Fig. 4.  They are measurement windows, not
% the analytical band-edge intervals (notably for panels c and f).
paperWindows = [ ...
    0.560 0.710; ...
    1.210 1.350; ...
    1.840 1.980; ...
    2.545 2.574; ...
    3.198 3.228; ...
    3.852 3.882];
publishedSigns = [1 1 -1 -1 -1 1];
nGapsToShow = size(paperWindows,1);

% Locate strict gaps independently of the display windows.
kNormScan = linspace(1e-6,4.05,40001);
kScan = kNormScan*k0;
bands = temporal_crystal_bands(kScan,epsCell,muCell,durations);
halfTrace = real(bands.halfTrace);
inGap = abs(halfTrace) > 1;
[gapStarts,gapEnds] = logical_segments(inGap);
if numel(gapStarts) < nGapsToShow
    error('The k scan found only %d gaps; six are required.',numel(gapStarts));
end
strictGapBounds = [kNormScan(gapStarts(1:nGapsToShow)).', ...
    kNormScan(gapEnds(1:nGapsToShow)).'];

kWindowData = cell(1,nGapsToShow);
phiData = cell(1,nGapsToShow);
strictGapMask = cell(1,nGapsToShow);
floquetMultiplier = cell(1,nGapsToShow);
phaseSigns = zeros(1,nGapsToShow);

fprintf(['Fig. 4: growing-Floquet-mode phase proxy.\n' ...
    '  The paper does not report enough FDTD phase-reference details ' ...
    'to identify the plotted ordinates uniquely.\n']);
for g = 1:nGapsToShow
    kNorm = linspace(paperWindows(g,1),paperWindows(g,2),601);
    kValues = kNorm*k0;
    gapMid = mean(strictGapBounds(g,:))*k0;
    [states,multipliers] = tracked_growing_branch(kValues, ...
        epsCell,muCell,durations,gapMid);

    directional = temporal_db_to_directional(states, ...
        epsReference,muReference);
    ratio = -directional(2,:)./directional(1,:);
    phi = unwrap(angle(ratio));

    localHalfTrace = zeros(size(kValues));
    for ik = 1:numel(kValues)
        U = temporal_crystal_monodromy(kValues(ik),epsCell, ...
            muCell,durations);
        localHalfTrace(ik) = real(trace(U)/2);
    end
    isStrictGap = abs(localHalfTrace) > 1;
    isInterior = abs(localHalfTrace) > 1+1e-5;
    if ~any(isInterior)
        error('Published window %d does not intersect its strict gap.',g);
    end

    % Select the continuous 2*pi branch whose gap interior has the Eq. (6)
    % sign.  This is a branch choice, not a fit to digitized figure data.
    target = publishedSigns(g)*pi/2;
    phi = phi + 2*pi*round((target-median(phi(isInterior)))/(2*pi));

    phaseSigns(g) = sign(median(phi(isInterior)));
    if phaseSigns(g) ~= publishedSigns(g)
        error('Gap %d phase sign is %+d; expected %+d.',g, ...
            phaseSigns(g),publishedSigns(g));
    end
    if any(publishedSigns(g)*phi(isInterior) < -2e-3)
        error('Gap %d proxy changes topological sign in its interior.',g);
    end
    strictPhi = phi(isStrictGap);
    if any(abs(diff(strictPhi)) > pi)
        error('Gap %d contains an unresolved in-gap 2*pi phase wrap.',g);
    end

    kWindowData{g} = kNorm;
    phiData{g} = phi;
    strictGapMask{g} = isStrictGap;
    floquetMultiplier{g} = multipliers;
    fprintf(['  gap %d: strict [%.6f, %.6f], displayed [%.3f, %.3f], ' ...
        'median(phi/pi)=%+.4f, sign=%+d\n'],g, ...
        strictGapBounds(g,1),strictGapBounds(g,2), ...
        paperWindows(g,1),paperWindows(g,2), ...
        median(phi(isInterior))/pi,phaseSigns(g));
end

fig = figure('Color','w','Position',[40 40 960 590]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
for g = 1:nGapsToShow
    ax = nexttile(tl);
    hold(ax,'on');
    kNorm = kWindowData{g};
    phiPi = phiData{g}/pi;
    isStrictGap = strictGapMask{g};

    [segmentStarts,segmentEnds] = logical_segments(isStrictGap);
    for s = 1:numel(segmentStarts)
        idx = segmentStarts(s):segmentEnds(s);
        plot(ax,kNorm(idx),phiPi(idx),'Color',[0 0.4470 0.7410], ...
            'LineWidth',2.0);
    end
    yline(ax,0,'Color',[0.60 0.60 0.60],'LineWidth',0.6);
    xlim(ax,paperWindows(g,:));
    visiblePhi = phiPi(isStrictGap);
    yPad = max(0.08,0.08*(max(visiblePhi)-min(visiblePhi)));
    ylim(ax,[min(-0.05,min(visiblePhi)-yPad), ...
        max(0.05,max(visiblePhi)+yPad)]);
    title(ax,sprintf('(%c)   Gap - %d',char('a'+g-1),g), ...
        'FontWeight','normal');
    if g > 3, xlabel(ax,'k/k_0'); end
    if mod(g-1,3) == 0, ylabel(ax,'Phase [1/\pi]'); end
    box(ax,'on');
    set(ax,'FontSize',9,'Layer','top');
    if isprop(ax,'Toolbar')
        ax.Toolbar.Visible = 'off';
    end
    disableDefaultInteractivity(ax);
end
title(tl,'Fig. 4 analytic proxy: in-gap growing Floquet branch', ...
    'FontWeight','normal','FontSize',11);

outputDir = fullfile(fileparts(mfilename('fullpath')),'output');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
outputFile = fullfile(outputDir,'fig4_relative_phase.png');
exportgraphics(fig,outputFile,'Resolution',250);
phaseDefinition = ['phi=arg(-E_minus/E_plus) at the eps1-centred ' ...
    'symmetry plane; growing Floquet branch; free-space basis'];
limitation = ['Published FDTD probe, Fourier gate/reference time, input ' ...
    'spectra, and exact termination phase are not specified; Eq. (6) ' ...
    'determines the in-gap sign but not a unique pointwise phase curve.'];
save(fullfile(outputDir,'fig4_data.mat'), ...
    'kNormScan','kScan','halfTrace','inGap','gapStarts','gapEnds', ...
    'strictGapBounds','paperWindows','kWindowData','phiData', ...
    'strictGapMask','floquetMultiplier','phaseSigns','publishedSigns', ...
    'phaseDefinition','limitation','eps1','eps2','T','k0', ...
    'epsCell','durations');
fprintf('Saved: %s\n',outputFile);
end

function [states,multipliers] = tracked_growing_branch(kValues, ...
    epsCell,muCell,durations,gapMid)
% Anchor on the growing gap eigenvector, then continue it by overlap.
nK = numel(kValues);
states = complex(zeros(2,nK));
multipliers = complex(zeros(1,nK));
[~,anchor] = min(abs(kValues-gapMid));

[V,D] = eig(temporal_crystal_monodromy(kValues(anchor), ...
    epsCell,muCell,durations));
lambda = diag(D);
[~,branch] = max(abs(lambda));
states(:,anchor) = V(:,branch)/norm(V(:,branch));
multipliers(anchor) = lambda(branch);

for ik = anchor+1:nK
    [V,D] = eig(temporal_crystal_monodromy(kValues(ik), ...
        epsCell,muCell,durations));
    V = V./vecnorm(V);
    [~,branch] = max(abs(states(:,ik-1)'*V));
    states(:,ik) = V(:,branch);
    multipliers(ik) = D(branch,branch);
end
for ik = anchor-1:-1:1
    [V,D] = eig(temporal_crystal_monodromy(kValues(ik), ...
        epsCell,muCell,durations));
    V = V./vecnorm(V);
    [~,branch] = max(abs(states(:,ik+1)'*V));
    states(:,ik) = V(:,branch);
    multipliers(ik) = D(branch,branch);
end
end

function [starts,ends] = logical_segments(mask)
d = diff([false mask false]);
starts = find(d == 1);
ends = find(d == -1)-1;
end
