function [zakPhases, info] = compute_zak_phases(varargin)
%COMPUTE_ZAK_PHASES  Zak phases of the binary photonic time crystal.
%
%   ZAK = COMPUTE_ZAK_PHASES() evaluates the first seven bands of the PTC
%   used in Lustig, Sharabi, and Segev, Optica 5, 1390 (2018).  The
%   calculation implements the closed Floquet-frequency Wilson loop that
%   discretizes Eq. (5) of the paper,
%
%       theta_m = int_BZ dOmega i <u_m,Omega | d_Omega u_m,Omega>,
%
%   with <u|v> = int_0^T epsilon(t) u*(t)v(t) dt.  For every Omega, the
%   conserved momentum k_m(Omega) is found from the exact binary-crystal
%   dispersion relation.  The displacement field is then reconstructed
%   throughout the temporal unit cell from its [D;B] Floquet eigenstate.
%
%   The Omega grid is half-open.  Its final Wilson link therefore includes
%   the temporal-Bloch sewing transformation
%
%       u_{Omega+2*pi/T}(t) = exp(i*2*pi*t/T) u_Omega(t),
%
%   which is essential for a gauge-invariant closed loop.  This differs
%   from an open product of monodromy eigenvectors along k; the latter is
%   not the Zak phase in Eq. (5).
%
%   The Zak phase depends on the temporal origin.  The labels printed in
%   the published Fig. 1(b), [0 0 pi 0 0 pi 0], correspond to placing the
%   time-inversion centre in the epsilon=1 segment.  Moving the origin by
%   T/2 to the epsilon=3 centre adds pi to every band.  The paper's Fig. 1
%   and prose use inconsistent epsilon_1/epsilon_2 labels, so this origin
%   convention is stated explicitly here.
%
%   Name-value options:
%     'InversionCenter' - 'low' (paper labels) or 'high'
%     'NumberOfBands'   - number of complete bands to evaluate (default 7)
%     'NOmega'          - points on the half-open temporal BZ (default 61)
%     'NTime'           - samples in one temporal period (default 301)
%     'KMax'            - maximum k/k0 used to locate bands (default 5)
%     'Nk'              - band-search samples (default 20001)
%     'ComputeBothOrigins' - independently evaluate both inversion centres
%                            (default true)
%     'MakePlots'       - create diagnostic plots
%     'SaveOutput'      - save diagnostics under output/
%     'Verbose'         - print the numerical audit table
%
%   [ZAK,INFO] also returns raw Wilson phases, band edges, link-overlap
%   diagnostics, and the temporal-origin convention used.

defaultInteractive = (nargout == 0);
p = inputParser;
p.FunctionName = mfilename;
addParameter(p,'InversionCenter','low', ...
    @(x)ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p,'NumberOfBands',7, ...
    @(x)isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
addParameter(p,'NOmega',61, ...
    @(x)isnumeric(x) && isscalar(x) && x >= 15 && x == round(x));
addParameter(p,'NTime',301, ...
    @(x)isnumeric(x) && isscalar(x) && x >= 81 && x == round(x));
addParameter(p,'KMax',5, ...
    @(x)isnumeric(x) && isscalar(x) && x > 0);
addParameter(p,'Nk',20001, ...
    @(x)isnumeric(x) && isscalar(x) && x >= 1001 && x == round(x));
addParameter(p,'ComputeBothOrigins',true, ...
    @(x)islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p,'MakePlots',defaultInteractive, ...
    @(x)islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p,'SaveOutput',defaultInteractive, ...
    @(x)islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p,'Verbose',true, ...
    @(x)islogical(x) || (isnumeric(x) && isscalar(x)));
parse(p,varargin{:});
opt = p.Results;

centreChoice = validatestring(char(opt.InversionCenter),{'low','high'});

% Paper parameters.  c=1 is used by the repository TMM convention.  Only
% k*T is relevant, so T=2 is the paper's physical value in femtoseconds
% and k0=2*pi/T supplies the plotted dimensionless momentum k/k0.
epsHigh = 3;
epsLow = 1;
T = 2;
k0 = 2*pi/T;
reciprocalOmega = 2*pi/T;

if strcmp(centreChoice,'low')
    epsCentre = epsLow;
    epsOther = epsHigh;
else
    epsCentre = epsHigh;
    epsOther = epsLow;
end

% The trace is invariant under the half-period origin shift.  Use the
% vectorized exact binary dispersion to locate complete pass bands.
kNorm = linspace(0,opt.KMax,opt.Nk);
kValues = kNorm*k0;
halfTrace = binary_half_trace(kValues,epsHigh,epsLow,T/2,T/2);
bandMask = abs(halfTrace) <= 1 + 5e-11;
gapMask = ~bandMask;
[bandStarts,bandEnds] = logical_segments(bandMask);
[gapStarts,gapEnds] = logical_segments(gapMask);

nBands = min(opt.NumberOfBands,numel(bandStarts));
if nBands < opt.NumberOfBands
    error('Only %d bands were found below k/k0=%g; requested %d.', ...
        nBands,opt.KMax,opt.NumberOfBands);
end

bandEdges = zeros(nBands,2);
for bandId = 1:nBands
    bandEdges(bandId,1) = refine_band_edge( ...
        bandStarts(bandId),'left',kNorm,halfTrace,k0, ...
        epsHigh,epsLow,T);
    bandEdges(bandId,2) = refine_band_edge( ...
        bandEnds(bandId),'right',kNorm,halfTrace,k0, ...
        epsHigh,epsLow,T);
end

% Midpoint sampling avoids defective transfer matrices exactly at a band
% edge while still forming a closed BZ through the sewing link below.
nOmega = opt.NOmega;
omegaStep = reciprocalOmega/nOmega;
omegaGrid = -pi/T + ((0:nOmega-1)+0.5)*omegaStep;

nTime = opt.NTime;
timeStep = T/nTime;
timeGrid = (0:nTime-1)*timeStep;
epsTime = epsOther*ones(size(timeGrid));
epsTime(timeGrid < T/4 | timeGrid >= 3*T/4) = epsCentre;
sewingPhase = exp(1i*reciprocalOmega*timeGrid(:));

zakPhases = zeros(1,nBands);
rawZakPhases = zeros(1,nBands);
wilsonLoops = complex(zeros(1,nBands));
minLinkMagnitude = zeros(1,nBands);
quantizationError = zeros(1,nBands);
kOfOmega = zeros(nBands,nOmega);
allLinks = cell(1,nBands);

for bandId = 1:nBands
    modes = complex(zeros(nTime,nOmega));
    kLeft = bandEdges(bandId,1);
    kRight = bandEdges(bandId,2);

    for omegaId = 1:nOmega
        omega = omegaGrid(omegaId);
        targetTrace = cos(omega*T);
        kNormHere = solve_band_momentum(targetTrace,kLeft,kRight, ...
            k0,epsHigh,epsLow,T);
        kHere = kNormHere*k0;
        kOfOmega(bandId,omegaId) = kNormHere;

        U = centred_monodromy(kHere,epsCentre,epsOther,T);
        [vectors,multipliers] = eig(U,'vector');
        targetMultiplier = exp(-1i*omega*T);
        [~,modeId] = min(abs(multipliers-targetMultiplier));
        stateAtCentre = vectors(:,modeId);

        displacement = centred_displacement_mode( ...
            kHere,stateAtCentre,timeGrid,epsCentre,epsOther,T);
        periodicMode = displacement(:).*exp(1i*omega*timeGrid(:));

        normSquared = timeStep*sum(epsTime(:).*abs(periodicMode).^2);
        if ~isfinite(normSquared) || normSquared <= 100*eps
            error('Floquet-mode normalization failed for band %d.',bandId);
        end
        modes(:,omegaId) = periodicMode/sqrt(normSquared);
    end

    links = complex(zeros(1,nOmega));
    for omegaId = 1:nOmega-1
        links(omegaId) = timeStep*sum(epsTime(:).* ...
            conj(modes(:,omegaId)).*modes(:,omegaId+1));
    end

    % Close -pi/T -> +pi/T in the same temporal-Bloch basis.
    sewnFirstMode = sewingPhase.*modes(:,1);
    links(end) = timeStep*sum(epsTime(:).* ...
        conj(modes(:,end)).*sewnFirstMode);

    if any(abs(links) < 1e-10)
        error('A Wilson link is singular in band %d.',bandId);
    end
    unitLinks = links./abs(links);
    wilson = prod(unitLinks);
    rawZak = mod(-angle(wilson),2*pi);

    % Inversion symmetry quantizes the result.  Quantize only after storing
    % the raw Wilson phase and its distance from the nearest invariant.
    if real(wilson) >= 0
        quantizedZak = 0;
    else
        quantizedZak = pi;
    end

    zakPhases(bandId) = quantizedZak;
    rawZakPhases(bandId) = rawZak;
    wilsonLoops(bandId) = wilson;
    minLinkMagnitude(bandId) = min(abs(links));
    quantizationError(bandId) = abs(angle( ...
        wilson/exp(-1i*quantizedZak)));
    allLinks{bandId} = links;
end

publishedLowCentre = [0 0 pi 0 0 pi 0];
if strcmp(centreChoice,'low')
    expectedZak = publishedLowCentre(1:min(nBands,numel(publishedLowCentre)));
else
    expectedZak = mod(publishedLowCentre( ...
        1:min(nBands,numel(publishedLowCentre)))+pi,2*pi);
end
matchesExpected = numel(expectedZak) == nBands && ...
    all(abs(zakPhases-expectedZak) < 1e-12);

% Evaluate the second admissible inversion centre with an independent
% Wilson loop.  A T/2 shift of the temporal origin changes every single-
% band Zak phase by pi.  Keeping both results makes the convention issue
% in the paper directly auditable: the labels drawn in Fig. 1(b) select
% the low-epsilon centre, whereas a literal reading of segment 1 with
% epsilon_1=3 selects the high-epsilon centre.
zakLowCentre = [];
zakHighCentre = [];
rawZakLowCentre = [];
rawZakHighCentre = [];
otherCentreInfo = struct([]);
if logical(opt.ComputeBothOrigins)
    if strcmp(centreChoice,'low')
        otherCentre = 'high';
    else
        otherCentre = 'low';
    end
    [otherCentreZak,otherCentreInfo] = compute_zak_phases( ...
        'InversionCenter',otherCentre, ...
        'NumberOfBands',opt.NumberOfBands, ...
        'NOmega',opt.NOmega,'NTime',opt.NTime, ...
        'KMax',opt.KMax,'Nk',opt.Nk, ...
        'ComputeBothOrigins',false, ...
        'MakePlots',false,'SaveOutput',false,'Verbose',false);

    if strcmp(centreChoice,'low')
        zakLowCentre = zakPhases;
        rawZakLowCentre = rawZakPhases;
        zakHighCentre = otherCentreZak;
        rawZakHighCentre = otherCentreInfo.rawZakPhases;
    else
        zakLowCentre = otherCentreZak;
        rawZakLowCentre = otherCentreInfo.rawZakPhases;
        zakHighCentre = zakPhases;
        rawZakHighCentre = rawZakPhases;
    end
else
    if strcmp(centreChoice,'low')
        zakLowCentre = zakPhases;
        rawZakLowCentre = rawZakPhases;
    else
        zakHighCentre = zakPhases;
        rawZakHighCentre = rawZakPhases;
    end
end

if ~isempty(zakLowCentre) && ~isempty(zakHighCentre)
    originShiftResidual = angle(exp(1i*( ...
        zakHighCentre-zakLowCentre-pi)));
    originShiftCheck = all(abs(originShiftResidual) < 1e-12);
else
    originShiftResidual = [];
    originShiftCheck = [];
end

if opt.Verbose
    fprintf('=== Closed Omega-Wilson Zak calculation ===\n');
    fprintf('Inversion centre: epsilon = %g (%s-permittivity segment)\n', ...
        epsCentre,centreChoice);
    fprintf('Band | k/k0 interval          | raw Zak/pi | Zak | min |link|\n');
    fprintf('-----|------------------------|------------|-----|-----------\n');
    for bandId = 1:nBands
        fprintf('%4d | [%8.6f, %8.6f] | %10.6f | %3s | %.6f\n', ...
            bandId,bandEdges(bandId,1),bandEdges(bandId,2), ...
            rawZakPhases(bandId)/pi,zak_label(zakPhases(bandId)), ...
            minLinkMagnitude(bandId));
    end
    fprintf('Sequence: [%s]\n',strjoin(arrayfun( ...
        @zak_label,zakPhases,'UniformOutput',false),', '));
    if matchesExpected
        fprintf('PASS: sequence matches the expected %s-centred convention.\n', ...
            centreChoice);
    else
        warning('Computed Zak sequence does not match the expected convention.');
    end
    if logical(opt.ComputeBothOrigins)
        fprintf(['Published Fig. 1 labels (low-epsilon centre):  ', ...
            '[%s]\n'],strjoin(arrayfun(@zak_label,zakLowCentre, ...
            'UniformOutput',false),', '));
        fprintf(['Literal segment-1, epsilon_1=3 centre:          ', ...
            '[%s]\n'],strjoin(arrayfun(@zak_label,zakHighCentre, ...
            'UniformOutput',false),', '));
        if originShiftCheck
            fprintf('PASS: the two independently computed sequences differ by pi band-by-band.\n');
        else
            warning('The two origin conventions do not differ by pi band-by-band.');
        end
    end
end

info = struct();
info.method = 'closed Omega Wilson loop with temporal-Bloch sewing';
info.inversionCenter = centreChoice;
info.epsCentre = epsCentre;
info.epsOther = epsOther;
info.epsHigh = epsHigh;
info.epsLow = epsLow;
info.T = T;
info.k0 = k0;
info.kNorm = kNorm;
info.halfTrace = halfTrace;
info.bandMask = bandMask;
info.gapMask = gapMask;
info.bandStarts = bandStarts;
info.bandEnds = bandEnds;
info.gapStarts = gapStarts;
info.gapEnds = gapEnds;
info.bandEdges = bandEdges;
info.omegaGrid = omegaGrid;
info.kOfOmega = kOfOmega;
info.rawZakPhases = rawZakPhases;
info.zakPhases = zakPhases;
info.wilsonLoops = wilsonLoops;
info.links = allLinks;
info.minLinkMagnitude = minLinkMagnitude;
info.quantizationError = quantizationError;
info.expectedZak = expectedZak;
info.matchesExpected = matchesExpected;
info.computeBothOrigins = logical(opt.ComputeBothOrigins);
info.zakLowCentre = zakLowCentre;
info.zakHighCentre = zakHighCentre;
info.rawZakLowCentre = rawZakLowCentre;
info.rawZakHighCentre = rawZakHighCentre;
info.publishedLabelZak = zakLowCentre;
info.literalSegment1Zak = zakHighCentre;
info.originShiftResidual = originShiftResidual;
info.originShiftCheck = originShiftCheck;
if logical(opt.ComputeBothOrigins)
    info.otherCentreQuantizationError = ...
        otherCentreInfo.quantizationError;
    info.otherCentreMinLinkMagnitude = ...
        otherCentreInfo.minLinkMagnitude;
else
    info.otherCentreQuantizationError = [];
    info.otherCentreMinLinkMagnitude = [];
end

fig = [];
if logical(opt.MakePlots)
    fig = make_diagnostic_figure(kNorm,halfTrace, ...
        bandStarts,bandEnds,gapStarts,gapEnds,zakPhases, ...
        rawZakPhases,T);
end

if logical(opt.SaveOutput)
    outputDir = fullfile(fileparts(mfilename('fullpath')),'output');
    if ~exist(outputDir,'dir'), mkdir(outputDir); end
    save(fullfile(outputDir,'compute_zak_phases.mat'),'zakPhases','info');
    if ~isempty(fig)
        outputFile = fullfile(outputDir,'compute_zak_phases.png');
        try
            exportgraphics(fig,outputFile,'Resolution',220);
        catch
            print(fig,outputFile,'-dpng','-r220');
        end
        if opt.Verbose, fprintf('Saved: %s\n',outputFile); end
    end
end
end

function h = binary_half_trace(k,epsA,epsB,durationA,durationB)
nA = sqrt(epsA);
nB = sqrt(epsB);
phaseA = k*durationA/nA;
phaseB = k*durationB/nB;
h = cos(phaseA).*cos(phaseB) - 0.5*(nA/nB+nB/nA).* ...
    sin(phaseA).*sin(phaseB);
end

function [starts,ends] = logical_segments(mask)
d = diff([false mask false]);
starts = find(d == 1);
ends = find(d == -1)-1;
end

function edge = refine_band_edge(index,side,kNorm,halfTrace,k0, ...
    epsHigh,epsLow,T)
switch side
    case 'left'
        if index == 1
            edge = kNorm(1);
            return;
        end
        bracket = [kNorm(index-1) kNorm(index)];
        targetSign = sign(halfTrace(index-1));
    case 'right'
        if index == numel(kNorm)
            edge = kNorm(end);
            return;
        end
        bracket = [kNorm(index) kNorm(index+1)];
        targetSign = sign(halfTrace(index+1));
    otherwise
        error('Unknown edge side.');
end
fun = @(q) binary_half_trace(q*k0,epsHigh,epsLow,T/2,T/2) ...
    - targetSign;
edge = fzero(fun,bracket);
end

function kNorm = solve_band_momentum(targetTrace,kLeft,kRight,k0, ...
    epsHigh,epsLow,T)
fun = @(q) binary_half_trace(q*k0,epsHigh,epsLow,T/2,T/2) ...
    - targetTrace;
leftValue = fun(kLeft);
rightValue = fun(kRight);
rootTolerance = 5e-12;
if abs(leftValue) < rootTolerance
    kNorm = kLeft;
elseif abs(rightValue) < rootTolerance
    kNorm = kRight;
elseif leftValue*rightValue > 0
    error('Band-momentum root is not bracketed in [%g,%g].',kLeft,kRight);
else
    kNorm = fzero(fun,[kLeft kRight]);
end
end

function U = centred_monodromy(k,epsCentre,epsOther,T)
halfCentre = segment_propagator(k,epsCentre,T/4);
other = segment_propagator(k,epsOther,T/2);
U = halfCentre*other*halfCentre;
end

function P = segment_propagator(k,epsr,duration)
n = sqrt(epsr);
phase = k*duration/n;
c = cos(phase);
s = sin(phase);
P = [c, -1i*n*s; -1i*s/n, c];
end

function displacement = centred_displacement_mode( ...
    k,stateAtCentre,timeGrid,epsCentre,epsOther,T)
displacement = complex(zeros(size(timeGrid)));

mask1 = timeGrid < T/4;
mask2 = timeGrid >= T/4 & timeGrid < 3*T/4;
mask3 = timeGrid >= 3*T/4;

[displacement(mask1),~] = advance_state( ...
    stateAtCentre,k,epsCentre,timeGrid(mask1));

stateQuarter = segment_propagator(k,epsCentre,T/4)*stateAtCentre;
[displacement(mask2),~] = advance_state( ...
    stateQuarter,k,epsOther,timeGrid(mask2)-T/4);

stateThreeQuarter = segment_propagator(k,epsOther,T/2)*stateQuarter;
[displacement(mask3),~] = advance_state( ...
    stateThreeQuarter,k,epsCentre,timeGrid(mask3)-3*T/4);
end

function [D,B] = advance_state(state,k,epsr,duration)
n = sqrt(epsr);
phase = k*duration/n;
c = cos(phase);
s = sin(phase);
D = c*state(1) - 1i*n*s*state(2);
B = c*state(2) - 1i*s*state(1)/n;
end

function label = zak_label(zak)
if abs(zak) < 0.25*pi
    label = '0';
else
    label = 'pi';
end
end

function fig = make_diagnostic_figure(kNorm,halfTrace, ...
    bandStarts,bandEnds,gapStarts,gapEnds,zakPhases,rawZak,T)
clamped = min(1,max(-1,halfTrace));
omegaT = acos(clamped);

fig = figure('Color','w','Position',[80 80 1050 430]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl);
hold(ax1,'on');
for gapId = 1:numel(gapStarts)
    x1 = kNorm(gapStarts(gapId));
    x2 = kNorm(gapEnds(gapId));
    patch(ax1,[x1 x2 x2 x1],[-pi -pi pi pi],[0.86 0.86 0.86], ...
        'EdgeColor','none');
end
for bandId = 1:numel(bandStarts)
    ids = bandStarts(bandId):bandEnds(bandId);
    plot(ax1,kNorm(ids), omegaT(ids),'b-','LineWidth',1.4);
    plot(ax1,kNorm(ids),-omegaT(ids),'b-','LineWidth',1.4);
end
for bandId = 1:min(numel(zakPhases),numel(bandStarts))
    x = 0.5*(kNorm(bandStarts(bandId))+kNorm(bandEnds(bandId)));
    text(ax1,x,0,zak_label(zakPhases(bandId)), ...
        'HorizontalAlignment','center','BackgroundColor','w');
end
xlabel(ax1,'k/k_0');
ylabel(ax1,'\Omega T');
xlim(ax1,[kNorm(1) kNorm(end)]);
ylim(ax1,[-pi pi]);
title(ax1,'Closed-\Omega Zak labels');
box(ax1,'on');

ax2 = nexttile(tl);
bar(ax2,1:numel(zakPhases),zakPhases/pi,0.55, ...
    'FaceColor',[0.20 0.45 0.75]);
hold(ax2,'on');
plot(ax2,1:numel(rawZak),rawZak/pi,'ro','MarkerFaceColor','r');
yline(ax2,0,'k:');
yline(ax2,1,'k:');
xlabel(ax2,'Band index');
ylabel(ax2,'Zak phase / \pi');
ylim(ax2,[-0.1 1.1]);
legend(ax2,'quantized','raw Wilson phase','Location','best');
title(ax2,sprintf('Temporal BZ: [-\\pi/T,\\pi/T), T=%g',T));
box(ax2,'on');
end
