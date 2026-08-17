function fig1_ptc_bands()
%FIG1_PTC_BANDS Reproduce Fig. 1 of Lustig et al., Optica 5, 1390 (2018).
%
% The paper plots the extended temporal-Brillouin-zone ordinate Omega*T,
% not only the two principal quasifrequency branches.  Momentum is shown as
% k/k0, with k0 = 2*pi/(T*c).  Zak phases require an explicit temporal
% origin.  The labels in the published panel correspond to putting t=0 in
% the middle of the low-permittivity segment, so the symmetric unit cell is
% [eps_low(T/4), eps_high(T/2), eps_low(T/4)].  A literal reading of the
% prose definition (segment 1 has epsilon_1=3 at t=0) instead chooses the
% high-permittivity centre and shifts every one-band Zak phase by pi.

rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(rootDir, 'startup_stm.m'));
addpath(fileparts(mfilename('fullpath')));

eps1 = 3;
eps2 = 1;
mu = 1;
T = 2;                         % fs (the normalization is immaterial here)
k0 = 2*pi/T;                   % c=1 in the TMM coordinate
epsCell = [eps2 eps1 eps2];
muCell = mu*ones(size(epsCell));
durations = [T/4 T/2 T/4];

% The published panel extends to k/k0=5 and displays seven Zak labels.
kNorm = linspace(0, 5, 6001);
k = kNorm*k0;
bands = temporal_crystal_bands(k, epsCell, muCell, durations);
halfTrace = real(bands.halfTrace);
inBand = abs(halfTrace) <= 1 + 5e-11;
inGap = ~inBand;

[bandStarts, bandEnds] = logical_segments(inBand);
[gapStarts, gapEnds] = logical_segments(inGap);
nBands = numel(bandStarts);
nGaps = numel(gapStarts);

% In a pass band, cos(Omega*T)=Tr(U)/2.  acos therefore gives the two
% extended-zone branches +/-Omega*T in [-pi,pi].
clampedTrace = min(1, max(-1, halfTrace));
omegaT = acos(clampedTrace);

% Compute, rather than assign, the topology.  COMPUTE_ZAK_PHASES evaluates
% the closed Omega-Wilson loop of Eq. (5), reconstructing the displacement
% mode throughout the temporal cell and applying the BZ sewing phase.
nZakLabels = min(7,nBands);
[zakPhases,zakInfo] = compute_zak_phases( ...
    'InversionCenter','low', ...
    'NumberOfBands',nZakLabels, ...
    'NOmega',61,'NTime',301, ...
    'KMax',5,'Nk',20001, ...
    'ComputeBothOrigins',true, ...
    'MakePlots',false,'SaveOutput',false,'Verbose',true);

% Paper values are used only as an acceptance test, never as plot data.
% Test both temporal-origin conventions so the epsilon-label ambiguity in
% panel (a) cannot silently change which topology is displayed.
zakReference = [0 0 pi 0 0 pi 0];
zakHighCentreReference = mod(zakReference+pi,2*pi);
if any(abs(zakPhases-zakReference(1:nZakLabels)) > 1e-12)
    error('Computed Zak sequence does not match the published Fig. 1 labels.');
end
if any(abs(zakInfo.zakLowCentre-zakReference(1:nZakLabels)) > 1e-12)
    error('Low-epsilon-centred Zak audit failed.');
end
if any(abs(zakInfo.zakHighCentre- ...
        zakHighCentreReference(1:nZakLabels)) > 1e-12)
    error('High-epsilon-centred Zak audit failed.');
end
if ~zakInfo.originShiftCheck
    error('The two inversion-origin Zak sequences do not differ by pi.');
end

fprintf('Fig. 1: %d bands and %d gaps found for 0 <= k/k0 <= 5.\n', ...
    nBands, nGaps);
for g = 1:min(7,nGaps)
    fprintf('  gap %d: k/k0 = [%.6f, %.6f]\n', g, ...
        kNorm(gapStarts(g)), kNorm(gapEnds(g)));
end

fig = figure('Color','w','Position',[60 60 1050 430]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

% ---------------------------------------------------------------------
% (a) Binary PTC.  The low-permittivity segment is centred at integer
% modulation periods, matching the temporal origin used for the Zak phase.
% ---------------------------------------------------------------------
axA = nexttile(tl);
u = linspace(-0.5,10,4201);    % u=t/T
epsProfile = eps1*ones(size(u));
epsProfile(cos(2*pi*u) >= 0) = eps2;
stairs(axA,u,epsProfile,'Color',[0 0.4470 0.7410],'LineWidth',1.7);
hold(axA,'on');
xline(axA,0,':','Color',[0.45 0.45 0.45],'LineWidth',0.8);
% text(axA,0.12,1.22,'time-inversion centre','FontSize',8, ...
%     'Color',[0.35 0.35 0.35]);
xlabel(axA,'t/T');
ylabel(axA,'\epsilon(t)');
yticks(axA,[eps2 eps1]);
yticklabels(axA,{'\epsilon_2','\epsilon_1'});
xlim(axA,[-0.5 10]);
ylim(axA,[0.75 3.25]);
title(axA,'(a)','FontWeight','normal');
box(axA,'on');
set(axA,'FontSize',10);

% ---------------------------------------------------------------------
% (b) Extended-zone dispersion.  Plot each connected pass band
% separately so MATLAB cannot draw spurious blue chords across gaps.
% ---------------------------------------------------------------------
axB = nexttile(tl);
hold(axB,'on');
for g = 1:nGaps
    x1 = kNorm(gapStarts(g));
    x2 = kNorm(gapEnds(g));
    patch(axB,[x1 x2 x2 x1],[-4 -4 4 4],[0.84 0.84 0.84], ...
        'EdgeColor','none','FaceAlpha',0.65);
end
for b = 1:nBands
    ids = bandStarts(b):bandEnds(b);
    plot(axB,kNorm(ids), omegaT(ids),'Color',[0.12 0.32 0.68], ...
        'LineWidth',1.8);
    plot(axB,kNorm(ids),-omegaT(ids),'Color',[0.12 0.32 0.68], ...
        'LineWidth',1.8);
end
for b = 1:nZakLabels
    km = 0.5*(kNorm(bandStarts(b))+kNorm(bandEnds(b)));
    if zakPhases(b) == 0
        label = '0';
    else
        label = '\pi';
    end
    text(axB,km,0,label,'HorizontalAlignment','center', ...
        'VerticalAlignment','middle','FontSize',10, ...
        'BackgroundColor','w','Margin',1);
end
xlabel(axB,'k/k_0');
ylabel(axB,'\Omega T');
xlim(axB,[0 5]);
ylim(axB,[-4 4]);
yticks(axB,-4:2:4);
title(axB,'(b)','FontWeight','normal');
box(axB,'on');
set(axB,'FontSize',10,'Layer','top');

outputDir = fullfile(fileparts(mfilename('fullpath')),'output');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
outputFile = fullfile(outputDir,'fig1_ptc_bands.png');
exportgraphics(fig,outputFile,'Resolution',250);
save(fullfile(outputDir,'fig1_data.mat'), ...
    'kNorm','k','k0','omegaT','halfTrace','inBand','inGap', ...
    'bandStarts','bandEnds','gapStarts','gapEnds', ...
    'zakPhases','zakInfo','zakReference','zakHighCentreReference', ...
    'eps1','eps2','T','epsCell','durations');
fprintf('Saved: %s\n',outputFile);
end

function [starts,ends] = logical_segments(mask)
d = diff([false mask false]);
starts = find(d == 1);
ends = find(d == -1)-1;
end
