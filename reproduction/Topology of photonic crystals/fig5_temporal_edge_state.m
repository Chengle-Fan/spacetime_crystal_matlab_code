function fig5_temporal_edge_state()
%FIG5_TEMPORAL_EDGE_STATE Reproduce Fig. 5 of Lustig et al. (2018).
%
% A narrow finite-bandwidth pulse is propagated in k space with the exact
% D/B propagator.  Keeping the full spectrum is essential: an exactly
% matched single-k domain-wall eigenstate decays forever to the right and
% cannot reproduce the re-growth shown in the paper.  The paper does not
% state this panel's pulse width, so the value is explicitly inferred below.

rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(rootDir,'startup_stm.m'));

epsA = 3;
epsB = 1;
T = 2;                              % fs
c0 = 0.299792458;                   % um/fs
k0 = 2*pi/T;                        % c=1 numerical coordinate
nLeft = 8;
nRight = 8;
interfacePeriod = nLeft;

% Time-inversion-symmetric cells: segment 1 is centred at the cell origin.
epsLeft = [epsA epsB epsA];
epsRight = [epsB epsA epsB];
muCell = [1 1 1];
durations = [T/4 T/2 T/4];

% Diagnose the ideal (single-k) domain-wall state, but do not use its norm
% as the pulse observable.
kNormProbe = linspace(0.53,0.73,1600);
mode = temporal_domain_wall_mode(kNormProbe*k0, ...
    epsLeft,muCell,durations,epsRight,muCell,durations,nLeft,nRight);

% Fig. 5 does not state a pulse wavelength or width (the 0.93 um / 45 fs
% values belong only to Fig. 2).  The physically distinguished centre is
% therefore the domain-wall matching mode found above.  The width below is
% inferred from the published finite-bandwidth envelope: after scaling the
% interface peak to 60, the last-period peak is approximately 55.
kCentre = mode.k;
kCentreNorm = kCentre/k0;
lambdaCentre = 2*pi*c0/kCentre;     % 0.987294 um
fwhmIntensity = 189;                % fs, inferred from published panel (b)
sigmaX = fwhmIntensity/sqrt(2*log(2));

% Uniform baseband grid for FFT reconstruction of max_z |D(z,t)|.
nK = 2048;
qMax = 0.8;
dq = 2*qMax/nK;
q = (-nK/2:nK/2-1)*dq;
kSpectrum = kCentre+q;
spectralAmplitude = exp(-((q*sigmaX/2).^2));

stepsPerPeriod = 200;               % all quarter-cell interfaces align
[timeFs,peakAmplitude] = wavepacket_peak_history( ...
    kSpectrum,spectralAmplitude,T,stepsPerPeriod, ...
    nLeft,nRight,epsA,epsB);
timePeriods = timeFs/T;

% Compare the carrier-resolved curve through a one-period local-peak
% envelope; sampling only at integer T can accidentally hit a carrier dip.
cycleTimes = 0:(nLeft+nRight);
cycleEnvelope = zeros(size(cycleTimes));
for j = 1:numel(cycleTimes)
    window = abs(timePeriods-cycleTimes(j)) <= 0.5;
    cycleEnvelope(j) = max(peakAmplitude(window));
end
idInterfaceCycle = interfacePeriod+1;
amplitudeScale = 60/cycleEnvelope(idInterfaceCycle);
peakAmplitude = amplitudeScale*peakAmplitude;
cycleEnvelope = amplitudeScale*cycleEnvelope;
interfacePeak = cycleEnvelope(idInterfaceCycle);
rightMinimum = min(cycleEnvelope(idInterfaceCycle:end));
finalAmplitude = cycleEnvelope(end);
fprintf(['Fig. 5 wavepacket: k_c/k0=%.6f, ideal edge k/k0=%.6f, ' ...
    'cycle peak at 8T=%.3f, right minimum=%.3f, final=%.3f.\n'], ...
    kCentreNorm,mode.k/k0,interfacePeak,rightMinimum,finalAmplitude);
if abs(interfacePeak-60) > 1e-10 || rightMinimum >= 0.35*interfacePeak || ...
        finalAmplitude < 0.75*interfacePeak || ...
        finalAmplitude > 1.05*interfacePeak
    error('The reconstructed pulse does not match the published edge envelope.');
end

fig = figure('Color','w','Position',[40 40 920 560]);
tl = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

% ---------------------------------------------------------------------
% (a) Two topologically distinct, cascaded centred cells.
% ---------------------------------------------------------------------
axA = nexttile(tl);
uStep = linspace(0,nLeft+nRight,6401);
epsStep = domain_wall_epsilon(uStep,nLeft,epsA,epsB);
stairs(axA,uStep,epsStep,'Color',[0 0.4470 0.7410],'LineWidth',1.4);
xline(axA,interfacePeriod,':','Color',[0.35 0.35 0.35]);
%text(axA,interfacePeriod+0.15,2.65,'interface','FontSize',8);
xlim(axA,[0 nLeft+nRight]);
ylim(axA,[0.75 3.25]);
yticks(axA,[epsB epsA]);
ylabel(axA,'\epsilon');
title(axA,'(a)','FontWeight','normal');
box(axA,'on');
set(axA,'FontSize',9);

% ---------------------------------------------------------------------
% (b) Actual pulse peak, not ||[D,B]|| of a pure eigenstate.
% ---------------------------------------------------------------------
axB = nexttile(tl);
plot(axB,timePeriods,peakAmplitude,'Color',[0 0.4470 0.7410], ...
    'LineWidth',1.5);
xline(axB,interfacePeriod,':','Color',[0.35 0.35 0.35]);
xlim(axB,[0 nLeft+nRight]);
ylim(axB,[0 1.06*max(peakAmplitude)]);
ylabel(axB,'|D|');
xlabel(axB,'Modulation Periods');
title(axB,'(b)','FontWeight','normal');
box(axB,'on');
set(axB,'FontSize',9);

% ---------------------------------------------------------------------
% (c) Periodic, inversion-symmetric smooth modulation.  tanh(cos) is
% continuous at both the half-period and period boundaries.
% ---------------------------------------------------------------------
axC = nexttile(tl);
uSmooth = linspace(0,9,2401);
beta = 2.4;
epsSmooth = 0.5*(epsA+epsB) + 0.5*(epsA-epsB)* ...
    tanh(beta*cos(2*pi*uSmooth))/tanh(beta);
plot(axC,uSmooth,epsSmooth,'Color',[0 0.4470 0.7410],'LineWidth',1.5);
xlim(axC,[0 9]);
ylim(axC,[0.75 3.25]);
yticks(axC,[epsB epsA]);
xlabel(axC,'Modulation Periods');
ylabel(axC,'\epsilon');
title(axC,'(c)','FontWeight','normal');
box(axC,'on');
set(axC,'FontSize',9);

outputDir = fullfile(fileparts(mfilename('fullpath')),'output');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
outputFile = fullfile(outputDir,'fig5_temporal_edge_state.png');
exportgraphics(fig,outputFile,'Resolution',250);
save(fullfile(outputDir,'fig5_data.mat'), ...
    'timeFs','timePeriods','peakAmplitude','interfacePeak','rightMinimum', ...
    'finalAmplitude','cycleTimes','cycleEnvelope','amplitudeScale', ...
    'uStep','epsStep','uSmooth','epsSmooth','mode', ...
    'kCentre','kCentreNorm','kSpectrum','spectralAmplitude', ...
    'lambdaCentre','fwhmIntensity','sigmaX','epsA','epsB','T', ...
    'epsLeft','epsRight','durations','nLeft','nRight');
fprintf('Saved: %s\n',outputFile);
end

function [timeFs,peakAmplitude] = wavepacket_peak_history( ...
    kSpectrum,spectrum,T,stepsPerPeriod,nLeft,nRight,epsA,epsB)
dt = T/stepsPerPeriod;
nSteps = (nLeft+nRight)*stepsPerPeriod;
timeFs = (0:nSteps)*dt;

D = complex(spectrum);
% At t=0 the left crystal is in epsA.  For its forward plane-wave
% components, D=epsA*E and B=H=sqrt(epsA)*E, hence B/D=1/sqrt(epsA).
B = complex(spectrum)/sqrt(epsA);
normalization = max(abs(ifft(ifftshift(D))));
peakAmplitude = zeros(size(timeFs));
peakAmplitude(1) = 1;

for step = 1:nSteps
    tMid = (step-0.5)*dt;
    uMid = tMid/T;
    epsNow = domain_wall_epsilon(uMid,nLeft,epsA,epsB);
    rootEps = sqrt(epsNow);
    phase = kSpectrum*dt/rootEps;
    c = cos(phase);
    s = sin(phase);
    Dold = D;
    Bold = B;
    D = c.*Dold-1i*rootEps*s.*Bold;
    B = -1i*s/rootEps.*Dold+c.*Bold;
    peakAmplitude(step+1) = ...
        max(abs(ifft(ifftshift(D))))/normalization;
end
end

function epsValue = domain_wall_epsilon(u,nLeft,epsA,epsB)
% epsA is centred at integer periods on the left; epsB is centred on the
% right.  This produces the required temporal phase slip at u=nLeft.
isLeft = u < nLeft;
highAtInteger = cos(2*pi*u) >= 0;
epsValue = epsB*ones(size(u));
epsValue(isLeft & highAtInteger) = epsA;
epsValue(~isLeft & ~highAtInteger) = epsA;
end
