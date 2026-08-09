function fig2_fdtd_simulations(forceRecompute)
%FIG2_FDTD_SIMULATIONS Reproduce Fig. 2 of Lustig et al. (2018).
%
% The FDTD coordinate is x=z/c0 in femtoseconds, so fdtd1d_db can retain
% c=1 while every physical time in the paper is entered directly in fs.
% The saved and plotted field is the electric displacement D, not E.

if nargin < 1
    forceRecompute = false;
end
validateattributes(forceRecompute,{'logical','numeric'},{'scalar'});
forceRecompute = logical(forceRecompute);

rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(rootDir,'startup_stm.m'));

outputDir = fullfile(fileparts(mfilename('fullpath')),'output');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

% -------------------------------------------------------------------------
% Literal parameters stated in the article.
% -------------------------------------------------------------------------
p.c0 = 0.299792458;                 % um/fs
p.eps1 = 3;
p.eps2 = 1;
p.epsBackground = 1;               % free space
p.mu = 1;
p.T = 2;                            % fs
p.tStart = 220;                     % fs
p.tEnd = 340;                       % fs (60 periods)
p.tStop = 500;                      % fs
p.nPeriods = round((p.tEnd-p.tStart)/p.T);
p.fwhm = 45;                        % fs, FWHM of the plotted |D| envelope
p.sigmaX = p.fwhm/(2*sqrt(log(2))); % exp[-(x/sigmaX)^2] amplitude
p.zInitial = 0;                     % um
p.zMin = -50;                       % padded computational domain
p.zMax = 190;
p.zPlot = [0 140];
p.dt = p.T/80;                      % 40 steps per half-period
p.courant = 0.8;
p.dx = p.dt/p.courant;              % x coordinate has units fs
p.recordEvery = 20;                 % output every 0.5 fs
p.modelVersion = 'fig2-literal60-D-boundary-v10-converged';

if abs(p.nPeriods*p.T-(p.tEnd-p.tStart)) > 10*eps(p.tEnd)
    error('The PTC window is not an integer number of periods.');
end

caseBand = struct('name','inband','label','Band propagation', ...
    'lambda',1.4,'clim',[-1.5 0]);
caseGap = struct('name','ingap','label','Gap propagation', ...
    'lambda',0.93,'clim',[0 10]);

dataBand = run_or_load_case(caseBand,p,outputDir,forceRecompute);
dataGap = run_or_load_case(caseGap,p,outputDir,forceRecompute);

% -------------------------------------------------------------------------
% The article states both N=60 and an approximately 20,000-fold maximum.
% Exact TMM exposes that these two numbers are mutually inconsistent.
% -------------------------------------------------------------------------
kGap = 2*pi*p.c0/caseGap.lambda;
% The finite Fig. 2 crystal turns on at the beginning of an eps1 segment.
% This phase gives the four published in-band branches.  It is a cyclic
% representation of the same bulk cell used for the Fig. 1 dispersion.
epsCell = [p.eps1 p.eps2];
durations = [p.T/2 p.T/2];
response20 = temporal_finite_crystal_response(kGap,epsCell,[1 1], ...
    durations,20,p.epsBackground,p.mu,[1;0]);
response60 = temporal_finite_crystal_response(kGap,epsCell,[1 1], ...
    durations,60,p.epsBackground,p.mu,[1;0]);
coherent20 = abs(response20.forward)+abs(response20.backward);
coherent60 = abs(response60.forward)+abs(response60.backward);
fprintf(['Paper consistency audit at lambda=0.93 um: exact coherent TMM ' ...
    'peak N=20 is %.5g (ln=%.3f), whereas N=60 is %.5g ' ...
    '(ln=%.3f).\n'],coherent20,log(coherent20), ...
    coherent60,log(coherent60));
fprintf(['The FDTD below keeps the literal 60-period interval; panel (b) ' ...
    'uses the published natural-log color range 0..10 and clips larger gain.\n']);

% -------------------------------------------------------------------------
% Figure: natural logarithm of normalized |D|, so the published upper gap
% limit 10 corresponds to exp(10)=2.20e4.
% -------------------------------------------------------------------------
fig = figure('Color','w','Position',[60 60 840 460]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
plot_case(nexttile(tl),dataBand,caseBand,p,'(a)');
plot_case(nexttile(tl),dataGap,caseGap,p,'(b)');

outputFile = fullfile(outputDir,'fig2_fdtd_simulations.png');
exportgraphics(fig,outputFile,'Resolution',250);

audit = struct();
audit.coherentTmm20 = coherent20;
audit.coherentTmm60 = coherent60;
audit.logCoherentTmm20 = log(coherent20);
audit.logCoherentTmm60 = log(coherent60);
audit.log10CoherentTmm20 = log10(coherent20);
audit.log10CoherentTmm60 = log10(coherent60);
audit.paperClaimedGain = 2e4;
audit.fdtdMaximumGain60 = dataGap.maximumGain;
audit.fdtdToTmm60Ratio = dataGap.maximumGain/coherent60;
audit.note = ['The paper''s 60-period duration and approximately 20000 ' ...
    'gain cannot both hold for its stated binary PTC parameters.'];
if ~isfinite(audit.fdtdToTmm60Ratio) || ...
        audit.fdtdToTmm60Ratio < 0.7 || audit.fdtdToTmm60Ratio > 1.3
    error('The 60-period gap FDTD gain does not agree with exact TMM.');
end
save(fullfile(outputDir,'fig2_audit.mat'),'p','caseBand','caseGap','audit');
fprintf('Saved: %s\n',outputFile);
end

function data = run_or_load_case(caseDef,p,outputDir,forceRecompute)
dataFile = fullfile(outputDir,['fig2_' caseDef.name '_data.mat']);
if ~forceRecompute && exist(dataFile,'file')
    cached = load(dataFile,'data');
    if isfield(cached,'data') && isfield(cached.data,'modelVersion') && ...
            strcmp(cached.data.modelVersion,p.modelVersion)
        candidate = cached.data;
        required = {'D','logAmplitude','initialPeak','maximumGain', ...
            'zPeakAtStart','measuredPeaksAtEnd', ...
            'measuredVisiblePeaksAtFinal'};
        cacheIsFinite = all(isfield(candidate,required)) && ...
            all(isfinite(candidate.D(:))) && ...
            all(isfinite(candidate.logAmplitude(:))) && ...
            all(isfinite([candidate.initialPeak candidate.maximumGain ...
                candidate.zPeakAtStart candidate.measuredPeaksAtEnd(:).' ...
                candidate.measuredVisiblePeaksAtFinal(:).']));
        if cacheIsFinite
            data = candidate;
            fprintf('Loaded validated Fig. 2 cache: %s\n',dataFile);
            return;
        end
        warning('Ignoring non-finite or incomplete Fig. 2 cache: %s',dataFile);
    end
end

fprintf('Running Fig. 2 %s FDTD (lambda=%.2f um)...\n', ...
    caseDef.name,caseDef.lambda);

% Use a padded domain; only z=0..140 um is retained for the paper panel.
nIntervals = ceil(((p.zMax-p.zMin)/p.c0)/p.dx);
x = (0:nIntervals)*p.dx+p.zMin/p.c0;
z = p.c0*x;
xH = x+p.dx/2;
k = 2*pi*p.c0/caseDef.lambda;
xInitial = p.zInitial/p.c0;
pulse = @(xq) exp(-((xq-xInitial)/p.sigmaX).^2).* ...
    exp(1i*k*(xq-xInitial));

% The problem is spatially homogeneous, so k is exactly conserved.  Build
% a periodic analytic-signal grid and explicitly retain the physical pulse
% spectrum.  The cutoff removes less than 5e-6 of the Gaussian amplitude.
Nx = numel(x);
if mod(Nx,2) == 0
    fftOrders = [0:Nx/2-1 -Nx/2:-1];
else
    fftOrders = [0:(Nx-1)/2 -(Nx-1)/2:-1];
end
kGrid = (2*pi/(Nx*p.dx))*fftOrders;
filterFlatWidth = 5/p.sigmaX;
filterStopWidth = 7/p.sigmaX;
distanceFromCarrier = abs(kGrid-k);
spectralFilterMask = zeros(size(kGrid));
spectralFilterMask(distanceFromCarrier <= filterFlatWidth) = 1;
transition = distanceFromCarrier > filterFlatWidth & ...
    distanceFromCarrier < filterStopWidth;
s = (distanceFromCarrier(transition)-filterFlatWidth)/ ...
    (filterStopWidth-filterFlatWidth);
leftBump = exp(-1./(1-s));
rightBump = exp(-1./s);
spectralFilterMask(transition) = leftBump./(leftBump+rightBump);
E0 = ifft(fft(pulse(x)).*spectralFilterMask);
Hhalf0 = ifft(fft(pulse(xH+p.dt/2)).*spectralFilterMask);

nSteps = round(p.tStop/p.dt);
events = ptc_interface_times(p.tStart,p.tEnd,p.T);
cfg = struct();
cfg.x = x;
cfg.dt = p.dt;
cfg.nSteps = nSteps;
cfg.epsFun = @(xq,tq) finite_ptc_epsilon(xq,tq,p);
cfg.muFun = @(xq,tq) ones(size(xq));
cfg.E0 = E0;
cfg.Hhalf0 = Hhalf0;
cfg.boundary = 'periodic';
cfg.recordEvery = p.recordEvery;
cfg.storeFields = false;
cfg.storeD = true;
cfg.probeIndices = round(linspace(1,numel(x),9));
% A time crystal exponentially amplifies every grid Fourier component that
% lies in any momentum gap.  Single-precision roundoff therefore becomes a
% visible, spatially uniform false background during 60 periods.
cfg.precision = 'double';
cfg.progressBar = true;
cfg.temporalInterfaces = events;
cfg.spectralFilterMask = spectralFilterMask;
cfg.stabilityTimes = [0 p.tStart p.tStart+p.T/4 ...
    p.tStart+3*p.T/4 p.tEnd p.tStop];

out = fdtd1d_db(cfg);
crop = z >= p.zPlot(1) & z <= p.zPlot(2);
Dcrop = single(out.D(:,crop));      % solve in double; compact the cache only
initialPeak = max(abs(Dcrop(1,:)));
logAmplitude = log(max(abs(Dcrop)/initialPeak,exp(-35)));

% Simple trajectory checks that catch the old unit/domain failure.
[~,idStart] = min(abs(out.t-p.tStart));
[~,idEnd] = min(abs(out.t-p.tEnd));
zCrop = z(crop);
[~,peakStartId] = max(abs(Dcrop(idStart,:)));
zPeakAtStart = zCrop(peakStartId);
expectedStart = p.zInitial+p.c0*p.tStart;
if abs(zPeakAtStart-expectedStart) > 3
    error(['Pulse centre at PTC onset is %.3f um; expected approximately ' ...
        '%.3f um. The physical-unit mapping is inconsistent.'], ...
        zPeakAtStart,expectedStart);
end
if max(abs(Dcrop(idStart,:))) < 0.75*initialPeak
    error('The pulse was substantially attenuated before the PTC began.');
end

data = struct();
data.modelVersion = p.modelVersion;
data.name = caseDef.name;
data.lambda = caseDef.lambda;
data.k = k;
data.z = zCrop;
data.t = out.t;
data.D = Dcrop;
data.logAmplitude = logAmplitude;
data.initialPeak = initialPeak;
data.zPeakAtStart = zPeakAtStart;
data.expectedZPeakAtStart = expectedStart;
data.peakAtStart = max(abs(Dcrop(idStart,:)));
data.peakAtEnd = max(abs(Dcrop(idEnd,:)));
data.maximumGain = max(abs(Dcrop(:)))/initialPeak;
data.sampledCourant = out.sampledCourant;
data.dxNumerical = out.dx;
data.dxMicrometre = p.c0*out.dx;
data.dt = out.dt;
data.temporalInterfaces = out.temporalInterfaces;
data.filterFlatWidth = filterFlatWidth;
data.filterStopWidth = filterStopWidth;
data.spectralFilterMask = spectralFilterMask;
data.parameters = p;

if strcmp(caseDef.name,'inband')
    groupSpeedX = floquet_group_speed(k,p);
    expectedEnd = expectedStart+[-1 1]*p.c0*groupSpeedX* ...
        (p.tEnd-p.tStart);
    actualEnd = window_peaks(zCrop,abs(Dcrop(idEnd,:)),expectedEnd,12);
    expectedFinal = reshape(expectedEnd(:)+ ...
        [-1 1]*p.c0*(p.tStop-p.tEnd),1,[]);
else
    groupSpeedX = 0;
    expectedEnd = [expectedStart expectedStart];
    actualEnd = window_peaks(zCrop,abs(Dcrop(idEnd,:)),expectedStart,12);
    expectedFinal = expectedStart+[-1 1]*p.c0*(p.tStop-p.tEnd);
end
[~,idFinal] = min(abs(out.t-p.tStop));
insidePlot = expectedFinal >= p.zPlot(1) & expectedFinal <= p.zPlot(2);
actualFinal = window_peaks(zCrop,abs(Dcrop(idFinal,:)), ...
    expectedFinal(insidePlot),12);
data.floquetGroupSpeedOverC = groupSpeedX;
data.expectedPeaksAtEnd = expectedEnd;
data.measuredPeaksAtEnd = actualEnd;
data.expectedVisiblePeaksAtFinal = expectedFinal(insidePlot);
data.measuredVisiblePeaksAtFinal = actualFinal;
if ~all(isfinite(Dcrop(:))) || ~all(isfinite(logAmplitude(:))) || ...
        ~all(isfinite([initialPeak data.maximumGain zPeakAtStart ...
            actualEnd(:).' actualFinal(:).']))
    error('The Fig. 2 FDTD produced non-finite fields or trajectory metrics.');
end
fprintf('  Ridge audit at tEnd expected %s um, measured %s um.\n', ...
    mat2str(expectedEnd,4),mat2str(actualEnd,4));
fprintf('  Ridge audit at tStop expected %s um, measured %s um.\n', ...
    mat2str(expectedFinal(insidePlot),4),mat2str(actualFinal,4));
if any(abs(actualEnd(:)-expectedEnd(1:numel(actualEnd)).') > 8) || ...
        any(abs(actualFinal(:)-expectedFinal(insidePlot).') > 8)
    error('The simulated pulse ridges deviate from the Floquet trajectory audit.');
end

if strcmp(caseDef.name,'inband') && data.maximumGain > 20
    error(['The in-band run developed nonphysical broadband gain ' ...
        '(maximum %.4g). Increase numerical precision or filtering.'], ...
        data.maximumGain);
end

save(dataFile,'data','-v7.3');
fprintf(['  Nx=%d, Nt=%d, dz=%.5f um, onset centre=%.3f um, ' ...
    'maximum |D|/|D0|=%.5g.\n'],numel(x),nSteps, ...
    data.dxMicrometre,zPeakAtStart,data.maximumGain);
fprintf('  Saved: %s\n',dataFile);
end

function plot_case(ax,data,caseDef,p,panelLabel)
imagesc(ax,data.z,data.t,data.logAmplitude);
axis(ax,'xy');
colormap(ax,turbo(256));
clim(ax,caseDef.clim);
hold(ax,'on');
zMarker = p.zPlot(2)-5;
plot(ax,[zMarker zMarker],[p.tStart p.tEnd],'w-','LineWidth',1.2);
plot(ax,zMarker,p.tStart,'w_','MarkerSize',8,'LineWidth',1.2);
plot(ax,zMarker,p.tEnd,'w_','MarkerSize',8,'LineWidth',1.2);
text(ax,zMarker-2,0.5*(p.tStart+p.tEnd),'\epsilon(t)', ...
    'Color','w','HorizontalAlignment','right','FontWeight','bold');
xlim(ax,p.zPlot);
ylim(ax,[0 p.tStop]);
xlabel(ax,'z [\mum]');
ylabel(ax,'t [fs]');
title(ax,sprintf('%s   %s',panelLabel,caseDef.label), ...
    'FontWeight','normal');
cb = colorbar(ax);
cb.Label.String = 'ln(|D|/max|D(t=0)|)';
set(ax,'FontSize',9,'Layer','top');
box(ax,'on');
end

function interfaces = ptc_interface_times(tStart,tEnd,T)
interfaces = tStart:T/2:tEnd;
end

function epsValue = finite_ptc_epsilon(x,t,p)
% Half-open PTC window [tStart,tEnd), beginning with a full eps1 half-cell.
if t < p.tStart || t >= p.tEnd
    scalarEps = p.epsBackground;
else
    tau = mod(t-p.tStart,p.T);
    if tau < p.T/2
        scalarEps = p.eps1;
    else
        scalarEps = p.eps2;
    end
end
epsValue = scalarEps*ones(size(x));
end

function speed = floquet_group_speed(k,p)
dk = 1e-5*max(1,abs(k));
probeK = [k-dk k+dk];
b = temporal_crystal_bands(probeK, ...
    [p.eps1 p.eps2],[1 1],[p.T/2 p.T/2]);
h = min(1,max(-1,real(b.halfTrace)));
omega = acos(h)/p.T;
speed = abs(diff(omega)/(2*dk));
end

function positions = window_peaks(z,amplitude,expected,halfWidth)
positions = nan(size(expected));
for j = 1:numel(expected)
    ids = find(abs(z-expected(j)) <= halfWidth);
    if isempty(ids), continue; end
    [~,localId] = max(amplitude(ids));
    positions(j) = z(ids(localId));
end
end
