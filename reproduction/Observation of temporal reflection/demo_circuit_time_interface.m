function results = demo_circuit_time_interface()
%DEMO_CIRCUIT_TIME_INTERFACE Kirchhoff/MNA reproduction of the experiment.
%
% This demo solves a Kirchhoff/MNA transmission-line ladder.  It does not
% call the repository's continuous-medium FDTD solver.  Thirty series RL
% cells and physical shunt capacitors approximate the paper's meandered
% line; all shunt loads are inserted at one grid-aligned instant.

paths = otr_setup();

%% Paper-scale circuit and asymmetric broadband input
cfg = struct();
cfg.physicsMode = 'ideal-capacitance';
cfg.nCells = 30;
cfg.d = 0.2080;
cfg.Z0 = 50;
cfg.epsEff = 8.36;
cfg.Cload = 82e-12;
cfg.Rseries = 0.15;
cfg.Gshunt = 0;
cfg.includeParasitics = true;
cfg.Cnearest = 2e-12;
cfg.CnextNearest = 1e-12;
model = otr_build_mna(cfg);

dt = 0.05e-9;
tEnd = 310e-9;
t = 0:dt:tEnd;
tSwitch = 105e-9;
assert(abs(t(round(tSwitch/dt)+1)-tSwitch) < 10*eps(tSwitch), ...
    'The switching event must lie exactly on the grid.');

% Two separated derivative-of-Gaussian wavelets create an unmistakable
% order code: a small pulse is launched before a large pulse.  Their zero
% DC content avoids slowly charging the shunt capacitors.
pulseWidth = 5.0e-9;
smallCentre = 20e-9;
largeCentre = 40e-9;
dog = @(time,centre,width) ((time-centre)/width).* ...
    exp(-0.5*((time-centre)/width).^2);
sourceWave = 0.42*dog(t,smallCentre,pulseWidth) + ...
    1.00*dog(t,largeCentre,pulseWidth);

scheduleSpec = struct('nSwitch',model.nSwitches, ...
    'eventTimes',tSwitch,'states',[0 1],'riseTime',0);
schedule = otr_switch_schedule(t,scheduleSpec);
sim = otr_simulate_circuit(model,t,sourceWave,schedule);
energy = otr_circuit_energy(model,sim);

%% Independent boundary-condition unit tests
% A one-cell ladder with prescribed pre-event voltage makes S14/S15
% directly testable without relying on pulse extraction.
testCfg = cfg;
testCfg.nCells = 1;
testCfg.includeParasitics = false;
testCfg.Rseries = 0;
testModel = otr_build_mna(testCfg);
testTime = [0 1e-12 2e-12];
testInitial = zeros(testModel.nState,1);
testInitial(testModel.index.voltage) = [0.2; 1.0];

onSchedule = otr_switch_schedule(testTime,struct( ...
    'nSwitch',testModel.nSwitches,'eventTimes',testTime(2), ...
    'states',[0 1]));
onResult = otr_simulate_circuit(testModel,testTime,0,onSchedule, ...
    struct('initialState',testInitial));
onEvent = onResult.eventRecords(1);
assert(onEvent.chargeResidual < 1e-11, ...
    'ON boundary failed conductor-charge conservation (S14).');
assert(onEvent.statePlus(testModel.index.lineCurrent) == ...
    onEvent.stateMinus(testModel.index.lineCurrent), ...
    'Line-inductor current changed across an instantaneous interface.');

offInitial = onEvent.statePlus;
offSchedule = otr_switch_schedule(testTime,struct( ...
    'nSwitch',testModel.nSwitches,'eventTimes',testTime(2), ...
    'states',[1 0]));
offResult = otr_simulate_circuit(testModel,testTime,0,offSchedule, ...
    struct('initialState',offInitial));
offEvent = offResult.eventRecords(1);
offVoltageJump = offEvent.statePlus(testModel.index.voltage)- ...
    offEvent.stateMinus(testModel.index.voltage);
assert(norm(offVoltageJump,inf) < 1e-13, ...
    'OFF boundary failed voltage continuity (S15).');

% The energy balance uses MNA midpoint power and explicitly records event
% energy.  Roundoff-sized residual is expected for the trapezoidal rule.
assert(energy.relativeBalanceResidual < 5e-8, ...
    'Circuit energy balance residual is unexpectedly large.');

%% Time-reversal ordering and polarity checks at a near-input probe
probeNode = 4;
probe = sim.voltage(probeNode,:);
oneWayToProbe = (probeNode-1)*model.parameters.cellDelay;
incidentCentres = [smallCentre largeCentre]+oneWayToProbe;

% Under k conservation, the loaded phase velocity is approximately
% v2=v1*sqrt(Cbase/(Cbase+Cload)).  A feature a time delta before the
% interface returns after delta*v1/v2, hence the reversed ordering.
capRatio = model.parameters.Ccell/(model.parameters.Ccell+ ...
    mean(model.parameters.Cload));
frequencyRatioLumped = sqrt(capRatio);
reflectionCentres = tSwitch + ...
    (tSwitch-incidentCentres)/frequencyRatioLumped;
reflectionAmplitude = local_peak_near(t,probe,reflectionCentres,10e-9);

assert(reflectionCentres(2) < reflectionCentres(1), ...
    'Predicted temporal-reflection feature ordering was not reversed.');
assert(abs(reflectionAmplitude(2)) > abs(reflectionAmplitude(1)), ...
    'The large-late input feature did not return before the small feature.');

% The signed correlation of the returned large feature with a time-reversed
% copy of the incident large feature must be negative for C_on>C_off.
correlationWindow = 10e-9;
incidentMask = abs(t-incidentCentres(2)) <= correlationWindow;
reflectedMask = abs(t-reflectionCentres(2)) <= ...
    correlationWindow/frequencyRatioLumped;
incidentSegment = probe(incidentMask);
reflectedSegment = probe(reflectedMask);
incidentAxis = linspace(-1,1,numel(incidentSegment));
reflectedAxis = linspace(-1,1,numel(reflectedSegment));
referenceReversed = interp1(-incidentAxis,incidentSegment, ...
    reflectedAxis,'linear',0);
signedTRCorrelation = sum(referenceReversed.*reflectedSegment)/sqrt( ...
    sum(referenceReversed.^2)*sum(reflectedSegment.^2));
assert(signedTRCorrelation < 0, ...
    'The temporal reflection did not have the expected negative polarity.');

%% Spatial-time voltage map and port/probe traces
outputDirectory = paths.output;

figureTime = figure('Color','w','Name','Kirchhoff/MNA time interface');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile;
imagesc(t*1e9,0:model.nCells,sim.voltage);
axis xy;
hold on;
xline(tSwitch*1e9,'w--','LineWidth',1.4);
hold off;
xlabel('Time (ns)');
ylabel('Circuit node');
title('Discrete MNA node voltages: forward pulse and temporal reflection');
colormap(gca,parula(256));
colorbar;

nexttile;
plot(t*1e9,sourceWave,'Color',[0.10 0.10 0.10],'LineWidth',1.1, ...
    'DisplayName','Thevenin source');
hold on;
plot(t*1e9,probe,'Color',[0.18 0.48 0.78],'LineWidth',1.2, ...
    'DisplayName',sprintf('Node %d voltage',probeNode));
xline(tSwitch*1e9,'Color',[0.20 0.65 0.25],'LineStyle','--', ...
    'DisplayName','Switch ON');
hold off;
grid on;
xlabel('Time (ns)');
ylabel('Voltage (V)');
title(sprintf(['Small-then-large input returns large-then-small; ' ...
    'signed TR correlation = %.3f'],signedTRCorrelation));
legend('Location','best');

nexttile;
yyaxis left;
plot(t*1e9,energy.stored*1e9,'LineWidth',1.1);
ylabel('Stored energy (nJ)');
yyaxis right;
stairs(t*1e9,mean(sim.weights,1),'Color',[0.18 0.62 0.25], ...
    'LineWidth',1.2);
ylabel('Mean switch weight');
xlabel('Time (ns)');
grid on;
title(sprintf('MNA energy balance relative residual %.2g', ...
    energy.relativeBalanceResidual));
drawnow;
print(figureTime,fullfile(outputDirectory,'circuit_time_interface.png'), ...
    '-dpng','-r180');
savefig(figureTime,fullfile(outputDirectory,'circuit_time_interface.fig'));

%% Fixed-k frequency translation from the actual circuit dispersion
% For a periodic series-L/shunt-C ladder, equal Bloch k before and after
% switching yields sin(omega_1*sqrt(LC)/2) / sin(...)=sqrt(C ratios).
% We evaluate the exact lumped ladder relation and compare its low-k limit.
k = linspace(0,0.85*pi/model.parameters.d,300);
Lcell = mean(model.parameters.Lseries);
Coff = model.parameters.Ccell;
Con = Coff+mean(model.parameters.Cload);
omegaOff = 2/sqrt(Lcell*Coff).*sin(k*model.parameters.d/2);
omegaOn = 2/sqrt(Lcell*Con).*sin(k*model.parameters.d/2);

figureDispersion = figure('Color','w','Name','Fixed-k frequency translation');
plot(k,omegaOff/(2*pi*1e6),'Color',[0.15 0.45 0.82], ...
    'LineWidth',1.7,'DisplayName','Switch OFF');
hold on;
plot(k,omegaOn/(2*pi*1e6),'Color',[0.85 0.24 0.25], ...
    'LineWidth',1.7,'DisplayName','Switch ON');
selected = round(0.55*numel(k));
plot([k(selected) k(selected)], ...
    [omegaOn(selected) omegaOff(selected)]/(2*pi*1e6), ...
    'k--','LineWidth',1.2,'DisplayName','Temporal transition (fixed k)');
hold off;
grid on;
xlabel('Bloch wavenumber k (rad/m)');
ylabel('Frequency (MHz)');
title(sprintf('Lumped MNA fixed-k redshift, f_{on}/f_{off}=%.3f', ...
    omegaOn(selected)/omegaOff(selected)));
legend('Location','northwest');
drawnow;
print(figureDispersion,fullfile(outputDirectory, ...
    'circuit_fixed_k_frequency_translation.png'),'-dpng','-r180');
savefig(figureDispersion,fullfile(outputDirectory, ...
    'circuit_fixed_k_frequency_translation.fig'));

fprintf('Discrete MNA time-interface demo completed.\n');
fprintf('  States: %d (%d node voltages, %d line currents)\n', ...
    model.nState,model.nNodes,model.nLineBranches);
fprintf('  Switch event: %.3f ns, grid error %.3g s\n', ...
    tSwitch*1e9,abs(schedule.eventTimes-tSwitch));
fprintf('  Lumped low-k frequency ratio: %.4f\n',frequencyRatioLumped);
fprintf('  Signed time-reversal correlation: %.4f\n',signedTRCorrelation);
fprintf('  Energy-balance relative residual: %.3g\n', ...
    energy.relativeBalanceResidual);
fprintf('  ON charge residual: %.3g; OFF voltage jump: %.3g V\n', ...
    onEvent.chargeResidual,norm(offVoltageJump,inf));
fprintf('  Figures saved in: %s\n',outputDirectory);

results = struct();
results.parameters = model.parameters;
results.time = t;
results.source = sourceWave;
results.probeNode = probeNode;
results.probe = probe;
results.switchTime = tSwitch;
results.frequencyRatioLumped = frequencyRatioLumped;
results.signedTRCorrelation = signedTRCorrelation;
results.energyBalanceResidual = energy.relativeBalanceResidual;
results.onChargeResidual = onEvent.chargeResidual;
results.offVoltageJump = norm(offVoltageJump,inf);
results.k = k;
results.frequencyOff = omegaOff/(2*pi);
results.frequencyOn = omegaOn/(2*pi);
save(fullfile(outputDirectory,'demo_circuit_time_interface.mat'), ...
    'results','-v7.3');
end

function peaks = local_peak_near(time,signal,centres,halfWidth)
peaks = zeros(size(centres));
for centreIndex = 1:numel(centres)
    mask = abs(time-centres(centreIndex)) <= halfWidth;
    segment = signal(mask);
    [~,maximumIndex] = max(abs(segment));
    peaks(centreIndex) = segment(maximumIndex);
end
end
