%RUN_TL_FDTD_FIELD Finite-chain voltage/current field simulation.
% I is a circuit proxy for local magnetic field and V for electric field;
% neither is a replacement for a full-wave PCB field solution.

clear; clc; close all;

physicalCfg = struct();
physicalCfg.topology = 'sspp';
physicalCfg.allowAssumptions = true;
physicalCfg.cellCount = 64;
model = tl_build_model(physicalCfg);

period = model.modulation.period;
stepsPerPeriod = max(256,ceil(period/(0.5*model.derived.maximumLeapfrogDt)));
dt = period/stepsPerPeriod;
simulationPeriods = 8;
pulseCenter = 0.8*period;
pulseWidth = 0.18*period;
carrierHz = model.modulation.fmHz/2;
sourceVoltage = @(time) exp(-0.5*((time-pulseCenter)/pulseWidth).^2).* ...
    cos(2*pi*carrierHz*(time-pulseCenter));

fdtdCfg = struct();
fdtdCfg.dt = dt;
fdtdCfg.nSteps = simulationPeriods*stepsPerPeriod;
fdtdCfg.recordEvery = 2;
fdtdCfg.boundaryType = 'matched';
fdtdCfg.modulationEnabled = true;
fdtdCfg.modulationStart = 2*period;
fdtdCfg.modulationEnd = Inf;
fdtdCfg.precision = 'single';
fdtdCfg.source = struct('type','thevenin','node',1, ...
    'impedance',model.ports.Zsource,'waveformFcn',sourceVoltage);
field = tl_fdtd1d(model,fdtdCfg);

figure('Color','w','Position',[100 100 1000 650]);
tiledlayout(2,1,'TileSpacing','compact');
nexttile;
imagesc(field.node.x/model.cell.a,field.node.t/period,real(field.node.V));
axis xy; colorbar;
xlabel('node position x/a');
ylabel('t/T');
title('Node voltage V (real part)');
nexttile;
if strcmp(model.kind,'crow')
    magneticChannel = field.resonator;
    currentLabel = 'Centered resonator current I_0 (real part)';
else
    magneticChannel = field.branch;
    currentLabel = 'Centered series current I (real part)';
end
imagesc(magneticChannel.x/model.cell.a,field.node.t/period, ...
    real(magneticChannel.IAtNodeTime));
axis xy; colorbar;
xlabel('branch position x/a');
ylabel('t/T');
title(currentLabel);

fprintf(['FDTD complete: dt*omegaMax=%.4g, records=%d, final relative ' ...
    'energy-ledger residual=%.3g.\n'],field.grid.stabilityNumber, ...
    numel(field.node.t),field.energy.relativeLedgerResidual(end));
