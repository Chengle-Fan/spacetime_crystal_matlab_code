%RUN_TL_TMM Temporal monodromy bands of the paper-seeded SSPP circuit.
% This TMM advances time at fixed Bloch k; it is not a spatial ABCD cascade.

clear; clc; close all;

physicalCfg = struct();
physicalCfg.topology = 'sspp';
physicalCfg.allowAssumptions = true;
physicalCfg.fmHz = 675e6;
model = tl_build_model(physicalCfg);

kScan = linspace(-pi/model.cell.a,pi/model.cell.a,181);
tmmCfg = struct('temporalSlices',1024);
tmmBands = tl_tmm_bands(model,kScan,tmmCfg);

figure('Color','w','Position',[100 100 900 620]);
tiledlayout(2,1,'TileSpacing','compact');
nexttile;
plot(tmmBands.kaOverPi,real(tmmBands.omegaFolded)/(2*pi*1e6), ...
    'LineWidth',1.3);
xlabel('ka/\pi');
ylabel('Re(f) (MHz)');
title('Transmission-line temporal TMM: first Floquet zone');
grid on; box on;
nexttile;
plot(tmmBands.kaOverPi,imag(tmmBands.omegaFolded)/(2*pi*1e6), ...
    'LineWidth',1.3);
xlabel('ka/\pi');
ylabel('Im(f) (MHz)');
grid on; box on;

fprintf('TMM complete: topology=%s, integration=%s, k points=%d.\n', ...
    model.kind,tmmBands.integrationName,numel(kScan));

