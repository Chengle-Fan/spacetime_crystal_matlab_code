%RUN_TL_PWE Floquet PWE bands of the paper-seeded SSPP circuit.
% Run this script from the Transmission line directory.  Values not given
% by the reference PDF are deliberately enabled as simulation assumptions.

clear; clc; close all;

physicalCfg = struct();
physicalCfg.topology = 'sspp';
physicalCfg.allowAssumptions = true;
physicalCfg.fmHz = 675e6;
model = tl_build_model(physicalCfg);

kScan = linspace(-pi/model.cell.a,pi/model.cell.a,181);
pweCfg = struct('Mtime',5,'Nt',256);
fourier = tl_pwe_fourier(model,pweCfg);
pweBands = tl_pwe_bands(fourier,model,kScan,pweCfg);

figure('Color','w','Position',[100 100 900 620]);
tiledlayout(2,1,'TileSpacing','compact');
nexttile;
plot(pweBands.kaOverPi,real(pweBands.omegaSelected)/(2*pi*1e6), ...
    'LineWidth',1.3);
xlabel('ka/\pi');
ylabel('Re(f) (MHz)');
title('Transmission-line Floquet PWE: selected physical branches');
grid on; box on;
nexttile;
plot(pweBands.kaOverPi,imag(pweBands.omegaSelected)/(2*pi*1e6), ...
    'LineWidth',1.3);
xlabel('ka/\pi');
ylabel('Im(f) (MHz)');
grid on; box on;

fprintf(['PWE complete: topology=%s, M=%d, Nt=%d.  Positive Im(omega) ' ...
    'means growth under exp(i*k*x-i*omega*t).\n'], ...
    model.kind,pweCfg.Mtime,pweCfg.Nt);

