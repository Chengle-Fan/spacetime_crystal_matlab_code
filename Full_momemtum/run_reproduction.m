%RUN_REPRODUCTION Reproduce the computable panels of Figs. 2 and 3.
% The numerical kernels are loaded from ../Transmission line; this folder
% contains only paper parameters, scans, diagnostics, and plotting code.

clear; clc; close all;

options = struct();
options.quick = false;
options.runFields = true;
options.runPhaseMap = true;
options.runFftBands = true;
options.saveOutputs = true;

results = reproduce_figures_2_3(options);
