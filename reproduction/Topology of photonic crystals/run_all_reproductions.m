function run_all_reproductions()
%RUN_ALL_REPRODUCTIONS  Run all reproduction scripts for Lustig et al. (2018).
%
%   This script sequentially executes all figure reproduction codes
%   for the published paper "Topological aspects of photonic time crystals,"
%   Optica 5, 1390-1395 (2018).  DOI: 10.1364/OPTICA.5.001390
%
%   Figure mapping (arXiv v1 → Optica published):
%     Old Fig 1(a,b)  →  New Fig 1(a,b): PTC schematic + band structure
%     Old Fig 1(c,d)  →  New Fig 2(a,b): FDTD in-band + in-gap
%     New             →  New Fig 3:      Conceptual schematic (not computational)
%     Old Fig 3(a-f)  →  New Fig 4(a-f): Relative phase in first 6 gaps
%     New             →  New Fig 5(a-c): Temporal topological edge states
%
%   Output figures and data are saved in the output/ subdirectory.
%
%   Usage:
%     run('run_all_reproductions')

% --- Add parent toolbox to path ---
rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(rootDir, 'startup_stm.m'));

% Get the directory of this script
scriptDir = fileparts(mfilename('fullpath'));

fprintf('==============================================================\n');
fprintf('  Reproduction of Lustig et al. (2018)\n');
fprintf('  "Topological aspects of photonic time crystals"\n');
fprintf('  Optica 5, 1390-1395 (2018)\n');
fprintf('  DOI: 10.1364/OPTICA.5.001390\n');
fprintf('==============================================================\n\n');

% -------------------------------------------------------------------------
% Figure 1: PTC schematic + band structure with Zak phases (fast)
% -------------------------------------------------------------------------
fprintf('\n%%%% Fig. 1: PTC schematic + Floquet band structure %%%%\n');
try
    run(fullfile(scriptDir, 'fig1_ptc_bands.m'));
    fprintf('✓ Fig. 1 completed.\n');
catch ME
    fprintf('✗ Fig. 1 failed: %s\n', ME.message);
end

% -------------------------------------------------------------------------
% Figure 2: FDTD simulations — in band + in gap (SLOW)
% -------------------------------------------------------------------------
fprintf('\n%%%% Fig. 2: FDTD simulations (SLOW, ~6-20 min total) %%%%\n');
fprintf('Set doFDTD=false in fig2_fdtd_simulations.m to skip FDTD and\n');
fprintf('load previously saved data.\n');
try
    run(fullfile(scriptDir, 'fig2_fdtd_simulations.m'));
    fprintf('✓ Fig. 2 completed.\n');
catch ME
    fprintf('✗ Fig. 2 failed: %s\n', ME.message);
end

% -------------------------------------------------------------------------
% Figure 4: Relative phase for first 6 gaps (moderate)
% -------------------------------------------------------------------------
fprintf('\n%%%% Fig. 4: Relative phase for first 6 momentum gaps %%%%\n');
try
    run(fullfile(scriptDir, 'fig4_relative_phase.m'));
    fprintf('✓ Fig. 4 completed.\n');
catch ME
    fprintf('✗ Fig. 4 failed: %s\n', ME.message);
end

% -------------------------------------------------------------------------
% Figure 5: Temporal topological edge states (moderate)
% -------------------------------------------------------------------------
fprintf('\n%%%% Fig. 5: Temporal topological edge states %%%%\n');
try
    run(fullfile(scriptDir, 'fig5_temporal_edge_state.m'));
    fprintf('✓ Fig. 5 completed.\n');
catch ME
    fprintf('✗ Fig. 5 failed: %s\n', ME.message);
end

% -------------------------------------------------------------------------
% Bonus: Zak phase detailed computation (fast)
% -------------------------------------------------------------------------
fprintf('\n%%%% Bonus: Detailed Zak phase computation (Wilson loop) %%%%\n');
try
    run(fullfile(scriptDir, 'compute_zak_phases.m'));
    fprintf('✓ Zak phase computation completed.\n');
catch ME
    fprintf('✗ Zak phase failed: %s\n', ME.message);
end

% -------------------------------------------------------------------------
% Summary
% -------------------------------------------------------------------------
fprintf('\n==============================================================\n');
fprintf('  All reproductions completed.\n');
fprintf('  Output files are in: %s\n', fullfile(scriptDir, 'output'));
fprintf('==============================================================\n');
end
