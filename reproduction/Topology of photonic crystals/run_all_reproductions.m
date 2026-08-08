function run_all_reproductions()
%RUN_ALL_REPRODUCTIONS  Run all reproduction scripts for Lustig et al. (2018).
%
%   This script sequentially executes all figure reproduction codes
%   for the paper "Topology of photonic time-crystals."
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
fprintf('  "Topology of photonic time-crystals"\n');
fprintf('  arXiv:1803.08731v1\n');
fprintf('==============================================================\n\n');

% -------------------------------------------------------------------------
% Figure 1(a): PTC schematic (fast)
% -------------------------------------------------------------------------
fprintf('\n%%%% Fig. 1(a): PTC schematic %%%%\n');
try
    run(fullfile(scriptDir, 'fig1a_ptc_schematic.m'));
    fprintf('✓ Fig. 1(a) completed.\n');
catch ME
    fprintf('✗ Fig. 1(a) failed: %s\n', ME.message);
end

% -------------------------------------------------------------------------
% Figure 1(b): Band structure with Zak phases (moderate)
% -------------------------------------------------------------------------
fprintf('\n%%%% Fig. 1(b): Floquet band structure %%%%\n');
try
    run(fullfile(scriptDir, 'fig1b_band_structure.m'));
    fprintf('✓ Fig. 1(b) completed.\n');
catch ME
    fprintf('✗ Fig. 1(b) failed: %s\n', ME.message);
end

% -------------------------------------------------------------------------
% Figure 1(c): FDTD — pulse in band (SLOW)
% -------------------------------------------------------------------------
fprintf('\n%%%% Fig. 1(c): FDTD — pulse in band (SLOW, ~3-10 min) %%%%\n');
fprintf('Proceeding with FDTD simulation...\n');
try
    run(fullfile(scriptDir, 'fig1c_fdtd_in_band.m'));
    fprintf('✓ Fig. 1(c) completed.\n');
catch ME
    fprintf('✗ Fig. 1(c) failed: %s\n', ME.message);
end

% -------------------------------------------------------------------------
% Figure 1(d): FDTD — pulse in gap (SLOW)
% -------------------------------------------------------------------------
fprintf('\n%%%% Fig. 1(d): FDTD — pulse in gap (SLOW, ~3-10 min) %%%%\n');
fprintf('Proceeding with FDTD simulation...\n');
try
    run(fullfile(scriptDir, 'fig1d_fdtd_in_gap.m'));
    fprintf('✓ Fig. 1(d) completed.\n');
catch ME
    fprintf('✗ Fig. 1(d) failed: %s\n', ME.message);
end

% -------------------------------------------------------------------------
% Figure 3: Relative phase for first 6 gaps (moderate)
% -------------------------------------------------------------------------
fprintf('\n%%%% Fig. 3: Relative phase for first 6 gaps %%%%\n');
try
    run(fullfile(scriptDir, 'fig3_relative_phase.m'));
    fprintf('✓ Fig. 3 completed.\n');
catch ME
    fprintf('✗ Fig. 3 failed: %s\n', ME.message);
end

% -------------------------------------------------------------------------
% Zak Phase: dedicated computation (fast)
% -------------------------------------------------------------------------
fprintf('\n%%%% Zak Phase: dedicated Wilson-loop computation %%%%\n');
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
