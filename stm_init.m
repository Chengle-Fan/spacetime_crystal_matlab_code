function rootDir = stm_init()
%STM_INIT Add all folders of the Space-Time Media toolbox to the MATLAB path.
%
%   rootDir = stm_init();
%
% This is the entry point for the STM (Space-Time Media) toolbox. It adds
% the core engine, TMM, FDTD, topology, examples, and tests directories to
% the MATLAB search path and creates the output/ folder if it does not exist.

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
addpath(fullfile(rootDir, 'core'));
addpath(fullfile(rootDir, 'tmm'));
addpath(fullfile(rootDir, 'fdtd'));
addpath(fullfile(rootDir, 'topology'));
addpath(fullfile(rootDir, 'examples'));
addpath(fullfile(rootDir, 'tests'));

outputDir = fullfile(rootDir, 'output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

fprintf('Space-Time Media (STM) toolbox ready.\n');
fprintf('Root: %s\n', rootDir);
end
