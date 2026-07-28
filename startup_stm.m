function rootDir = startup_stm()
%STARTUP_STM Add all folders in this package to the MATLAB path.
%
%   rootDir = startup_stm();

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
addpath(fullfile(rootDir, 'core'));
addpath(fullfile(rootDir, 'tmm'));
addpath(fullfile(rootDir, 'fdtd'));
addpath(fullfile(rootDir, 'topology'));
addpath(fullfile(rootDir, 'demos'));
addpath(fullfile(rootDir, 'tests'));

outputDir = fullfile(rootDir, 'output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

fprintf('Space-time media MATLAB research package v2 ready.\n');
fprintf('Root: %s\n', rootDir);
end
