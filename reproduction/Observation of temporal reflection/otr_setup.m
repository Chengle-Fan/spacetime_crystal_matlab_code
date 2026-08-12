function paths = otr_setup()
%OTR_SETUP Add this reproduction and the parent toolbox to the MATLAB path.
%
%   paths = OTR_SETUP() is safe to call from any working directory.  All
%   generated figures and MAT files are kept below this reproduction
%   directory, as requested; no source file outside it is modified.

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fileparts(scriptDir));
parentReady = exist('fdtd1d_db','file') == 2 && ...
    exist('temporal_crystal_monodromy','file') == 2;
if ~parentReady && exist(fullfile(rootDir,'startup_stm.m'),'file')
    addpath(rootDir);
    startup_stm();
end
% Add this folder last so the reproduction's paper-specific helpers retain
% priority if the parent package later acquires a same-named function.
addpath(scriptDir,'-begin');

paths.root = rootDir;
paths.reproduction = scriptDir;
paths.data = fullfile(scriptDir,'data','Time_Interface_Source_Data');
paths.output = fullfile(scriptDir,'output');
if ~exist(paths.output,'dir')
    mkdir(paths.output);
end
end
