function run_all_reproductions(forceFDTD)
%RUN_ALL_REPRODUCTIONS Generate the corrected published Figs. 1, 2, 4, 5.
%
%   run_all_reproductions()       reuses validated Fig. 2 cache files
%   run_all_reproductions(true)   recomputes both Fig. 2 FDTD cases

if nargin < 1
    forceFDTD = false;
end
validateattributes(forceFDTD,{'logical','numeric'},{'scalar'});
forceFDTD = logical(forceFDTD);

rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(rootDir,'startup_stm.m'));
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

fprintf('==============================================================\n');
fprintf(' Lustig, Sharabi & Segev, Optica 5, 1390 (2018)\n');
fprintf(' Corrected reproduction of published Figs. 1, 2, 4 and 5\n');
fprintf('==============================================================\n');

jobs = { ...
    'Fig. 1: bands and Zak labels', @() fig1_ptc_bands(); ...
    'Fig. 2: D/B-Yee FDTD',        @() fig2_fdtd_simulations(forceFDTD); ...
    'Fig. 4: first six gap phases',@() fig4_relative_phase(); ...
    'Fig. 5: temporal edge pulse', @() fig5_temporal_edge_state()};

for j = 1:size(jobs,1)
    fprintf('\n--- %s ---\n',jobs{j,1});
    tic;
    jobs{j,2}();
    fprintf('Completed in %.2f s.\n',toc);
end

fprintf('\nAll corrected reproductions completed.\n');
fprintf('Output directory: %s\n',fullfile(scriptDir,'output'));
end
