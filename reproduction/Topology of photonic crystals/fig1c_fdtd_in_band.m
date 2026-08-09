function fig1c_fdtd_in_band()
%FIG1C_FDTD_IN_BAND Compatibility entry point for published Fig. 2(a).
%
% The old standalone script duplicated an incorrect normalized-unit FDTD
% setup.  Keep one validated implementation so the legacy arXiv figure name
% cannot silently regenerate the obsolete result.

warning('fig1c_fdtd_in_band:renamed', ...
    ['The published panel is Fig. 2(a). Running the validated combined ' ...
    'Fig. 2 reproduction (cached data are reused when compatible).']);
fig2_fdtd_simulations(false);
end
