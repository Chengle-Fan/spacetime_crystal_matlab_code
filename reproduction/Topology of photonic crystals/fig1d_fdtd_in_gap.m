function fig1d_fdtd_in_gap()
%FIG1D_FDTD_IN_GAP Compatibility entry point for published Fig. 2(b).
%
% The old standalone script duplicated an incorrect normalized-unit FDTD
% setup.  Keep one validated implementation so the legacy arXiv figure name
% cannot silently regenerate the obsolete result.

warning('fig1d_fdtd_in_gap:renamed', ...
    ['The published panel is Fig. 2(b). Running the validated combined ' ...
    'Fig. 2 reproduction (cached data are reused when compatible).']);
fig2_fdtd_simulations(false);
end
