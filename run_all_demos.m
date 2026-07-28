function run_all_demos(runHeavy)
%RUN_ALL_DEMOS Run the demonstrations in a sensible order.
%
%   run_all_demos(false) runs the fast demonstrations.
%   run_all_demos(true) also runs the paper-quality Fig. 2 calculation and
%   the longer FDTD and convergence calculations.

if nargin < 1
    runHeavy = false;
end

startup_stm();

demo02_complex_frequency_and_momentum_gaps();
demo03_ptc_bands_pwe_vs_tmm();
demo04_temporal_multilayer_tmm();
demo05_fdtd_temporal_interface();
demo08_fhs_chern_thouless_pump();
demo09_finite_ptc_order_and_phase();
demo10_temporal_domain_wall_mode();
demo11_coherent_time_interface();

if runHeavy
    demo01_reproduce_fig2_stpwe('paper');
    demo06_fdtd_spacetime_wavepacket();
    demo07_zak_phase_stpwe();
    demo12_ptc_convergence_audit();
else
    fprintf(['\nHeavy demonstrations skipped. Run run_all_demos(true), or ', ...
        'call demo01/demo06/demo07/demo12 separately.\n']);
end
end
