function stm_run_examples(runHeavy)
%STM_RUN_EXAMPLES Run the toolbox examples in a sensible order.
%
%   stm_run_examples(false) runs the fast examples.
%   stm_run_examples(true) also runs the paper-quality Fig. 2 calculation
%   and the longer FDTD simulation.

if nargin < 1
    runHeavy = false;
end

stm_init();

example_complex_gaps();
example_ptc_pwe_vs_tmm();
example_tmm_multilayer();
example_fdtd_interface();
example_chern_pump();

if runHeavy
    example_floquet_bands_and_fields('paper');
    example_fdtd_wavepacket();
    example_zak_phase();
else
    fprintf(['\nHeavy examples skipped. Run stm_run_examples(true), or ', ...
        'call example_floquet_bands_and_fields / example_fdtd_wavepacket / ', ...
        'example_zak_phase separately.\n']);
end
end
