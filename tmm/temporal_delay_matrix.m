function DM = temporal_delay_matrix( ...
    omegaInitial, nInitial, nSlab, duration)
%TEMPORAL_DELAY_MATRIX Phase accumulation inside one temporal slab.
%
% The code follows the exp(+i*omega*t) convention used by Ramaccia et al.
% The conserved wavevector implies omegaSlab=(nInitial/nSlab)*omegaInitial.

omegaSlab = (nInitial/nSlab)*omegaInitial;
phase = omegaSlab*duration;
DM = diag([exp(1i*phase), exp(-1i*phase)]);
end
