function spectrum = temporal_tmm_spectrum(omegaInitial, ...
    epsInitial, muInitial, epsSlabs, muSlabs, durations, ...
    epsFinal, muFinal)
%TEMPORAL_TMM_SPECTRUM Forward/backward electric-field coefficients.

omegaInitial = omegaInitial(:).';
forward = complex(zeros(size(omegaInitial)));
backward = complex(zeros(size(omegaInitial)));
determinant = complex(zeros(size(omegaInitial)));

for q = 1:numel(omegaInitial)
    TM = temporal_multilayer_tmm(omegaInitial(q), ...
        epsInitial, muInitial, epsSlabs, muSlabs, durations, ...
        epsFinal, muFinal);
    output = TM*[1;0];
    forward(q) = output(1);
    backward(q) = output(2);
    determinant(q) = det(TM);
end

spectrum.omegaInitial = omegaInitial;
spectrum.forward = forward;
spectrum.backward = backward;
spectrum.detTM = determinant;
end
