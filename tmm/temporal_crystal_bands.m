function bands = temporal_crystal_bands(kValues, epsSequence, ...
    muSequence, durations)
%TEMPORAL_CRYSTAL_BANDS Floquet quasifrequencies from a temporal monodromy.
%
% If lambda=exp(-i*omegaF*T), then omegaF=i*log(lambda)/T. The principal
% logarithm places Re(omegaF) in the first temporal Brillouin zone.

kValues = kValues(:).';
T = sum(durations);
omegaF = complex(zeros(2,numel(kValues)));
lambda = complex(zeros(2,numel(kValues)));
halfTrace = complex(zeros(size(kValues)));

for ik = 1:numel(kValues)
    U = temporal_crystal_monodromy(kValues(ik), epsSequence, ...
        muSequence, durations);
    lam = eig(U);
    omega = 1i*log(lam)/T;
    [~, order] = sortrows([real(omega), imag(omega)], [1 2]);
    omegaF(:,ik) = omega(order);
    lambda(:,ik) = lam(order);
    halfTrace(ik) = trace(U)/2;
end

bands.k = kValues;
bands.omegaF = omegaF;
bands.lambda = lambda;
bands.halfTrace = halfTrace;
bands.T = T;
bands.Omega = 2*pi/T;
end
