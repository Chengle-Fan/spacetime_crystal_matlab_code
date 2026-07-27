function omegaFolded = stpwe_fold_frequency(omega, Omega)
%STPWE_FOLD_FREQUENCY Fold Re(omega) into [-Omega/2,Omega/2).

omegaFolded = mod(real(omega) + Omega/2, Omega) - Omega/2 ...
    + 1i*imag(omega);
end
