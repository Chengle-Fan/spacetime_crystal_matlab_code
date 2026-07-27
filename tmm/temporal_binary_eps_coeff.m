function coeff = temporal_binary_eps_coeff(m, epsA, epsB, dutyA)
%TEMPORAL_BINARY_EPS_COEFF Fourier coefficient of a binary time crystal.
%
% epsilon(t)=epsA for 0<=t<dutyA*T and epsB otherwise, with
% epsilon(t)=sum_m epsilon_m exp(-i*m*Omega*t).

if m == 0
    coeff = dutyA*epsA + (1-dutyA)*epsB;
else
    coeff = (epsA-epsB)*(exp(1i*2*pi*m*dutyA)-1) ...
        /(1i*2*pi*m);
end
end
