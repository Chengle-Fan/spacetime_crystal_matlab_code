function [upper,lower] = tl_static_omega(ka,Ls,S,C,L0,Cblock)
%TL_STATIC_OMEGA Lossless LC dispersion, including finite DC-block capacitance.
% Inputs may be compatible arrays for component-tolerance sampling.
acousticSquared = 4*sin(ka/2).^2./((Ls+2*S.*cos(ka)).*C);
sumSquared = acousticSquared+1./(L0.*C)+1./(L0.*Cblock);
productSquared = acousticSquared./(L0.*Cblock);
upperSquared = (sumSquared+sqrt(max( ...
    sumSquared.^2-4*productSquared,0)))/2;
upper = sqrt(upperSquared);
lower = sqrt(productSquared./max(upperSquared,realmin));
end
