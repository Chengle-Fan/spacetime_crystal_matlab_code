function fBands = stpwe_static_bands(kValues, epsCoeff0, Nspace, g, c0, nBands)
%STPWE_STATIC_BANDS Static 1-D photonic bands in the same first-order form.

n = (-Nspace:Nspace).';
Ns = numel(n);
G = diag(n*g);
I = eye(Ns);
Z = zeros(Ns);
Ceps = complex(zeros(Ns));

for row = 1:Ns
    for col = 1:Ns
        Ceps(row,col) = epsCoeff0(n(row)-n(col));
    end
end

B = [Z, I; Ceps, Z];
fBands = nan(nBands, numel(kValues));
for ik = 1:numel(kValues)
    K = kValues(ik)*I + G;
    omega = eig([K, Z; Z, K], B);
    good = isfinite(omega) & abs(imag(omega)) < 1e-8 & real(omega) >= -1e-9;
    omega = sort(real(omega(good)));
    take = min(nBands, numel(omega));
    % Since Lambda=2*pi/g, fBar=omega*Lambda/(2*pi*c0).
    fBands(1:take,ik) = omega(1:take)/(g*c0);
end
end
