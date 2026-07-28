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
    % At k = 0 the matrix K = diag(n*g) has a zero on the n = 0
    % diagonal entry, making A = [K, Z; Z, K] singular.  The
    % generalized eigenvalue problem then produces two ω ≈ 0
    % eigenvalues (one from the E-block nullspace, one from the
    % H-block).  Only one is physical — the Γ-point crossing of
    % the forward-traveling wave.  The duplicate shifts all
    % subsequent band indices at k = 0, creating an artificial
    % "spike" where the plotted row jumps between distinct bands.
    idxZero = abs(omega) < 1e-8;
    if sum(idxZero) > 1
        omega = omega(~idxZero);       % remove all near-zero
        omega = [0; omega];            % keep exactly one
    end
    take = min(nBands, numel(omega));
    % Since Lambda=2*pi/g, fBar=omega*Lambda/(2*pi*c0).
    fBands(1:take,ik) = omega(1:take)/(g*c0);
end
end
