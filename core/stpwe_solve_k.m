function sol = stpwe_solve_k(sys, omega)
%STPWE_SOLVE_K Fixed-omega ST-PWE eigenproblem (Park-Min Eq. 7).

OW = omega*sys.I + sys.W;
A = [-sys.G, OW*sys.Cmu; ...
    OW*sys.Ceps, -sys.G];

[R, D, L] = eig(A);
k = diag(D);

for q = 1:numel(k)
    overlap = L(:,q)'*R(:,q);
    if isfinite(overlap) && abs(overlap) > 1e-12
        R(:,q) = R(:,q)/overlap;
    end
end

sol.omega = omega;
sol.k = k;
sol.R = R;
sol.L = L;
sol.A = A;
end
