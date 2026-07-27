function sol = stpwe_solve_omega(sys, k)
%STPWE_SOLVE_OMEGA Fixed-k ST-PWE eigenproblem (Park-Min Eq. 6).
%
% The generalized problem is A*R = omega*B*R. Left eigenvectors obey
% L'*A = omega*L'*B and are normalized so L_i'*B*R_i = 1 whenever the
% eigenpair is not self-orthogonal.

K = k*sys.I + sys.G;
A = [K, -sys.W*sys.Cmu; ...
    -sys.W*sys.Ceps, K];
B = sys.Bomega;

[R, D, L] = eig(A, B);
omega = diag(D);

for q = 1:numel(omega)
    overlap = L(:,q)'*B*R(:,q);
    if isfinite(overlap) && abs(overlap) > 1e-12
        R(:,q) = R(:,q)/overlap;
    else
        nr = norm(R(:,q));
        nl = norm(L(:,q));
        if nr > 0
            R(:,q) = R(:,q)/nr;
        end
        if nl > 0
            L(:,q) = L(:,q)/nl;
        end
    end
end

S = sys.S;
mask0 = (sys.mList == 0);
E = R(1:S,:);
H = R(S+1:end,:);
numerator = sum(abs(E(mask0,:)).^2 + abs(H(mask0,:)).^2, 1);
denominator = sum(abs(E).^2 + abs(H).^2, 1);
m0Weight = real(numerator./max(denominator, eps));

sol.k = k;
sol.omega = omega;
sol.R = R;
sol.L = L;
sol.A = A;
sol.B = B;
sol.m0Weight = m0Weight(:);
end
