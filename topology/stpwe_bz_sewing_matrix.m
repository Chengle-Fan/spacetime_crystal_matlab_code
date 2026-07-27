function sewing = stpwe_bz_sewing_matrix(sys)
%STPWE_BZ_SEWING_MATRIX Join k=-g/2 to the equivalent state at k=+g/2.
%
% With E=exp(i*k*x)*sum_n u_n exp(i*n*g*x), the same physical state obeys
% u_n(k+g)=u_{n+1}(k). The finite Fourier cutoff drops only the outermost
% coefficient; convergence requires negligible weight at that cutoff.

S = sys.S;
shift = zeros(S);
for row = 1:S
    for col = 1:S
        if sys.mList(row) == sys.mList(col) ...
                && sys.nList(col) == sys.nList(row)+1
            shift(row,col) = 1;
        end
    end
end
sewing = blkdiag(shift,shift);
end
