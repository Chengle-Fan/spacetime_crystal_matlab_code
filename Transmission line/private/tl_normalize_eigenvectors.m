function [V,W,conditionNumber] = tl_normalize_eigenvectors(V,W)
%TL_NORMALIZE_EIGENVECTORS Unit left/right vectors and eigenvalue sensitivity.
nv = vecnorm(V);
nw = vecnorm(W);
if any(~isfinite([nv nw])) || any([nv nw] == 0)
    error('The eigensolver returned an invalid left or right eigenvector.');
end
V = V./nv;
W = W./nw;
conditionNumber = (1./max(abs(sum(conj(W).*V,1)),realmin)).';
end
