function sys = stpwe_build_system(epsCoeff, muCoeff, Nspace, Mtime, g, Omega)
%STPWE_BUILD_SYSTEM Build the convolution matrices for 1-D ST-PWE.
%
% epsCoeff(m,n) and muCoeff(m,n) return coefficients in
%   p(x,t) = sum_{m,n} p_{m,n} exp(i*n*g*x - i*m*Omega*t).
%
% Nspace keeps n=-Nspace:Nspace and Mtime keeps m=-Mtime:Mtime.

if nargin < 2 || isempty(muCoeff)
    muCoeff = @(m,n) double(m == 0 && n == 0);
end

[NN, MM] = ndgrid(-Nspace:Nspace, -Mtime:Mtime);
nList = NN(:);
mList = MM(:);
S = numel(nList);

Ceps = complex(zeros(S));
Cmu = complex(zeros(S));
for row = 1:S
    for col = 1:S
        dm = mList(row) - mList(col);
        dn = nList(row) - nList(col);
        Ceps(row,col) = epsCoeff(dm,dn);
        Cmu(row,col) = muCoeff(dm,dn);
    end
end

sys.Nspace = Nspace;
sys.Mtime = Mtime;
sys.nList = nList;
sys.mList = mList;
sys.S = S;
sys.g = g;
sys.Omega = Omega;
sys.Ceps = Ceps;
sys.Cmu = Cmu;
sys.G = diag(nList*g);
sys.W = diag(mList*Omega);
sys.I = eye(S);
sys.Z = zeros(S);
sys.Bomega = [sys.Z, Cmu; Ceps, sys.Z];
end
