function result = pwe_bands(fourier, pweCfg, mode)
%PWE_BANDS Unified spacetime PWE band solver (fixed-k or fixed-omega).
%
% Solves the plane-wave-expansion eigenproblems of a 1D, scalar,
% non-dispersive spacetime crystal in normalized units (c0 = eps0 = mu0 = 1,
% mu_r = 1 by default). The field vector R stacks the S displacement-flux
% harmonics D above the S magnetic-flux harmonics B, so every eigenproblem
% has 2*S unknowns and yields 2*S eigenvalues per scan point.
%
%   mode = 'omega': fixed Bloch wavenumber k -> quasifrequencies omega.
%       Generalized eigenproblem  A*R = omega*B*R  (Park-Min Eq. 6) with
%           K = k*eye(S) + G,  G = diag(nList*g),  W = diag(mList*Omega),
%           A = [K, -W*Cmu; -W*Ceps, K],
%           B = [zeros(S), Cmu; Ceps, zeros(S)].
%
%   mode = 'k': fixed frequency omega -> Bloch wavenumbers k.
%       Standard eigenproblem  A*R = k*R  (Park-Min Eq. 7) with
%           OW = omega*eye(S) + W,
%           A  = [-G, OW*Cmu; OW*Ceps, -G].
%
% The raw spectrum is always preserved; folding (when enabled) maps the real
% part into the first Floquet Brillouin zone but never overwrites omegaRaw /
% kRaw.
%
% INPUTS
%   fourier : struct from pwe_fourier (provides epsCoeff/muCoeff, geometry).
%   pweCfg  : struct with required fields Lambda, T, Nspace, Mtime (must
%             match fourier) and, depending on mode, a scan vector:
%               'omega' -> kScan     (row of wavenumbers)
%               'k'     -> omegaScan (row of frequencies)
%             plus optional fields:
%               returnEigenvectors : logical, default false.
%               fold               : logical, default true.
%   mode    : 'omega' or 'k'.
%
% OUTPUT (struct)
%   result.mode/.Lambda/.T/.g/.Omega/.Nspace/.Mtime/.S/.nList/.mList.
%   mode 'omega':
%       kScan (1 x nK), omegaRaw (2S x nK, unfiltered), omega (2S x nK,
%       folded), m0Weight (nK x 2S central-harmonic weight), R and L
%       (2S x 2S x nK, only if returnEigenvectors, else []).
%   mode 'k':
%       omegaScan (1 x nK), kRaw (2S x nK, unfiltered), k (2S x nK, folded),
%       R and L (2S x 2S x nK, only if returnEigenvectors, else []).

if nargin < 3 || isempty(mode)
    mode = 'omega';
end
mode = lower(mode);
if ~ismember(mode, {'omega', 'k'})
    error('mode must be ''omega'' (fixed-k) or ''k'' (fixed-omega).');
end

% Truncation must match the sampled coefficient table.
if ~isfield(pweCfg, 'Nspace') || ~isfield(pweCfg, 'Mtime')
    error('pweCfg.Nspace and pweCfg.Mtime are required.');
end
Nspace = pweCfg.Nspace;
Mtime  = pweCfg.Mtime;
if Nspace ~= fourier.Nspace || Mtime ~= fourier.Mtime
    error('pweCfg.Nspace/Mtime must match the fourier truncation.');
end

Lambda = fourier.Lambda;
T       = fourier.T;
g       = fourier.g;
Omega   = fourier.Omega;

returnEigenvectors = isfield(pweCfg, 'returnEigenvectors') ...
    && logical(pweCfg.returnEigenvectors);
if isfield(pweCfg, 'fold')
    doFold = logical(pweCfg.fold);
else
    doFold = true;
end

% Harmonic index lists (column) and convolution matrices.
[NN, MM] = ndgrid(-Nspace:Nspace, -Mtime:Mtime);
nList = NN(:);                          % column
mList = MM(:);                          % column
S = numel(nList);

Ceps = complex(zeros(S));
Cmu  = complex(zeros(S));
for row = 1:S
    for col = 1:S
        dm = mList(row) - mList(col);
        dn = nList(row) - nList(col);
        Ceps(row,col) = fourier.epsCoeff(dm, dn);
        Cmu(row,col)  = fourier.muCoeff(dm, dn);
    end
end

G = diag(nList*g);
W = diag(mList*Omega);

% Common header fields.
result.mode   = mode;
result.Lambda = Lambda;
result.T      = T;
result.g      = g;
result.Omega  = Omega;
result.Nspace = Nspace;
result.Mtime  = Mtime;
result.S      = S;
result.nList  = nList;
result.mList  = mList;

if strcmp(mode, 'omega')
    if ~isfield(pweCfg, 'kScan') || isempty(pweCfg.kScan)
        error('pweCfg.kScan (row) is required for mode ''omega''.');
    end
    kScan = pweCfg.kScan(:).';
    nK = numel(kScan);

    omegaRaw = complex(zeros(2*S, nK));
    omega    = complex(zeros(2*S, nK));
    m0Weight = zeros(nK, 2*S);
    R = []; L = [];
    if returnEigenvectors
        R = complex(zeros(2*S, 2*S, nK));
        L = complex(zeros(2*S, 2*S, nK));
    end

    m0 = (mList == 0);                  % logical mask of the m=0 rows
    for ik = 1:nK
        K = kScan(ik)*eye(S) + G;
        A = [K, -W*Cmu; -W*Ceps, K];
        B = [zeros(S), Cmu; Ceps, zeros(S)];
        [Rv, Dv, Lv] = eig(A, B);
        wRaw = diag(Dv);                % raw complex spectrum, never filtered
        omegaRaw(:, ik) = wRaw;

        % Central-harmonic (m=0) weight per eigenvector.
        E = Rv(1:S, :);
        H = Rv(S+1:end, :);
        num = sum(abs(E(m0,:)).^2 + abs(H(m0,:)).^2, 1);
        den = sum(abs(E).^2 + abs(H).^2, 1);
        m0Weight(ik, :) = real(num ./ max(den, eps));

        if doFold
            omega(:, ik) = fold_frequency(wRaw, Omega);
        else
            omega(:, ik) = wRaw;
        end
        if returnEigenvectors
            R(:, :, ik) = Rv;
            L(:, :, ik) = Lv;
        end
    end

    result.kScan    = kScan;
    result.omegaRaw = omegaRaw;
    result.omega    = omega;
    result.m0Weight = m0Weight;
    result.R        = R;
    result.L        = L;
else
    if ~isfield(pweCfg, 'omegaScan') || isempty(pweCfg.omegaScan)
        error('pweCfg.omegaScan (row) is required for mode ''k''.');
    end
    omegaScan = pweCfg.omegaScan(:).';
    nK = numel(omegaScan);

    kRaw = complex(zeros(2*S, nK));
    k    = complex(zeros(2*S, nK));
    R = []; L = [];
    if returnEigenvectors
        R = complex(zeros(2*S, 2*S, nK));
        L = complex(zeros(2*S, 2*S, nK));
    end

    for iw = 1:nK
        OW = omegaScan(iw)*eye(S) + W;
        A = [-G, OW*Cmu; OW*Ceps, -G];
        [Rv, Dv, Lv] = eig(A);
        kkRaw = diag(Dv);               % raw complex spectrum, never filtered
        kRaw(:, iw) = kkRaw;

        if doFold
            k(:, iw) = fold_wavenumber(kkRaw, g);
        else
            k(:, iw) = kkRaw;
        end
        if returnEigenvectors
            R(:, :, iw) = Rv;
            L(:, :, iw) = Lv;
        end
    end

    result.omegaScan = omegaScan;
    result.kRaw      = kRaw;
    result.k         = k;
    result.R         = R;
    result.L         = L;
end
end

% -------------------------------------------------------------------------
function omegaFolded = fold_frequency(omega, Omega)
%FOLD_FREQUENCY Fold Re(omega) into [-Omega/2, Omega/2).
omegaFolded = mod(real(omega) + Omega/2, Omega) - Omega/2 ...
    + 1i*imag(omega);
end

% -------------------------------------------------------------------------
function kFolded = fold_wavenumber(k, g)
%FOLD_WAVENUMBER Fold Re(k) into [-g/2, g/2).
kFolded = mod(real(k) + g/2, g) - g/2 + 1i*imag(k);
end
