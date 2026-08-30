function result = pwe_bands(fourier, pweCfg, mode)
%PWE_BANDS Unified spacetime PWE band solver (fixed-k or fixed-omega).
%
% Solves the plane-wave-expansion eigenproblems of a 1D, scalar,
% non-dispersive spacetime crystal in normalized units (c0 = eps0 = mu0 = 1,
% mu_r = 1 by default). The field vector R stacks the S electric-field
% harmonics E above the S magnetic-field harmonics H, so every eigenproblem
% has 2*S unknowns and yields 2*S eigenvalues per scan point.
%
%   mode = 'omega': fixed Bloch wavenumber k -> quasifrequencies omega.
%       Generalized eigenproblem  A*R = omega*B*R  with
%           K = k*eye(S) + G,  G = diag(nList*g),  W = diag(mList*Omega),
%           A = [K, -W*Cmu; -W*Ceps, K],
%           B = [zeros(S), Cmu; Ceps, zeros(S)].
%       The row blocks give K*E = (omega*I + W)*Cmu*H and
%       K*H = (omega*I + W)*Ceps*E, so R = [E; H].
%
%   mode = 'k': fixed frequency omega -> Bloch wavenumbers k.
%       Standard eigenproblem  A*R = k*R  with
%           OW = omega*eye(S) + W,
%           A  = [-G, OW*Cmu; OW*Ceps, -G],
%       and the same state R = [E; H].
%
% The raw spectrum is always preserved; folding (when enabled) maps the real
% part into the first Floquet Brillouin zone but never overwrites omegaRaw /
% kRaw. The central-harmonic weight is the fraction of |E|^2 + |H|^2 carried
% by the m = 0 temporal-DC harmonics (all n with m = 0, across both the E and
% H blocks); it is returned in BOTH modes with the shape nScan x 2S.
%
% INPUTS
%   fourier : struct from pwe_fourier (the SOLE source of Lambda/T/g/Omega
%             and of epsCoeff/muCoeff; pweCfg must carry matching values).
%   pweCfg  : struct with required fields Nspace, Mtime (nonnegative finite
%             integers, must match fourier) and, depending on mode, a scan
%             vector of finite real values:
%               'omega' -> kScan     (row of wavenumbers)
%               'k'     -> omegaScan (row of frequencies)
%             plus optional fields:
%               returnEigenvectors : scalar logical, default false.
%               fold               : scalar logical, default true.
%               Lambda, T          : optional; if present they must equal
%                                    fourier.Lambda / fourier.T.
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
%       m0Weight (nK x 2S), R and L (2S x 2S x nK, only if
%       returnEigenvectors, else []).

%==========================================================================
% Validate mode, fourier structure, and config
%==========================================================================
if nargin < 3 || isempty(mode)
    mode = 'omega';
end
mode = lower(mode);
if ~ismember(mode, {'omega', 'k'})
    error('mode must be ''omega'' (fixed-k) or ''k'' (fixed-omega).');
end
if ~isstruct(fourier)
    error('fourier must be a struct from pwe_fourier.');
end
requiredFourier = {'Nspace','Mtime','Lambda','T','g','Omega','epsCoeff','muCoeff'};
missing = requiredFourier(~isfield(fourier, requiredFourier));
if ~isempty(missing)
    error('fourier is missing required field(s): %s.', strjoin(missing, ', '));
end

% Truncation must be a nonnegative finite integer and match the sampled table.
validate_pwe_truncation(pweCfg, 'Nspace');
validate_pwe_truncation(pweCfg, 'Mtime');
Nspace = pweCfg.Nspace;
Mtime  = pweCfg.Mtime;
if Nspace ~= fourier.Nspace || Mtime ~= fourier.Mtime
    error('pweCfg.Nspace/Mtime must match the fourier truncation.');
end

% pweCfg.Lambda/T are redundant with fourier; if supplied they must agree
% (P3-01). fourier is the only data source actually used.
if isfield(pweCfg, 'Lambda') && ~isempty(pweCfg.Lambda) && ...
        pweCfg.Lambda ~= fourier.Lambda
    error('pweCfg.Lambda does not match fourier.Lambda; fourier is the sole source.');
end
if isfield(pweCfg, 'T') && ~isempty(pweCfg.T) && pweCfg.T ~= fourier.T
    error('pweCfg.T does not match fourier.T; fourier is the sole source.');
end

Lambda = fourier.Lambda;
T       = fourier.T;
g       = fourier.g;
Omega   = fourier.Omega;

returnEigenvectors = logical_opt(pweCfg, 'returnEigenvectors', false);
doFold             = logical_opt(pweCfg, 'fold', true);

%==========================================================================
% Harmonic index lists and convolution matrices
%==========================================================================
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

% Central-harmonic mask: all rows with temporal order m = 0 (any n).
m0 = (mList == 0);

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

%==========================================================================
% mode 'omega': fixed k -> quasifrequencies
%==========================================================================
if strcmp(mode, 'omega')
    kScan = validate_scan(pweCfg, 'kScan', 'mode ''omega''');
    nK = numel(kScan);

    omegaRaw = complex(zeros(2*S, nK));
    omega    = complex(zeros(2*S, nK));
    m0Weight = zeros(nK, 2*S);
    R = []; L = [];
    if returnEigenvectors
        R = complex(zeros(2*S, 2*S, nK));
        L = complex(zeros(2*S, 2*S, nK));
    end

    for ik = 1:nK
        K = kScan(ik)*eye(S) + G;
        A = [K, -W*Cmu; -W*Ceps, K];
        B = [zeros(S), Cmu; Ceps, zeros(S)];
        [Rv, Dv, Lv] = eig(A, B);
        wRaw = diag(Dv);                % raw complex spectrum, never filtered
        omegaRaw(:, ik) = wRaw;

        % Central-harmonic (m=0) weight per eigenvector.
        m0Weight(ik, :) = central_weight(Rv, S, m0);

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
%==========================================================================
% mode 'k': fixed frequency -> Bloch wavenumbers
%==========================================================================
    omegaScan = validate_scan(pweCfg, 'omegaScan', 'mode ''k''');
    nK = numel(omegaScan);

    kRaw = complex(zeros(2*S, nK));
    k    = complex(zeros(2*S, nK));
    m0Weight = zeros(nK, 2*S);
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

        % Same m=0 central-harmonic weight contract as mode 'omega'.
        m0Weight(iw, :) = central_weight(Rv, S, m0);

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
    result.m0Weight  = m0Weight;
    result.R         = R;
    result.L         = L;
end
end

% -------------------------------------------------------------------------
function weight = central_weight(Rv, S, m0)
%CENTRAL_WEIGHT Fraction of |E|^2 + |H|^2 carried by the m = 0 rows.
E = Rv(1:S, :);
H = Rv(S+1:end, :);
num = sum(abs(E(m0,:)).^2 + abs(H(m0,:)).^2, 1);
den = sum(abs(E).^2 + abs(H).^2, 1);
weight = real(num ./ max(den, eps));
end

% -------------------------------------------------------------------------
function val = logical_opt(pweCfg, fieldName, defaultVal)
%LOGICAL_OPT Scalar-logical option with strict scalar 0/1/logical validation.
if isfield(pweCfg, fieldName) && ~isempty(pweCfg.(fieldName))
    x = pweCfg.(fieldName);
    if ~islogical(x) && ~(isnumeric(x) && isscalar(x) && (x == 0 || x == 1))
        error('pweCfg.%s must be a scalar logical (true/false or 0/1).', ...
            fieldName);
    end
    val = logical(x);
else
    val = defaultVal;
end
end

% -------------------------------------------------------------------------
function scan = validate_scan(pweCfg, fieldName, modeName)
%VALIDATE_SCAN A scan vector must be nonempty with all finite real values.
if ~isfield(pweCfg, fieldName) || isempty(pweCfg.(fieldName))
    error('pweCfg.%s (row) is required for %s.', fieldName, modeName);
end
scan = pweCfg.(fieldName)(:).';
if ~isreal(scan) || ~all(isfinite(scan))
    error('pweCfg.%s must contain only finite real values.', fieldName);
end
end

% -------------------------------------------------------------------------
function ord = validate_pwe_truncation(pweCfg, fieldName)
%VALIDATE_PWE_TRUNCATION A truncation (Nspace/Mtime) must be a nonneg integer.
if ~isfield(pweCfg, fieldName)
    error('pweCfg.%s is required.', fieldName);
end
ord = pweCfg.(fieldName);
if ~isscalar(ord) || ~isreal(ord) || ~isfinite(ord) || ord < 0 || ...
        ord ~= round(ord)
    error('pweCfg.%s must be a nonnegative finite integer scalar.', fieldName);
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
