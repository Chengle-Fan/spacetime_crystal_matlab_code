function result = pwe_bands(fourier, kScan)
%PWE_BANDS Complex quasifrequency bands of a photonic time crystal.
%
% For each real Bloch wavenumber k, solve the spacetime PWE generalized
% eigenproblem for the field harmonics [E;H]:
%
%   [kI  -W        ] [E] = omega [0         I] [E]
%   [-W*Cepsilon kI] [H]         [Cepsilon 0] [H],
%
% where W=diag(m*Omega), epsilon is periodic in time and mu_r=1. The many raw
% roots contain Floquet copies separated by Omega. One representative of each
% of the two physical branches is selected in a shifted Floquet interval, then
% folded into the standard first zone [-Omega/2,Omega/2).
%
% INPUTS
%   fourier : output of pwe_fourier
%   kScan   : nonempty vector of finite real wavenumbers
%
% OUTPUT
%   result.k     : 1 x nK scan vector
%   result.omega : 2 x nK complex quasifrequencies in the first Floquet zone

if nargin ~= 2 || ~isstruct(fourier) || ~isscalar(fourier)
    error('Use pwe_bands(fourier,kScan).');
end
required = {'T','Omega','Mtime','mList','Cepsilon'};
missing = required(~isfield(fourier,required));
if ~isempty(missing)
    error('fourier is missing field(s): %s.',strjoin(missing,', '));
end
if ~isnumeric(kScan) || isempty(kScan) || ~isvector(kScan) || ...
        ~isreal(kScan) || any(~isfinite(kScan))
    error('kScan must be a nonempty vector of finite real wavenumbers.');
end
kScan = kScan(:).';

T = fourier.T;
Omega = fourier.Omega;
Mtime = fourier.Mtime;
if ~isnumeric(T) || ~isscalar(T) || ~isreal(T) || ~isfinite(T) || T <= 0 || ...
        ~isnumeric(Omega) || ~isscalar(Omega) || ~isreal(Omega) || ...
        ~isfinite(Omega) || abs(Omega-2*pi/T) > 64*eps(max(1,Omega)) || ...
        ~isnumeric(Mtime) || ~isscalar(Mtime) || Mtime < 0 || ...
        Mtime ~= round(Mtime)
    error('The Fourier geometry is inconsistent; rerun pwe_fourier.');
end

mList = (-Mtime:Mtime).';
S = numel(mList);
Cepsilon = fourier.Cepsilon;
if ~isequal(fourier.mList,mList) || ~isnumeric(Cepsilon) || ...
        ~isequal(size(Cepsilon),[S S]) || any(~isfinite(Cepsilon),'all')
    error('The Fourier harmonic ordering or convolution matrix is invalid.');
end

I = eye(S);
Z = zeros(S);
W = diag(mList*Omega);
B = [Z,I;Cepsilon,Z];
omega = complex(zeros(2,numel(kScan)));

for ik = 1:numel(kScan)
    K = kScan(ik)*I;
    A = [K,-W;-W*Cepsilon,K];
    raw = eig(A,B);

    % Floquet copies differ by integer Omega. Selecting the two roots nearest
    % Omega/4 is equivalent to choosing the shifted representative interval
    % [-Omega/4,3*Omega/4), whose boundaries avoid the physical zone center
    % and edge where momentum gaps normally form. No eigenvector weight or
    % external reference solver is used.
    [~,order] = sort(abs(real(raw)-Omega/4));
    selected = raw(order(1:2));
    omega(:,ik) = fold_frequency(selected,Omega);
end

result.k = kScan;
result.omega = omega;
end

% -------------------------------------------------------------------------
function folded = fold_frequency(omega,Omega)
folded = mod(real(omega)+Omega/2,Omega)-Omega/2 + 1i*imag(omega);
end
