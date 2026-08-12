function [beta,ZB,isPass,info] = otr_bloch_dispersion(f,state,model,p)
%OTR_BLOCH_DISPERSION Bloch wavenumber and impedance (SI Eqs. S2--S3).
%   [BETA,ZB,ISPASS,INFO] = OTR_BLOCH_DISPERSION(F,STATE,MODEL,P)
%   evaluates a periodic chain of OTR_ABCD_UNIT_CELL cells. BETA has units
%   rad/m and ZB ohms. First-zone BETA is chosen with Im(BETA)<=0 for the
%   exp(+j*omega*t) / forward-decay convention used in the paper.
%
%   For the ideal loaded line this is strictly
%      cos(beta*d)=cos(theta)-b*sin(theta)/2                 (S2)
%      ZB = +/- B*Z0/sqrt(A^2-1)                            (S3)
%   where A,B are the normalized entries of Eq. S1. The sign of ZB is
%   selected for positive real part in a lossless pass band. REALISTIC
%   permits complex beta/ZB and uses the same periodic ABCD definition.
%
%   isPass is true when |(A+D)/2|<=1 for a lossless cell; for a lossy cell
%   it flags points with small Bloch attenuation relative to phase.

if nargin < 2 || isempty(state), state = 'on'; end
if nargin < 3 || isempty(model), model = 'ideal'; end
if nargin < 4 || isempty(p), p = otr_parameters(); end
[~,a] = otr_abcd_unit_cell(f,state,model,p);
Mn = a.normalizedMatrix;
fv = f(:).';
n = numel(fv);
traceHalf = complex(zeros(1,n));
betaV = complex(zeros(1,n));
zbV = complex(zeros(1,n));
detV = complex(zeros(1,n));

for k = 1:n
    A = Mn(1,1,k); B = Mn(1,2,k); D = Mn(2,2,k);
    traceHalf(k) = 0.5*(A+D);
    detV(k) = det(Mn(:,:,k));
    q = acos(traceHalf(k))/p.d;
    if imag(q) > 0, q = -q; end
    if real(q) < 0, q = -q; end
    betaV(k) = q;

    % Eq. S3 is printed as B*Z0/sqrt(A^2-1).  In a pass band both B and
    % sqrt(A^2-1) are imaginary; choose the square-root branch that makes
    % the forward Bloch impedance positive-real.
    den = sqrt(A.^2-1);
    if abs(den) <= 64*eps(max(1,abs(A)))
        z = NaN;
    else
        z = B*p.Z0/den;                    % S3
        if real(z) < 0 || (abs(real(z)) < 1e-12*max(1,abs(z)) && imag(z) < 0)
            z = -z;
        end
    end
    zbV(k) = z;
end

tol = 256*eps;
lossless = max(abs(real(a.normalizedAdmittance(:)))) < 1e-11;
if lossless
    passV = abs(real(traceHalf)) <= 1+tol & abs(imag(traceHalf)) <= 1e-10;
else
    passV = abs(imag(betaV))*p.d <= max(2e-2,0.05*abs(real(betaV))*p.d);
end

shape = size(f);
beta = reshape(betaV,shape);
ZB = reshape(zbV,shape);
isPass = reshape(passV,shape);
info = struct('frequency',f,'traceHalf',reshape(traceHalf,shape), ...
    'determinant',reshape(detV,shape),'theta',a.theta,'b',a.b, ...
    'normalizedAdmittance',a.normalizedAdmittance, ...
    'state',a.state,'model',a.model,'firstBrillouinEdge',pi/p.d);
end
