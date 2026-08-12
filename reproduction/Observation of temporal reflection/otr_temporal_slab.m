function [S,paths] = otr_temporal_slab(Zouter,Zslab,tau,omegaSlab,varargin)
%OTR_TEMPORAL_SLAB Two-interface temporal slab and its four paths.
%   [S,PATHS] = OTR_TEMPORAL_SLAB(ZOUTER,ZSLAB,TAU,OMEGASLAB) models the
%   experimental OFF-ON-OFF slab. The first facet uses S14/S16 (charge
%   conservation), the second S15/S17 (voltage conservation). OMEGASLAB is
%   the angular frequency during the slab and TAU its duration in seconds.
%   Scalar impedances may accompany an arbitrary equal-size frequency array.
%
%   S = Joff*diag(exp(+j*phi),exp(-j*phi))*Jon, phi=omegaSlab*tau.
%   PATHS contains the four elementary amplitudes cited below SI Eq. S17:
%     TT = Ton*Toff*exp(+j phi), RR = Ron*Roff*exp(-j phi),
%     RT = Ron*Toff*exp(-j phi), TR = Ton*Roff*exp(+j phi).
%   Thus S(1,1)=TT+RR (time refraction) and S(2,1)=RT+TR (reflection).
%
%   Optional 'FirstLaw' and 'SecondLaw' override the default microscopic
%   laws, which also supports the inverted ON-OFF-ON control experiment.

if ~(isnumeric(tau) && isreal(tau) && isscalar(tau) && isfinite(tau) && tau>=0)
    error('otr:slab:Duration','tau must be a nonnegative finite scalar [s].');
end
if ~(isnumeric(omegaSlab) && isreal(omegaSlab) && ~isempty(omegaSlab) && ...
        all(isfinite(omegaSlab(:))) && all(omegaSlab(:)>=0))
    error('otr:slab:Frequency','omegaSlab must be nonnegative finite [rad/s].');
end
firstLaw = 'charge'; secondLaw = 'voltage';
if mod(numel(varargin),2) ~= 0
    error('otr:slab:Options','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = lower(char(varargin{k})); value = varargin{k+1};
    switch name
        case 'firstlaw', firstLaw = value;
        case 'secondlaw', secondLaw = value;
        otherwise, error('otr:slab:Options','Unknown option ''%s''.',name);
    end
end

[J1,T1,R1,i1] = otr_temporal_interface(Zouter,Zslab,firstLaw);
[J2,T2,R2,i2] = otr_temporal_interface(Zslab,Zouter,secondLaw);
% Establish one common frequency shape. Each quantity may be scalar, but
% nonscalar inputs must describe the same frequency samples.
arrays = {omegaSlab,Zouter,Zslab};
sz = size(omegaSlab);
hasArray = ~isscalar(omegaSlab);
for k = 1:numel(arrays)
    if ~isscalar(arrays{k})
        if ~hasArray, sz = size(arrays{k}); hasArray = true;
        elseif ~isequal(size(arrays{k}),sz)
            error('otr:slab:Size', ...
                'Nonscalar impedances and omegaSlab must have equal sizes.');
        end
    end
end
if ~hasArray, sz = size(omegaSlab); end
omegaSlab = omegaSlab+zeros(sz);
phi = omegaSlab*tau;
eplus = exp(1i*phi); eminus = exp(-1i*phi);

T1 = T1+zeros(sz); R1 = R1+zeros(sz);
T2 = T2+zeros(sz); R2 = R2+zeros(sz);
TT = T1.*T2.*eplus;
RR = R1.*R2.*eminus;
RT = R1.*T2.*eminus;
TR = T1.*R2.*eplus;
totalForward = TT+RR;
totalReflected = RT+TR;

n = numel(phi);
S = complex(zeros(2,2,n));
for k = 1:n
    P = diag([eplus(k),eminus(k)]);
    m1 = J1(:,:,min(k,size(J1,3)));
    m2 = J2(:,:,min(k,size(J2,3)));
    S(:,:,k) = m2*P*m1;
end
if isscalar(phi), S = S(:,:,1); end

paths = struct('phase',phi,'TT',TT,'RR',RR,'RT',RT,'TR',TR, ...
    'timeRefracted',totalForward,'timeReflected',totalReflected, ...
    'firstMatrix',J1,'secondMatrix',J2,'firstInterface',i1, ...
    'secondInterface',i2,'identityResidualForward', ...
    reshape(S(1,1,:),sz)-totalForward,'identityResidualReflected', ...
    reshape(S(2,1,:),sz)-totalReflected);
end
