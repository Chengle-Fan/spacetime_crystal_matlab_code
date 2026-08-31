function fourier = tl_pwe_fourier(model,pweCfg)
%TL_PWE_FOURIER Fourier coefficients for a time-modulated TL state model.
%
% The conservative Q/Phi formulation makes the periodic coefficient enter
% through 1/C(t).  Only k-independent temporal data are computed here;
% tl_pwe_bands later combines them with Bloch incidence matrices.

if nargin ~= 2 || ~isstruct(model) || ~isscalar(model) || ...
        ~isstruct(pweCfg) || ~isscalar(pweCfg)
    error('Use tl_pwe_fourier(model,pweCfg).');
end
requiredModel = {'modulation','bulk','functions','snapshot'};
missing = requiredModel(~isfield(model,requiredModel));
if ~isempty(missing)
    error('model is missing field(s): %s.',strjoin(missing,', '));
end
Mtime = read_integer(pweCfg,'Mtime',0);
Nt = read_integer(pweCfg,'Nt',1);
if Nt < 4*Mtime+1
    error('pweCfg.Nt must satisfy Nt >= 4*Mtime+1 = %d.',4*Mtime+1);
end
T = model.modulation.period;
Omega = model.modulation.OmegaRadPerSec;
t = (0:Nt-1)*T/Nt;
C = model.functions.capacitanceBulk(t);
if ~isequal(size(C),size(t)) || ~isreal(C) || any(~isfinite(C)) || any(C <= 0)
    error('The bulk capacitance callback returned invalid values.');
end
inverseC = 1./C;

differenceOrders = -2*Mtime:2*Mtime;
inverseCCoeff = complex(zeros(size(differenceOrders)));
for index = 1:numel(differenceOrders)
    order = differenceOrders(index);
    inverseCCoeff(index) = mean(inverseC.*exp(1i*order*Omega*t));
end

mList = (-Mtime:Mtime).';
dm = mList-mList.';
Cinverse = inverseCCoeff(dm+2*Mtime+1);
scale = max(1,norm(Cinverse,'fro'));
if norm(Cinverse-Cinverse','fro')/scale > 1e-11
    error('The inverse-capacitance convolution matrix is not Hermitian.');
end
Cinverse = (Cinverse+Cinverse')/2;
[~,flag] = chol(Cinverse);
if flag ~= 0
    error(['The inverse-capacitance convolution matrix is not positive ' ...
        'definite; increase Nt or inspect the capacitance waveform.']);
end

fourier.T = T;
fourier.Omega = Omega;
fourier.Mtime = Mtime;
fourier.Nt = Nt;
fourier.mList = mList;
fourier.differenceOrders = differenceOrders;
fourier.inverseCCoeff = inverseCCoeff;
fourier.Cinverse = Cinverse;
fourier.capacitanceMinimumSampled = min(C);
fourier.capacitanceMaximumSampled = max(C);
fourier.stateDimension = model.bulk.stateDimension;
fourier.modelKind = model.kind;
fourier.modelSnapshot = model.snapshot;
end

% -------------------------------------------------------------------------
function value = read_integer(cfg,name,minimumValue)
if ~isfield(cfg,name)
    error('pweCfg.%s is required.',name);
end
value = cfg.(name);
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value ~= round(value) || value < minimumValue
    error('pweCfg.%s must be an integer not smaller than %d.', ...
        name,minimumValue);
end
end
