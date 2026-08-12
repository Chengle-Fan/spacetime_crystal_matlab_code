function p = otr_parameters(varargin)
%OTR_PARAMETERS Parameters of the switched transmission-line experiment.
%   P = OTR_PARAMETERS() returns the article/Methods/SI values in SI units.
%   P = OTR_PARAMETERS('name',value,...) overrides selected scalar fields.
%   P = OTR_PARAMETERS(S) applies all fields of scalar structure S.
%
%   参数出处 / Sources
%   ------------------
%   Methods: Z0, d, epsEff, losses, ideal-switch resistance, series-RLC
%   load and parasitic capacitive coupling.  Main text and SI: 30 cells,
%   theta(100 MHz)=72 deg, b(100 MHz)=2.753 and 3 ns switching time.
%
%   No toolbox is required.  MATLAB R2019b and later are supported.

p = struct();
p.c0 = 299792458;               % vacuum light speed [m/s]
p.d = 0.208;                    % physical unit-cell length [m]
p.Z0 = 50;                      % unloaded TL characteristic impedance [ohm]
p.epsEff = 8.36;                % effective relative permittivity
p.Cload = 82e-12;               % switched shunt load [F]
p.switchRon = 1;                % switch ON resistance [ohm]
p.Rpar = 4;                     % load parasitic series resistance [ohm]
p.Lpar = 8e-9;                  % grounding-via series inductance [H]
p.nCells = 30;                  % fabricated sample cell count
p.Cnearest = 2e-12;             % nearest-neighbour coupling [F]
p.CnextNearest = 1e-12;         % second-nearest-neighbour coupling [F]
p.theta100 = 72*pi/180;         % design electrical length at 100 MHz [rad]
p.b100 = 2.753;                 % normalized shunt susceptance at 100 MHz
p.switchingTime = 3e-9;         % measured 10--90 % switching time [s]
p.fReference = 100e6;           % SI design reference [Hz]
p.alpha100 = 0.5;               % ADS TLINP attenuation A [dB/m] at 100 MHz
p.lossTangent = 0.0019;         % substrate loss tangent

% These derived quantities expose both the reported 72-degree design phase
% model and the Methods epsEff microstrip estimate.  S1--S3 use thetaScale;
% the physical L'/C' values use epsEff, matching the MNA implementation.
p.phaseVelocity = 2*pi*p.fReference*p.d/p.theta100;
p.phaseVelocityMethods = p.c0/sqrt(p.epsEff);
p.ClinePerLength = 1/(p.Z0*p.phaseVelocityMethods);
p.LlinePerLength = p.Z0/p.phaseVelocityMethods;
p.ClineCell = p.ClinePerLength*p.d;
p.LlineCell = p.LlinePerLength*p.d;

if nargin == 1 && isstruct(varargin{1})
    p = localApplyStruct(p,varargin{1});
elseif nargin > 0
    if mod(nargin,2) ~= 0
        error('otr:parameters:NameValue', ...
            'Overrides must be supplied as name/value pairs.');
    end
    for k = 1:2:nargin
        name = varargin{k};
        if ~(ischar(name) || (isstring(name) && isscalar(name)))
            error('otr:parameters:BadName','Parameter names must be text.');
        end
        p = localSet(p,char(name),varargin{k+1});
    end
end

localValidate(p);

% Recompute quantities depending on overridable primitive fields.
p.phaseVelocity = 2*pi*p.fReference*p.d/p.theta100;
p.phaseVelocityMethods = p.c0/sqrt(p.epsEff);
p.ClinePerLength = 1/(p.Z0*p.phaseVelocityMethods);
p.LlinePerLength = p.Z0/p.phaseVelocityMethods;
p.ClineCell = p.ClinePerLength*p.d;
p.LlineCell = p.LlinePerLength*p.d;
p.sampleLength = p.nCells*p.d;
p.thetaScale = p.theta100/p.fReference;   % theta(f)=thetaScale*f
p.bScale = p.b100/p.fReference;           % ideal b(f)=bScale*f
p.Rload = p.switchRon + p.Rpar;
p.ZonDesign = p.Z0/2;                     % SI design target near 50 MHz
p.frequencyRatioDesign = p.ZonDesign/p.Z0;
end

function p = localApplyStruct(p,s)
if ~isscalar(s)
    error('otr:parameters:StructScalar','Override structure must be scalar.');
end
names = fieldnames(s);
derived = {'phaseVelocity','phaseVelocityMethods','ClinePerLength', ...
    'LlinePerLength','ClineCell','LlineCell','sampleLength','thetaScale', ...
    'bScale','Rload','ZonDesign','frequencyRatioDesign'};
for k = 1:numel(names)
    % Passing a previously returned parameter structure is intentionally
    % idempotent: ignore its derived fields and recompute them below.
    if ~any(strcmpi(names{k},derived))
        p = localSet(p,names{k},s.(names{k}));
    end
end
end

function p = localSet(p,name,value)
names = fieldnames(p);
idx = find(strcmpi(name,names),1);
if isempty(idx)
    error('otr:parameters:UnknownParameter', ...
        'Unknown OTR parameter ''%s''.',name);
end
p.(names{idx}) = value;
end

function localValidate(p)
positive = {'c0','d','Z0','epsEff','Cload','Lpar','nCells', ...
    'Cnearest','CnextNearest','theta100','b100','switchingTime', ...
    'fReference','phaseVelocity','phaseVelocityMethods', ...
    'ClinePerLength','LlinePerLength','ClineCell','LlineCell'};
nonnegative = {'switchRon','Rpar','alpha100','lossTangent'};
for k = 1:numel(positive)
    v = p.(positive{k});
    if ~(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v > 0)
        error('otr:parameters:BadValue','%s must be a positive finite scalar.',positive{k});
    end
end
for k = 1:numel(nonnegative)
    v = p.(nonnegative{k});
    if ~(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v >= 0)
        error('otr:parameters:BadValue','%s must be a nonnegative finite scalar.',nonnegative{k});
    end
end
if p.nCells ~= round(p.nCells)
    error('otr:parameters:CellCount','nCells must be a positive integer.');
end
end
