function [M,info] = otr_abcd_unit_cell(f,state,model,p)
%OTR_ABCD_UNIT_CELL ABCD matrix of one switched TLM cell (SI Eq. S1).
%   M = OTR_ABCD_UNIT_CELL(F,STATE,MODEL,P) returns a 2-by-2-by-N array.
%   F is frequency in Hz. STATE is 'off' or 'on'; MODEL is:
%     'ideal'     - SI design shunt, b(f)=2.753*f/(100 MHz);
%     'realistic' - switchRon+Rpar, Lpar and Cload in series.
%   P defaults to OTR_PARAMETERS().  OFF means no shunt load.
%
%   Equation S1 is written in normalized ABCD variables [V; Z0*I].
%   This routine returns the physical matrix [V; I], hence B is multiplied
%   by Z0 and C divided by Z0. INFO.normalizedMatrix is exactly Eq. S1.
%
%   单元由半段无载传输线、并联支路、半段传输线组成。理想模型逐项严格
%   等价于补充材料式 (S1)；realistic 模型保留相同 T 网络而令负载为串联 RLC。

if nargin < 2 || isempty(state), state = 'on'; end
if nargin < 3 || isempty(model), model = 'ideal'; end
if nargin < 4 || isempty(p), p = otr_parameters(); end
localFrequency(f);
state = localChoice(state,{'off','on'},'state');
model = localChoice(model,{'ideal','realistic'},'model');
localParameters(p);

shape = size(f);
fv = f(:).';
omega = 2*pi*fv;
theta = p.thetaScale*fv;

if strcmp(state,'off')
    yNorm = complex(zeros(size(fv)));
elseif strcmp(model,'ideal')
    % The SI ideal model is specified by the independently reported design
    % susceptance b(100 MHz)=2.753.  It is not omega*Cload*Z0: the latter
    % belongs to the practical 82 pF circuit and equals 2.576 at 100 MHz.
    yNorm = 1i*p.bScale*fv;
else
    Zload = p.switchRon + p.Rpar + 1i*omega*p.Lpar + 1./(1i*omega*p.Cload);
    yNorm = p.Z0./Zload;
end

% General symmetric T-network: half line, shunt load, half line. With
% yNorm=j*b this reduces exactly to every entry printed in Eq. S1.
c = cos(theta);
s = sin(theta);
A = c + 0.5i*yNorm.*s;
Bnorm = 1i*s + 0.5*yNorm.*(c-1);
Cnorm = 1i*s + 0.5*yNorm.*(c+1);
D = A;

n = numel(fv);
M = complex(zeros(2,2,n));
Mn = complex(zeros(2,2,n));
for k = 1:n
    Mn(:,:,k) = [A(k),Bnorm(k);Cnorm(k),D(k)];
    M(:,:,k) = [A(k),p.Z0*Bnorm(k);Cnorm(k)/p.Z0,D(k)];
end
if isscalar(f), M = M(:,:,1); Mn = Mn(:,:,1); end

info = struct();
info.frequency = reshape(fv,shape);
info.omega = reshape(omega,shape);
info.theta = reshape(theta,shape);
info.normalizedAdmittance = reshape(yNorm,shape);
info.b = reshape(imag(yNorm),shape);
info.normalizedMatrix = Mn;
info.state = state;
info.model = model;
info.referenceConvention = '[V; I] physical, [V; Z0*I] normalized';
end

function localFrequency(f)
if ~(isnumeric(f) && isreal(f) && ~isempty(f) && all(isfinite(f(:))) && all(f(:) > 0))
    error('otr:abcd:Frequency','f must contain positive finite frequencies in Hz.');
end
end

function out = localChoice(in,allowed,label)
if ~(ischar(in) || (isstring(in) && isscalar(in)))
    error('otr:abcd:Choice','%s must be text.',label);
end
out = lower(strtrim(char(in)));
if ~any(strcmp(out,allowed))
    error('otr:abcd:Choice','Unknown %s ''%s''.',label,out);
end
end

function localParameters(p)
need = {'Z0','Cload','switchRon','Rpar','Lpar','thetaScale','bScale'};
if ~isstruct(p) || ~isscalar(p) || ~all(isfield(p,need))
    error('otr:abcd:Parameters','p must be an OTR_PARAMETERS structure.');
end
end
