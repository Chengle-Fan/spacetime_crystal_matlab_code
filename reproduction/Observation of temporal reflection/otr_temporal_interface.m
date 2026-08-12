function [J,T,R,info] = otr_temporal_interface(Zbefore,Zafter,law,varargin)
%OTR_TEMPORAL_INTERFACE Temporal scattering with the SI microscopic laws.
%   [J,T,R,INFO] = OTR_TEMPORAL_INTERFACE(ZBEFORE,ZAFTER,LAW) returns the
%   amplitude map [V_after+; V_after-] = J*[V_before+; V_before-].
%
%   LAW='charge'  : connecting capacitance, C2*v(0+)=C1*v(0-) (S14),
%                   and coefficients in Eq. S16 ('on').
%   LAW='voltage' : removing capacitance, v(0+)=v(0-) (S15), and
%                   coefficients in Eq. S17 ('off').
%   LAW='auto' chooses charge when |Zafter|<|Zbefore| and voltage otherwise.
%   Aliases 'on'/'S16' and 'off'/'S17' are accepted. T=J(1,1), R=J(2,1).
%
%   For dispersive calculations ZBEFORE/ZAFTER may be equal-size arrays;
%   J is 2-by-2-by-N. The SI identities assume a nonmagnetic TL so that
%   C_before/C_after=(Z_after/Z_before)^2.
%
%   Optional 'CapacitanceRatio',Cbefore/Cafter checks/overrides that ratio.

if nargin < 3 || isempty(law), law = 'auto'; end
if ~(isnumeric(Zbefore) && isnumeric(Zafter) && ~isempty(Zbefore) && ...
        all(isfinite(Zbefore(:))) && all(isfinite(Zafter(:))) && ...
        all(real(Zbefore(:))>0) && all(real(Zafter(:))>0))
    error('otr:interface:Impedance','Impedances must have positive finite real parts.');
end
if ~(isscalar(Zbefore) || isscalar(Zafter) || isequal(size(Zbefore),size(Zafter)))
    error('otr:interface:Size','Impedances must be scalar or equal-size arrays.');
end
if isscalar(Zbefore), Zbefore = Zbefore+zeros(size(Zafter)); end
if isscalar(Zafter), Zafter = Zafter+zeros(size(Zbefore)); end

capRatio = (Zafter./Zbefore).^2;
if mod(numel(varargin),2) ~= 0
    error('otr:interface:Options','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = lower(char(varargin{k})); value = varargin{k+1};
    switch name
        case 'capacitanceratio'
            if ~(isnumeric(value) && all(isfinite(value(:))) && all(real(value(:))>0))
                error('otr:interface:Capacitance','CapacitanceRatio must be positive and finite.');
            end
            capRatio = value+zeros(size(Zbefore));
        otherwise
            error('otr:interface:Options','Unknown option ''%s''.',name);
    end
end

law = localLaw(law,Zbefore,Zafter);
if strcmp(law,'charge')
    % S16: Ton=1/2(C1/C2+sqrt(C1/C2)), Ron=difference.
    g = sqrt(capRatio);
    T = 0.5*(capRatio+g);
    R = 0.5*(capRatio-g);
    equation = 'S14/S16'; direction = 'on';
else
    % S17: Toff=1/2(1+sqrt(C1/C2)), Roff=difference.
    g = sqrt(capRatio);
    T = 0.5*(1+g);
    R = 0.5*(1-g);
    equation = 'S15/S17'; direction = 'off';
end

% Reciprocity of the real field gives the second column by swapping the
% positive/negative-frequency incident amplitudes: symmetric temporal J.
n = numel(T);
J = complex(zeros(2,2,n));
tv = T(:); rv = R(:);
for k = 1:n, J(:,:,k) = [tv(k),rv(k);rv(k),tv(k)]; end
if isscalar(T), J = J(:,:,1); end

info = struct('law',law,'direction',direction,'equation',equation, ...
    'Zbefore',Zbefore,'Zafter',Zafter,'capacitanceRatio',capRatio, ...
    'frequencyRatio',sqrt(capRatio), ...
    'boundaryScale',T+R,'difference',T-R);
end

function law = localLaw(law,zb,za)
if ~(ischar(law) || (isstring(law) && isscalar(law)))
    error('otr:interface:Law','law must be text.');
end
law = lower(strtrim(char(law)));
switch law
    case {'charge','on','s14','s16','connect'}
        law = 'charge';
    case {'voltage','off','s15','s17','disconnect'}
        law = 'voltage';
    case 'auto'
        down = abs(za) < abs(zb);
        if all(down(:)), law = 'charge';
        elseif all(~down(:)), law = 'voltage';
        else
            error('otr:interface:AutoArray', ...
                'AUTO cannot mix increasing and decreasing impedance in one call.');
        end
    otherwise
        error('otr:interface:Law','Unknown microscopic law ''%s''.',law);
end
end
