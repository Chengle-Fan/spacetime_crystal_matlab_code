function out = otr_synthetic_wavepacket(options)
%OTR_SYNTHETIC_WAVEPACKET Analytic wave-packet scattering at time interfaces.
%
% This is a spectral solution of the ideal, spatially homogeneous circuit
% limit used in SI S5.  It complements (but does not replace) the discrete
% Kirchhoff/MNA circuit solver.  Each spatial Fourier component conserves k.
%
% Options (all optional):
%   direction       'on' (OFF->ON), 'off' (ON->OFF), or 'slab'
%   Zoff/Zon        line impedances (default 50/25 ohm)
%   fInitial        carrier before the first interface (default 50 MHz)
%   x, t            SI-unit grids
%   x0, sigma       initial Gaussian centre and amplitude width
%   tSwitch         first time interface
%   tau             temporal-slab duration
%   initialVelocity phase velocity before the first interface

if nargin < 1, options = struct(); end
Zoff = get_option(options,'Zoff',50);
Zon = get_option(options,'Zon',25);
direction = lower(get_option(options,'direction','on'));
fInitial = get_option(options,'fInitial',50e6);
vOff = get_option(options,'initialVelocity',1.04e8);
switch direction
    case {'on','slab'}
        vInitial = vOff;
        vOther = vOff*Zon/Zoff;
        firstLaw = 'charge';
        Zbefore = Zoff;
        Zafter = Zon;
    case 'off'
        vInitial = vOff*Zon/Zoff;
        vOther = vOff;
        firstLaw = 'voltage';
        Zbefore = Zon;
        Zafter = Zoff;
    otherwise
        error('direction must be ''on'', ''off'', or ''slab''.');
end

lambda = vInitial/fInitial;
x = get_option(options,'x',linspace(0,8*lambda,1201));
tSwitch = get_option(options,'tSwitch',40e-9);
t = get_option(options,'t',linspace(0,140e-9,701));
x0 = get_option(options,'x0',2.2*lambda);
sigma = get_option(options,'sigma',0.65*lambda);
tau = get_option(options,'tau',25e-9);

x = x(:).'; t = t(:);
Nx = numel(x); dx = mean(diff(x));
kGrid = 2*pi*[0:floor((Nx-1)/2), -ceil((Nx-1)/2):-1]/(Nx*dx);
k0 = 2*pi*fInitial/vInitial;
u0 = exp(-((x-x0)/sigma).^2).*exp(1i*k0*(x-x0));
A0 = fft(u0);

% Avoid the unphysical |k| cusp only at the numerically negligible DC bin.
wInitial = vInitial*abs(kGrid);
wOther = vOther*abs(kGrid);
[~,Tfirst,Rfirst] = otr_temporal_interface_local( ...
    Zbefore,Zafter,firstLaw,direction);

field = complex(zeros(numel(t),Nx));
forward = field; backward = field;
for it = 1:numel(t)
    tq = t(it);
    if tq < tSwitch
        Af = A0.*exp(-1i*wInitial*tq);
        Ab = zeros(size(Af));
    elseif strcmp(direction,'slab') && tq >= tSwitch+tau
        [~,Ton,Ron] = otr_temporal_interface_local(Zoff,Zon,'charge','on');
        [~,Toff,Roff] = otr_temporal_interface_local(Zon,Zoff,'voltage','off');
        pPlus = exp(-1i*wOther*tau);
        pMinus = exp(1i*wOther*tau);
        Af = A0.*exp(-1i*wInitial*tSwitch).*( ...
            Toff*Ton.*pPlus + Roff*Ron.*pMinus).* ...
            exp(-1i*wInitial*(tq-tSwitch-tau));
        Ab = A0.*exp(-1i*wInitial*tSwitch).*( ...
            Roff*Ton.*pPlus + Toff*Ron.*pMinus).* ...
            exp(1i*wInitial*(tq-tSwitch-tau));
    else
        elapsed = tq-tSwitch;
        Af = Tfirst*A0.*exp(-1i*wInitial*tSwitch).*exp(-1i*wOther*elapsed);
        Ab = Rfirst*A0.*exp(-1i*wInitial*tSwitch).*exp(1i*wOther*elapsed);
    end
    forward(it,:) = ifft(Af);
    backward(it,:) = ifft(Ab);
    field(it,:) = forward(it,:)+backward(it,:);
end

out.x = x; out.t = t; out.field = field;
out.forward = forward; out.backward = backward; out.initial = u0;
out.k = kGrid; out.k0 = k0; out.omegaInitial = wInitial;
out.omegaOther = wOther; out.fRatio = vOther/vInitial;
out.Tfirst = Tfirst; out.Rfirst = Rfirst;
out.tSwitch = tSwitch; out.tau = tau; out.direction = direction;
end

function [M,T,R] = otr_temporal_interface_local(Zbefore,Zafter,law,direction)
% Local fallback makes this helper usable while keeping the public API thin.
if exist('otr_temporal_interface','file') == 2
    [M,T,R] = otr_temporal_interface(Zbefore,Zafter,law);
    return;
end
switch lower(law)
    case 'charge'
        ratio = (Zafter/Zbefore)^2;
        rootRatio = Zafter/Zbefore;
        T = 0.5*(ratio+rootRatio); R = 0.5*(ratio-rootRatio);
    case 'voltage'
        rootRatio = Zafter/Zbefore;
        T = 0.5*(1+rootRatio); R = 0.5*(1-rootRatio);
    otherwise
        error('Unknown law for %s transition.',direction);
end
M = [T R;R T];
end

function value = get_option(options,name,defaultValue)
if isfield(options,name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end
