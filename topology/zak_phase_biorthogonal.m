function [zak, info] = zak_phase_biorthogonal( ...
    rightStates, leftStates, metric, sewing)
%ZAK_PHASE_BIORTHOGONAL Gauge-invariant single-band Wilson loop.
%
% rightStates(:,j) and leftStates(:,j) are generalized right/left
% eigenvectors. The link is L_j'*metric*R_{j+1}. The k grid must exclude
% the duplicate right endpoint. sewing maps the first state to k+G.

if nargin < 4 || isempty(sewing)
    sewing = eye(size(rightStates,1));
end
if ~isequal(size(rightStates),size(leftStates))
    error('rightStates and leftStates must have identical sizes.');
end

Nk = size(rightStates,2);
links = complex(zeros(1,Nk));
for ik = 1:Nk-1
    links(ik) = leftStates(:,ik)'*metric*rightStates(:,ik+1);
end
links(Nk) = leftStates(:,Nk)'*metric*sewing*rightStates(:,1);

if any(abs(links) < 1e-12)
    error(['A Wilson link is numerically singular. The band may touch ', ...
        'another band or the Fourier cutoff may be too small.']);
end
unitLinks = links./abs(links);
wilsonLoop = prod(unitLinks);
zak = wrap_to_pi(-angle(wilsonLoop));

info.links = links;
info.unitLinks = unitLinks;
info.wilsonLoop = wilsonLoop;
info.minimumLinkMagnitude = min(abs(links));
end

function value = wrap_to_pi(value)
value = mod(value+pi,2*pi)-pi;
end
