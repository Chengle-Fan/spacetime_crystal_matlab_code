function [chern, curvature, links] = fhs_chern_number( ...
    rightStates, leftStates, metric)
%FHS_CHERN_NUMBER Fukui-Hatsugai-Suzuki Chern number on a periodic grid.
%
% rightStates has size [Nbasis, Noccupied, Nk, Np]. For a Hermitian
% problem, omit leftStates and metric. For a generalized/non-Hermitian
% problem, supply biorthogonal leftStates and the overlap metric so each
% link is det(L_i'*metric*R_j).
%
% Both grid directions are closed periodically. Duplicate endpoints must
% not be included.

if nargin < 2 || isempty(leftStates)
    leftStates = rightStates;
end
if nargin < 3 || isempty(metric)
    metric = eye(size(rightStates,1));
end
if ndims(rightStates) ~= 4 || ~isequal(size(rightStates),size(leftStates))
    error(['rightStates and leftStates must have identical size ', ...
        '[Nbasis,Noccupied,Nk,Np].']);
end
[~,~,Nk,Np] = size(rightStates);
linkK = complex(zeros(Nk,Np));
linkP = complex(zeros(Nk,Np));

for ik = 1:Nk
    ikNext = mod(ik,Nk)+1;
    for ip = 1:Np
        ipNext = mod(ip,Np)+1;
        overlapK = leftStates(:,:,ik,ip)'*metric ...
            *rightStates(:,:,ikNext,ip);
        overlapP = leftStates(:,:,ik,ip)'*metric ...
            *rightStates(:,:,ik,ipNext);
        valueK = det(overlapK);
        valueP = det(overlapP);
        if abs(valueK) < 1e-13 || abs(valueP) < 1e-13
            error('A FHS link is singular at grid point (%d,%d).',ik,ip);
        end
        linkK(ik,ip) = valueK/abs(valueK);
        linkP(ik,ip) = valueP/abs(valueP);
    end
end

curvature = zeros(Nk,Np);
for ik = 1:Nk
    ikNext = mod(ik,Nk)+1;
    for ip = 1:Np
        ipNext = mod(ip,Np)+1;
        plaquette = linkK(ik,ip)*linkP(ikNext,ip) ...
            /(linkK(ik,ipNext)*linkP(ik,ip));
        curvature(ik,ip) = angle(plaquette);
    end
end

chern = sum(curvature(:))/(2*pi);
links.k = linkK;
links.parameter = linkP;
end
