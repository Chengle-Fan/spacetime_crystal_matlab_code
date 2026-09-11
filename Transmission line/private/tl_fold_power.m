function [folded,orders] = tl_fold_power(raw,rawOrders,nFold)
%TL_FOLD_POWER Sum physical DFT bins modulo Omega along the first dimension.
% Integer-bin mapping works for centered/unshifted and odd/even input grids.
% Every remaining dimension (probe, channel or k) is retained independently.
shape = size(raw);
mapping = sparse(mod(rawOrders(:),nFold)+1,(1:shape(1)).', ...
    1,nFold,shape(1));
orders = (-floor(nFold/2):ceil(nFold/2)-1).';
folded = mapping(mod(orders,nFold)+1,:)*reshape(raw,shape(1),[]);
folded = reshape(folded,[nFold,shape(2:end)]);
end
