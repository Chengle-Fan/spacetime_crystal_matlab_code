function value = stpwe_lookup_coefficient(table, m, n)
%STPWE_LOOKUP_COEFFICIENT Read a sampled coefficient; return zero outside.

im = find(table.mOrders == m,1);
in = find(table.nOrders == n,1);
if isempty(im) || isempty(in)
    value = 0;
else
    value = table.values(im,in);
end
end
