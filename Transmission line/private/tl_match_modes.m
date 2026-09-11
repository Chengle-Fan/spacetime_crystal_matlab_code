function order = tl_match_modes(overlap)
%TL_MATCH_MODES Exact assignment for the two-to-four-state TL bulk models.
n = size(overlap,1);
candidate = perms(1:n);
indices = (1:n)+n*(candidate-1);
[~,best] = max(sum(overlap(indices),2));
order = candidate(best,:);
end
