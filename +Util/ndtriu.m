function [ indices ] = ndtriu( size, k )
%NDTRIU Get 'upper triangular' indices for an nd-array of the given size

if nargin < 2, k=0; end

ndims = length(size);

if ndims == 1
    indices = 1:size(1);
else
    indices = arrayfun(@(i) is_unique_index(size,i,k), 1:prod(size));
end

end

function [ unq ] = is_unique_index( size, i, k )
% return true iff the given nd-index should be used. 
%
% we keep any nondecreasing subscripts
sub = cell(1,length(size));
[sub{:}] = ind2sub(size, i);
idx_diff = diff(cell2mat(sub));
unq = all(idx_diff >= 0);
unq = unq && min(idx_diff) >= k;

end