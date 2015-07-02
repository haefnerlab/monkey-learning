function [ indices ] = ndtriu( size )
%NDTRIU Get 'upper triangular' indices for an nd-array of the given size

indices = arrayfun(@(i) is_unique_index(size,i), 1:prod(size));

end

function [ unq ] = is_unique_index( size, i )
% return true iff the given nd-index should be used. 
%
% we keep any nondecreasing subscripts
sub = cell(1,length(size));
[sub{:}] = ind2sub(size, i);
unq = all(diff(cell2mat(sub)) >= 0);

end