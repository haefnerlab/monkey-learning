function [ vv ] = ndouter( v, n )
%NDOUTER take 'n-dimensional' outer product of a vector v with itself. If v
%has length m, output is m x m x m x ... x m (n times total)

grids = cell(1,n);

[grids{:}] = ndgrid(v);

vv = grids{1};

for i=2:n
    vv = vv .* grids{i};
end

end

