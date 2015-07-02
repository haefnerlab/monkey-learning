function [ output ] = nanfndim( fn, X, dim )
%fundim apply fn (which takes in 1D vector and outputs a scalar) along dim 
% of X, ignoring NaN values

if nargin < 3, dim = 1; end

if dim > ndims(X), error('X does not have enough dimensions'); end

% here begins strange code that does fn() on an arbitratry dim
%
% high-level summary: X is indexed using 1D indices and a 'stride' such
% that we pick out the elements along dim.
%
% note that output has same size as X, but with sz(dim) = 1. This means
% we can use sub2ind(Xsize, ind2sub(OUTsize, i)) to do the 'hard' part
% of finding starting indices for each slice along dim
sz = size(X);

stride = prod(sz(1:dim-1)); % spacing between elements on dim
size_along_dim = sz(dim);

osz = sz;
osz(dim) = 1;
output = zeros(osz);

% ind2sub uses varargout, which we will capture in this cell array
idxs_cell = cell(1,ndims(output));

for oi=1:numel(output)
    [idxs_cell{:}] = ind2sub(osz, oi);
    x_idx = sub2ind(sz, idxs_cell{:});
    end_x_idx = x_idx + stride * (size_along_dim - 1);

    X_vec = X(x_idx:stride:end_x_idx);
    output(oi) = fn(X_vec(~isnan(X_vec)));
end

end
