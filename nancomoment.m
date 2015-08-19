function [ output, counts, indices ] = nancomoment( X, order, symmetries, minimum_count, minimum_value )
%NANCOMOMENT get the order-th co-moment of X, ignoring NaN values, where each
% column of X is a variable and each row an observation
%
% output is an ndarray with 'order' dimensions, each of size #variables.
%
% [output, counts, indices] = nancomoment( X, order, [symmetries, [minimum_count, [minimum_value]]] )
%
% output[i,j,k,...] is the expected (i.e. mean) value of
%   (X(:,i)-u(i))*(X(:,j)-u(j))*(X(:,k)-u(k))*...
% where 'u' is the mean of each variable, and we discard any (i,j,k,...)
% pairs where one or more of the measurements is NaN, and the mean of the
% measurements is above minimum_value
%
% counts[i,j,k,...] is the number of non-nan pairs found. Where counts is 
% zero, output is NaN

if order <= 0, error('order of moments must be greater than 0'); end

if nargin < 3, symmetries=false; end
if nargin < 4, minimum_count=1; end
if nargin < 5, minimum_value=-inf; end

if order == 1
    output = nanmean(X, 1); % built-in function is faster for means
else
    osize = size(X,2) * ones(1, order);
    nddots = zeros(osize);
    counts = zeros(osize);
    
    indices = 1:numel(nddots);
    if symmetries, indices = find(ndtriu(size(nddots))); end
    
    % subtract mean from X
    Xzero = X - repmat(nanmean(X), size(X,1), 1);
    
    % ind2sub uses varargout, which we will capture in this cell array
    nd_idxs_cell = cell(1,order);
    
    for i=indices
        [nd_idxs_cell{:}] = ind2sub(osize, i);
        idxs = cell2mat(nd_idxs_cell);
        
        vecs = Xzero(:,idxs);
        
        % find which indices are not NaN for all vecs across observations
        all_valid_observations = all(~isnan(vecs), 2);
        
        % also only keep those where the mean (across vecs) is >
        % minimum_value
        all_valid_observations = all_valid_observations & mean(X(:,idxs),2) > minimum_value;
        
        % count how many data points survived the not-nan and min-value
        % filters
        counts(i) = sum(all_valid_observations);
        
        % compute order-th moment
        % like a dot product with <order>-many vectors instead of just 2
        nddots(i) = sum(prod(vecs(all_valid_observations, :), 2), 1);
    end
    
    indices(counts(indices) < minimum_count) = [];
    nddots(counts < minimum_count) = NaN;
    output = nddots ./ counts;
end

end

