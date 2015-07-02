function [ output, counts ] = nancomoment( X, order )
%NANCOMOMENT get the order-th co-moment of X, ignoring NaN values, where each
% column of X is a variable and each row an observation
%
% output is an ndarray with 'order' dimensions, each of size #variables.
%
% output[i,j,k,...] is the expected (i.e. mean) value of
%   (X(:,i)-u(i))*(X(:,j)-u(j))*(X(:,k)-u(k))*...
% where 'u' is the mean of each variable, and we discard any (i,j,k,...)
% pairs where one or more of the measurements is NaN.
%
% counts[i,j,k,...] is the number of non-nan pairs found. Where counts is 
% zero, output is NaN

if order <= 0, error('order of moments must be greater than 0'); end

if order == 1
    output = nanmean(X, 1); % built-in function is faster for means
else
    osize = size(X,2) * ones(1, order);
    nddots = zeros(osize);
    counts = zeros(osize);
    
    % subtract mean from X
    X = X - repmat(nanmean(X), size(X,1), 1);
    
    % ind2sub uses varargout, which we will capture in this cell array
    nd_idxs_cell = cell(1,order);
    
    % TODO use symmetries; generalize 'upper triangular' to 'upper
    % N-dimensional pyramid' of the ndarray output
    for i=1:numel(nddots)
        [nd_idxs_cell{:}] = ind2sub(osize, i);
        idxs = cell2mat(nd_idxs_cell);
        
        vecs = X(:,idxs);
        
        % find which indices are not NaN for all <order>-many vecs across
        % observations. prod() is used like a logical AND here
        all_not_nan_observations = logical(prod(~isnan(vecs), 2));
        
        counts(i) = sum(all_not_nan_observations);
        
        % like a dot product with <order>-many vectors instead of just 2
        nddots(i) = sum(prod(vecs(all_not_nan_observations, :), 2), 1);
    end
    
    nddots(counts == 0) = NaN;
    output = nddots ./ counts;
end

end

