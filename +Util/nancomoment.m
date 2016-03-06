function [ output, counts, indices, all_indices ] = nancomoment( X, order, symmetries, sub_mean, norm_variance, minimum_count, minimum_value )
%NANCOMOMENT get the order-th co-moment of X, ignoring NaN values, where each
% column of X is a variable and each row an observation
%
% output is an ndarray with 'order' dimensions, each of size #variables.
%
% [output, counts, indices] = nancomoment( X, order, [symmetries, [norm_variance, [minimum_count, [minimum_value]]]] )
%
% output[i,j,k,...] is the expected (i.e. mean) value of
%   (X(:,i)-u(i))*(X(:,j)-u(j))*(X(:,k)-u(k))*...
% where 'u' is the mean of each variable, and we discard any (i,j,k,...)
% pairs where one or more of the measurements is NaN, and the mean of the
% measurements is above minimum_value. If norm_variance is true, the output
% is divided by sqrt(var(i)*var(j)*var(k)*...) (i.e. in moment=2 case, this
% gives correlation rather than covariance)
%
% counts[i,j,k,...] is the number of non-nan observation involving i,j, and
% k. Where counts is less than minimum_count, output is NaN
%
% all_indices is a vector of indices into the output where a value was
%   computed. If 'symmetries' is true, this includes at most one of each
%   unique combination of variables (IE find(triu(output)) in the 2d case)
%
% indices is all_indices but with any 'invalid' outputs removed. By 'valid'
%   we mean that it is not NaN, and the minimum_count and minimum_value
%   requirements are satisfied

if order <= 0, error('order of moments must be greater than 0'); end

if nargin < 3, symmetries=false; end
if nargin < 4, sub_mean=true; end
if nargin < 5, norm_variance=false; end
if nargin < 6, minimum_count=1; end
if nargin < 7, minimum_value=-inf; end

stddevs = sqrt(nanvar(X,1));
means = nanmean(X,1);

if order == 1
    output = nanmean(X, 1); % built-in function is faster for means
else
    osize = size(X,2) * ones(1, order);
    nddots = zeros(osize);
    counts = zeros(osize);
    
    all_indices = 1:numel(nddots);
    if symmetries, all_indices = find(Util.ndtriu(size(nddots))); end
    
    % ind2sub uses varargout, which we will capture in this cell array
    nd_idxs_cell = cell(1,order);
    
    for i=all_indices
        [nd_idxs_cell{:}] = ind2sub(osize, i);
        idxs = cell2mat(nd_idxs_cell);
        
        vecs = X(:,idxs);
        
        % find which indices are not NaN for all vecs across observations
        all_valid_observations = all(~isnan(vecs), 2);
        
        % also only keep those where the mean (across vecs) is >
        % minimum_value
        % TODO - allow geometric mean in params
        all_valid_observations = all_valid_observations & mean(vecs,2) > minimum_value;
        
        % count how many data points survived the not-nan and min-value
        % filters
        counts(i) = sum(all_valid_observations);
        vecs = vecs(all_valid_observations,:);
        
        % subtract mean to get moment around the center
        if sub_mean, vecs = vecs - ones(counts(i),1) * means(idxs); end
        
        % compute order-th moment
        % like a dot product with <order>-many vectors instead of just 2
        nddots(i) = sum(prod(vecs, 2), 1);
        
        if norm_variance
            % normalize by sqrt(var i * var j * var k ...)
            nddots(i) = nddots(i) / prod(stddevs(idxs));
        end
    end
    
    indices = all_indices;
    indices(counts(indices) < minimum_count) = [];
    nddots(counts < minimum_count) = NaN;
    output = nddots ./ counts;
end

end

