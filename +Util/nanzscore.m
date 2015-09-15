function [ Z ] = nanzscore( X )
%NANZSCORE nan-safe version of zscore

expand = ones(size(X,1), 1);

% subtract from each column its (nan)mean
Z = X - (expand * nanmean(X));
% divide each column by its (nan)variance
Z = Z ./ (expand * sqrt(nanvar(X,1)));

end

