function [ output ] = nanmoment( X, order, dim )
%NANMOMENT same interface as the moment() function, but like nanmean(), it
%just ignores NaN values

if nargin < 3, dim = 1; end

if ~any(isnan(X))
    output = moment(X, order, dim);
else
    output = nanfndim(@(vec) moment(vec, order), X, dim);
end

end

