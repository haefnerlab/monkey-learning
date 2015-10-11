function [ M, L, U ] = meanci( data, interval )
%MEANCI compute mean and lower/upper confidence intervals, generally
%handling NaN values gracefully
%
% [mean,lower,upper] = MEANCI(data, [interval]) where data has measurements
% on rows and variables across columns, returns a row vector for each of
% mean, lower, and upper. interval defaults to 0.95

if nargin < 2, interval = 0.95; end

if ~(interval > 0 && interval < 1)
    error('confidence interval must be in (0,1)');
end

M = nanmean(data,1);
dsorted = sort(data,1); % note that sort() puts NaNs at the end
indexes = 1:size(data,1);

for var=size(data,2):-1:1
    % consider only non-nan values for this var in computing the confidence
    % interval (1 thru n, where n is #non-nan)
    n = sum(~isnan(data(:,var)));
    
    if n == 0
        L(var) = NaN;
        U(var) = NaN;
    else
        lo_idx = n * (1-interval)/2;
        hi_idx = n - lo_idx;

        % get data values at lo_idx and hi_idx, which in general won't be
        % integers, we interpolate in the data for them.
        L(var) = interp1(indexes, dsorted(:,var), lo_idx, 'pchip');
        U(var) = interp1(indexes, dsorted(:,var), hi_idx, 'pchip');
    end
end

end

