function [ optim, curve, quality ] = fitVonMises( orientations, spikeCounts )
%FITVONMISES ML fit a von Mises tuning curve (see vonMises()) to
% (orientation,count) data. Assumed poisson variability
%
% inputs:
% - orientations: vector of orientations (i.e. stimuli)
% - spikeCounts: vector of spike counts of same length as orientations
%
% return values:
% - params: struct with r_0, r_max, k, theta_pref that define the curve
% - curve: a function handle where curve(orientation) gives our predicted
%          value
% - quality: log-likelihood of data with MLE fit of Poisson distributed
%            counts

N = numel(spikeCounts);

% treat both as column vectors
orientations = reshape(orientations(1:N), N, 1);
spikeCounts  = reshape(spikeCounts,  N, 1);

% remove data points with NaN values
valid_values = ~isnan(spikeCounts) & ~isinf(spikeCounts) & ~isinf(orientations);
orientations = orientations(valid_values);
spikeCounts = spikeCounts(valid_values);

fn_to_minimize = @(params) -sum(log_likelihood(orientations, spikeCounts, ...
    params(1), params(2), params(3), params(4)));

% TODO - intelligent initialization
optim = fminsearch(fn_to_minimize, [0,1,1,pi/2], ...
    optimset('MaxFunEvals', 4000, 'MaxIter', 2000, 'FunValCheck','on'));

curve = @(o) vonMises(o, optim(1), optim(2), optim(3), optim(4));

% we were minimizing negative log likelihood. max likelihood (and the
% corresponding value) is negative of the minimized fn at the solution
quality = -fn_to_minimize(optim);

end

function ll = log_likelihood(o, c, r_0, r_max, k, theta_pref)
% poisson-likelihood of seeing (orientation o, count c) pair under given
% parameterization of vonMises

dist_mean = vonMises(o, r_0, r_max, k, theta_pref);

% Poisson = lambda^k e^-lambda / k!
% so, log Poisson = k log(lambda) - lambda - log(k!)
% where k is our count, and lambda is the mean value as expected from VM
ll = c .* log(dist_mean) - dist_mean - log(factorial(c));

if isnan(ll) | isinf(ll)
    fprintf('wtf?\n');
end

end