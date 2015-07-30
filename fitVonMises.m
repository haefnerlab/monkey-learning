function [ best_params, curve, best_map, worst_params, worst_map ] = fitVonMises( orientations, spikeCounts, use_map, n_init )
%FITVONMISES ML or MAP fit a von Mises tuning curve (see vonMises()) to
% (orientation,count) data. Assumed poisson variability and some
% mostly-arbitrary priors on parameters
%
% inputs:
% - orientations: vector of orientations (i.e. stimuli)
% - spikeCounts: vector of spike counts of same length as orientations
%
% return values:
% -  best_params: struct with r_0, r_max, k, theta_pref that define the curve
% - curve: a function handle where curve(orientation) gives our predicted
%          value
% - quality: log-posterior of data with MAP fit of Poisson distributed
%            counts

if nargin < 3
    use_map = true;
end

if nargin < 4
    n_init = 20;
end

N = numel(spikeCounts);

% treat both as column vectors
orientations = reshape(orientations(1:N), N, 1);
spikeCounts  = reshape(spikeCounts,  N, 1);

% remove data points with NaN values
valid_values = ~isnan(spikeCounts) & ~isinf(spikeCounts) & ~isinf(orientations);
orientations = orientations(valid_values);
spikeCounts = spikeCounts(valid_values);

% Priors for sampling model parameters
count_mean = mean(spikeCounts);
count_min = min(spikeCounts);
count_max = max(spikeCounts);
orientation_spacing = min(diff(sort(unique(orientations))));

r_0_mean = count_mean;
r_max_mean = 0;
r_max_dev = sqrt((count_max - count_min));
k_mean = 2*orientation_spacing;
k_dev = 10;

init_r_0   = @() exprnd(r_0_mean);
init_r_max = @() normrnd(r_max_mean, r_max_dev);
init_k     = @() normrnd(k_mean, k_dev);
init_th    = @() rand * 180;

log_p_r_0   = @(r_0)   -explike(r_0_mean, r_0);
log_p_r_max = @(r_max) -normlike([r_max_mean, r_max_dev], r_max);
log_p_k     = @(k)     -normlike([k_mean, k_dev], k);
log_p_th    = @(th)    log(1.0/180);

if use_map
    fn_to_maximize = @(params) sum(log_likelihood(orientations, spikeCounts, params)) ...
        + log_p_r_0(params(1)) ...
        + log_p_r_max(params(2)) ...
        + log_p_k(params(3)) ...
        + log_p_th(params(4));
else
    % Max Likelihood version is MAP sans priors, and otherwise everything
    % is identical (still sample initial values n_init times, etc)
    fn_to_maximize = @(params) sum(log_likelihood(orientations, spikeCounts, params));
end

fn_to_minimize = @(params) -fn_to_maximize(params);

% posterior is non-convex, so we sample initialization n_init times and
% keep track of which had the best MAP score
% (and for the sake of teaching and debugging, also the worst one)
best_map = -inf;
worst_map = inf;
for run=1:n_init

    params = [init_r_0(), init_r_max(), init_k(), init_th()];
    optim = fminsearch(fn_to_minimize, params, ...
        optimset('MaxFunEvals', 4000, 'MaxIter', 2000));

    score = fn_to_maximize(optim);
    if ~isreal(score)
        fprintf('wtf?\n');
    end
    if score > best_map
        best_map = score;
        best_params = optim;
        fprintf('best updated to %f\n', best_map);
    end
    
    if score < worst_map
        worst_map = score;
        worst_params = optim;
        fprintf('worst updated to %f\n', worst_map);
    end
end

% we were minimizing negative log likelihood. max likelihood (and the
% corresponding value) is negative of the minimized fn at the solution
curve = @(o) vonMises(o, best_params);

end

function ll = log_likelihood(o, c, params)
% poisson-likelihood of seeing (orientation o, count c) pair under given
% parameterization of vonMises

dist_mean = vonMises(o, params);

% Poisson = lambda^k e^-lambda / k!
% so, log Poisson = k log(lambda) - lambda - log(k!)
% where k is our count, and lambda is the mean value as expected from VM
ll = c .* log(dist_mean) - dist_mean - log(factorial(c));

end