function [ populations ] = Compute_fPrime_stimulus_means( populations, recompute )
%Compute_fPrime compute f' (change in mean response with change in
%stimulus) and return modified populations with new .fprime_stimulus_mean vector
%
%Also computes sensitivity; the p-value that each f' value is different
%from zero is stored in .fprime_pvalue

if nargin < 2, recompute = false; end

if ~isfield(populations, 'fprime_stimulus_mean') || ~isfield(populations, 'fprime_pvalue') || recompute
    populations = arrayfun(@Compute_fPrime_Single_Pop, populations);
end

end

function [pop] = Compute_fPrime_Single_Pop( pop )

pop.fprime_stimulus_mean = zeros(1,length(pop.cellnos));
pop.fprime_pvalue = zeros(1,length(pop.cellnos));

for n_idx=1:length(pop.cellnos)
    % todo - Poisson variability?
    valid = ~isnan(pop.spikeRates(n_idx, :));
    stats = regstats(pop.spikeRates(n_idx,valid)', pop.condVec(valid), ...
        'linear', 'tstat');
    % linear fit with r(s) = stats.beta(2)*s + stats.beta(1)
    pop.fprime_stimulus_mean(n_idx) = stats.tstat.beta(2);
    pop.fprime_pvalue = stats.tstat.pval(2);
end

end

