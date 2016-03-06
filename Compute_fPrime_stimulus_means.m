function [ populations ] = Compute_fPrime_stimulus_means( populations, recompute )
%Compute_fPrime compute f' (change in mean response with change in
%stimulus) and return modified populations with new .fprime_stimulus_mean vector

if nargin < 2, recompute = false; end

if(~isfield(populations, 'spikeCounts_choiceA')) || recompute
    populations = Split_Conditions( populations );
end

populations = arrayfun(@Compute_fPrime_Single_Pop, populations);

end

function [pop] = Compute_fPrime_Single_Pop( pop )

for n_idx=1:length(pop.cellnos)
    % todo - Poisson variability?
    valid = ~isnan(pop.spikeRates(n_idx, :));
    coeffs = polyfit(pop.condVec(valid), pop.spikeRates(n_idx,valid)', 1);
    % linear fit with r(s) = coeffs(1)*s + coeffs(0)
    pop.fprime_stimulus_mean(n_idx) = coeffs(1);
end

end

