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

% get means, ignore NaN values
%
% note that this is the same as using our helper function 
%   nanfndim(@mean, counts, 2)
% but is built-in and slightly faster
meanA = nanmean(pop.spikeCounts_stimA, 2);
meanB = nanmean(pop.spikeCounts_stimB, 2);

pop.fprime_stimulus_mean = (meanA - meanB)';

end

