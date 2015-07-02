function [ populations ] = Compute_fPrime( populations )
%Compute_fPrime compute f' (change in mean response with change in
%stimulus) and return modified populations with new .fprime vector

if(~isfield(populations, 'spikeCounts_choiceA'))
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
meanA = nanmean(pop.spikeCounts_choiceA, 2);
meanB = nanmean(pop.spikeCounts_choiceB, 2);

pop.fprime = (meanA - meanB)';

end

