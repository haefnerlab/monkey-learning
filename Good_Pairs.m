function pairs = Good_Pairs(pop, params)
%GOOD_PAIRS finds indexes (i,j,k,..) where
% 1) they form a unique set (i.e. 1,2 and 2,1 not both included) and
% 2) all neurons i, j, and k have well-defined tuning.
%
%A neuron is 'well tuned' if it has significant variation in the fixation
%condition as measured by a p-value from 'anova' or 'fprime_pvalue'. Which
%one is used is controlled by params.exclusion_rule where neurons with a
%value above params.exclusion_threshold is thrown out.
%
% "pairs" is a bit of a misnomer; with moment=3 it's triples, etc.

if params.diagonal, k = 0;
else k = 1; end

% part 1: get unique pairs
n_neurons = length(pop.cellnos);
unq_pairs = find(Util.ndtriu(n_neurons*ones(1,params.moment), k));

% part 2: keep only neurons that are 'well tuned' as defined either by ANOVA test
% or by t-test of fprime being different from zero.
well_tuned = (pop.(params.exclusion_rule) < params.exclusion_threshold);

% now find where "all cells i,j,k,... had good tuning" by multiplying the
% logical array well_tuned by itself using a generalization of outer 
% product to n-dimensions (and there are 'moment' dimensions here)
well_tuned_idxs = find(Util.ndouter(well_tuned, params.moment));

% end result is the intersection of the two sets of indexes
pairs = intersect(unq_pairs, well_tuned_idxs);

end