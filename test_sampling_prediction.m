function params = test_sampling_prediction(e, params)

% generate predictions from sampling simulation
deltas = [0 0.005 0.02 0.08]; % what was actually used in sampling model

if nargin < 2
    params = New_Parameters();
end

params.monkey = 'SIM';
params.min_pairs = 1;
params.min_value = 7;

if ~any(params.sampling_delta == deltas)
    error('sampling model was not done on delta=%f', params.sampling_delta);
end

% load sampling data
pops_task = Convert_Sampling_Output(e(deltas==params.sampling_delta,:), params.sampling_npops, params.sampling_nneurons, size(e{1}.X, 1), 5);

% do our usual preprocessing (note anything involving pops_fix is
% impossible)
pops_task = Split_Conditions(pops_task);
pops_task = Compute_fPrime_stimulus_means(pops_task);

analyze_scatter_moments(...
    params, ...
    pops_task);

test_plot_noisecorrelations;
title(sprintf('Predicted noise correlations, delta=%.2f', params.sampling_delta));

end