function params = test_sampling_prediction(delta, n, params)
% generate predictions from sampling simulation

% load preprocessed sampling data
filename = fullfile('data', 'SIM', sprintf('sample_d%.3f_n%d.mat', delta, n));
if ~exist(filename, 'file')
    error('no sampling model prediction d=%f n=%d', delta, n);
end
data = load(filename);
pops_task = data.pops_task;

if nargin < 3, params = New_Parameters(); end

params.monkey = 'SIM';
% nothing NaN in sampling model; keep it all
params.min_pairs = 1;
params.min_value = 7;

% do our usual preprocessing (note anything involving pops_fix is
% impossible)
pops_task = Split_Conditions(pops_task);
pops_task = Compute_fPrime_stimulus_means(pops_task);

analyze_scatter_moments(params, pops_task);

test_plot_noisecorrelations;
title(sprintf('Predicted noise correlations, delta=%.2f', params.sampling_delta));

end