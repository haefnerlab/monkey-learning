% generate predictions from sampling simulation
deltas = [0 0.005 0.02 0.08];

% load sampling data
pops_task = Convert_Sampling_Output(e(4,:), 20, 10, size(e{1}.X, 1), 5);

% do our usual preprocessing (note anything involving pops_fix is
% impossible)
pops_task = Split_Conditions(pops_task);
pops_task = Compute_fPrime_stimulus_means(pops_task);

analyze_scatter_moments(...
    New_Parameters('monkey', 'SIM', 'moment', 2, 'min_pairs', 1), ...
    pops_task);