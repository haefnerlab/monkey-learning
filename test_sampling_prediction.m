function params = test_sampling_prediction(delta, n, params)
% generate predictions from sampling simulation

% load preprocessed sampling data
filename = fullfile('data', sprintf('SIM_d%.3f_v%d', delta, n), 'preprocessed.mat');
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

scatter_correlation(pops_task, params);

Vis.plot_noisecorrelations(pops_task, params);
title(sprintf('Predicted noise correlations, delta=%.2f', params.sampling_delta));

end

function f = scatter_correlation(pops_task, params)

f = figure();
colors = hsv(length(pops_task));

for p_idx = 1:length(pops_task)
    pop = pops_task(p_idx);
    n_neurons = length(pop.cellnos);
    [noise_correlations,~,indices] = Util.nancomoment(pop.spikeCounts_stim0', 2, false, true, true, params.min_pairs);
    correlations = noise_correlations(indices);
    variances = nanvar(pop.spikeCounts_stim0',1);
    
    fprime_moment = Util.nancomoment(pop.fprime_stimulus_mean, 2, false, false);
    fprime_moment_norm = fprime_moment ./ sqrt(variances * variances');
    choice_triggered_delta_means = (nanmean(pop.spikeCounts_choiceA,2)-nanmean(pop.spikeCounts_choiceB,2))';
    
    subplot(1,3,1);
    hold on;
    scatter(pop.fprime_stimulus_mean.*pop.fprime_stimulus_mean, variances, 8, colors(p_idx,:));
    
    subplot(1,3,2);
    hold on;
    scatter(fprime_moment_norm(indices), correlations, 8, colors(p_idx,:));
%     [r,p] = corr(pop.fprime_moment_norm(indices), choice_triggered_delta_means);
%     title(sprintf('r = %.3e, p = %.3e', r,p));
    
    subplot(1,3,3);
    hold on;
    [sort_tuning, sort_idxs] = sort(pop.tuning);
    plot(sort_tuning, variances(sort_idxs), 'o', 'Color', colors(p_idx, :));

end

subplot(1,3,1);
title('var_i vs f_i^2');
xlabel('f_i^2');
ylabel('var_i');
subplot(1,3,2);
xlabel('f_i''f_j''/sqrt(f_if_j)');
ylabel('corr');

end