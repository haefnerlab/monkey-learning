function analyze_scatter_moments( params, figpath )
%ANALYZE_SCATTER_MOMENTS compares statistical moments of f' tuning curves
% and 'choice-triggered' distributions (when there is no stimulus)
%
% Comparisons are done for 1st- and 2nd-order predictions
%
% params.params.min_pairs: min # n-tuples of not-NaN trials for the moment to be considered valid 
% params.params.min_rates: min avg spike rate (all trials) for a neuron to be included

%% Load and preprocess
[pops_task, pops_fix, full_pops_task, full_pops_fix] = Load_Preprocess(params);

colors = hsv(length(pops_task));

%% compare first moment (diff in means)

figure();
params.moment = 1;

% concatenation of all fprimes and CTdMs for getting correlations
n_pops = length(pops_task);
all_fprimes = cell(1,n_pops);
all_ctdms   = cell(1,n_pops);

hold on;
for p_idx=1:n_pops
    pop = pops_task(p_idx);
    % bugfix: what is called 'A' or 'B' is not consistently the positive
    % direction of the stimulus axis 's'. This causes a large fraction of
    % the choice_triggered_delta_means to have the wrong sign, which
    % completely ruins the correlation.. We correct for it here.
    signal_trials = pop.condVec > 0;
    % dot product between 'correctChoice' and 'sign(condVec)' will be
    % positive if they align and negative if they don't.
    choice_sign = sign(pop.correctChoice(signal_trials)' * sign(pop.condVec(signal_trials)));
    choice_triggered_delta_means = (nanmean(pop.spikeRates_choiceA,2)-nanmean(pop.spikeRates_choiceB,2))';
    choice_triggered_delta_means = choice_sign * choice_triggered_delta_means;
    % normalize by standard deviation
    variances = nanvar(pop.spikeRates_stim0,1,2)';
    choice_triggered_delta_means = choice_triggered_delta_means ./ sqrt(variances);
    fprime = pop.fprime_stimulus_mean ./ sqrt(variances);
    % FILTER for well-tuned neurons and minimum rate
    selections = Good_Pairs(pop, params);
    fprime = fprime(selections);
    choice_triggered_delta_means = choice_triggered_delta_means(selections);
    % scatter plot this population
    scatter(fprime, choice_triggered_delta_means, 5, colors(p_idx,:));
    
    all_fprimes{p_idx} = fprime;
    all_ctdms{p_idx} = choice_triggered_delta_means;
end
hold off;

% concatenate together across all populations to get total correlation
all_fprimes = horzcat(all_fprimes{:});
all_ctdms = horzcat(all_ctdms{:});

valid_entries = ~isnan(all_ctdms);
disp(sum(valid_entries) / length(valid_entries));
[R, P] = corrcoef(all_fprimes(valid_entries), all_ctdms(valid_entries));

axis square;
title(sprintf('%s :: 1st-order prediction :: Correlation = %.4f :: p=%.2e', params.monkey, R(2), P(2)));
xlabel('f'' / sigma_i');
ylabel('?_choice r_i / sigma_i');
if nargin > 1
    savefig(fullfile(figpath, 'scatter1.fig'));
end

%% compare second moment (noise correlations vs f'f')

figure();
params.moment = 2;

% concatenation of all fprimes and CTdMs for getting correlations
fpfp  = cell(1,n_pops);
corrs = cell(1,n_pops);

hold on;
for p_idx=1:n_pops
    pop = pops_task(p_idx);
    if params.verbose, fprintf('\tPopulation %d of %d (%d neurons)\n', p_idx, n_pops, length(pop.cellnos)); end;
    
    % compare zero-signal correlations to fp_i fp_j / (sigma_i sigma_j)
    pop_corrs = Util.nancomoment(pop.spikeRates_stim0', 2, true, true, true, params.min_pairs);
    
    variances = nanvar(pop.spikeRates_stim0,1,2);
    sigma_ij = sqrt(variances * variances');
    pop_fpfp = Util.ndouter(pop.fprime_stimulus_mean', 2) ./ sigma_ij;
    selections = Good_Pairs(pop, params);
    % flatten to 1d array of pairwise stats
    pop_corrs = pop_corrs(selections);
    pop_fpfp = pop_fpfp(selections);

    scatter(pop_fpfp, pop_corrs, 5, colors(p_idx,:));

    corrs{p_idx} = pop_corrs;
    fpfp{p_idx}  = pop_fpfp;
end
hold off;

% concatenation of all fprimes and CTdMs for getting correlations
fpfp  = horzcat(fpfp{:});
corrs = horzcat(corrs{:});

valid_entries = ~isnan(corrs);
[R, P] = corrcoef(fpfp(valid_entries), corrs(valid_entries));

axis square;
title(sprintf('%s :: 2nd-order prediction :: Correlation = %.4f :: p=%.2e', params.monkey, R(2), P(2)));
xlabel('fp_i fp_j / sig_i sig_j');
ylabel('noise correlations');
if nargin > 1
    savefig(fullfile(figpath, 'scatter2.fig'));
end
    
end