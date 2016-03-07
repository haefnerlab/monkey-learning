function analyze_scatter_moments( params, figpath )
%ANALYZE_SCATTER_MOMENTS compares statistical moments of f' tuning curves
% and 'choice-triggered' distributions (when there is no stimulus)
%
% Comparisons are done for 1st- and 2nd-order predictions
%
% params.params.min_pairs: min # n-tuples of not-NaN trials for the moment to be considered valid 
% params.params.min_rates: min avg # spikes in the n-tuple for a trial to be counted

%% Load and preprocess
[pops_task, pops_fix, full_pops_task, full_pops_fix] = Load_Preprocess(params);

colors = hsv(length(pops_task));

%% compare first moment (diff in means)

figure();
subplot(1,2,1);

% concatenation of all fprimes and CT?Ms for getting correlations
n_all_fprimes = sum(cellfun(@length, {pops_task.fprime_stimulus_mean}));
all_fprimes = zeros(1,n_all_fprimes);
all_ctdms   = zeros(1,n_all_fprimes);
all_selections = zeros(1,n_all_fprimes);
i = 1;

hold on;
for pi=1:length(pops_task)
    pop = pops_task(pi);
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
    % FILTER: tuned to orientation and >minimum spike/sec
    selections = (pop.(params.exclusion_rule) < params.exclusion_threshold) & nanmean(pop.spikeRates_stim0,2)' > params.min_rates;
    scatter(fprime(selections), choice_triggered_delta_means(selections), 5, colors(pi,:));
    
    next_i = i+length(fprime);
    all_fprimes(i:next_i-1) = fprime;
    all_ctdms(i:next_i-1) = choice_triggered_delta_means;
    all_selections(i:next_i-1) = selections;
    i = next_i;
end
hold off;

valid_entries = ~isnan(all_ctdms) & all_selections;
[R, P] = corrcoef(all_fprimes(valid_entries), all_ctdms(valid_entries));

axis square;
title(sprintf('%s :: 1st-order prediction :: Correlation = %.4f :: p=%.2e', params.monkey, R(2), P(2)));
xlabel('f'' / sigma_i');
ylabel('?_choice r_i / sigma_i');
if nargin > 1
    savefig(fullfile(figpath, 'scatter1.fig'));
end

%% compare second moment (noise correlations vs f'f')

subplot(1,2,2);

% concatenation of all fprimes and CT?Ms for getting correlations
corrs = [];
fpfp  = [];
selections = [];

hold on;
for pi=1:length(pops_task)
    pop = pops_task(pi);
    if params.verbose, fprintf('\tPopulation %d of %d (%d neurons)\n', pi, length(pops_task), length(pop.cellnos)); end;
    
    % compare zero-signal correlations to fp_i fp_j / (sigma_i sigma_j)
    variances = nanvar(pop.spikeRates_stim0,1,2);
    sigma_ij = sqrt(variances * variances');
    pop_corrs = Util.nancomoment(pop.spikeRates_stim0', 2, true, true, true, params.min_pairs, params.min_rates);
    pop_fpfp = Util.ndouter(pop.fprime_stimulus_mean', 2) ./ sigma_ij;
    pop_select = logical(Util.ndouter(pop.(params.exclusion_rule) < (params.exclusion_threshold), 2));
    % ignore diagonal
    for i=1:length(variances), pop_select(i,i) = false; end
    % flatten to 1d array of pairwise stats
    pop_corrs = pop_corrs(:);
    pop_fpfp = pop_fpfp(:);
    pop_select = pop_select(:);

    scatter(pop_fpfp(pop_select), pop_corrs(pop_select), 5, colors(pi,:));

    corrs = vertcat(corrs, pop_corrs);
    fpfp = vertcat(fpfp, pop_fpfp);
    selections = vertcat(selections, pop_select);
end
hold off;

valid_entries = ~isnan(corrs) & selections;
[R, P] = corrcoef(fpfp(valid_entries), corrs(valid_entries));

axis square;
title(sprintf('%s :: 2nd-order prediction :: Correlation = %.4f :: p=%.2e', params.monkey, R(2), P(2)));
xlabel('fp_i fp_j / sig_i sig_j');
ylabel('noise correlations');
if nargin > 1
    savefig(fullfile(figpath, 'scatter2.fig'));
end
    
end