function analyze_scatter_moments( params )
%ANALYZE_SCATTER_MOMENTS compares statistical moments of f' tuning curves
% and 'choice-triggered' distributions (when there is no stimulus)
%
% Comparisons are done up to the params.moment moment
%
% params.params.min_pairs: min # n-tuples of not-NaN trials for the moment to be considered valid 
% params.params.min_rates: min avg # spikes in the n-tuple for a trial to be counted

%% Load and preprocess
[pops_task, pops_fix] = Load_Preprocess(params);

colors = hsv(length(pops_task));

%% compare first moment (diff in means)

figure();
Util.subplotsquare(params.moment, 1);

% concatenation of all fprimes and CT?Ms for getting correlations
n_all_fprimes = sum(cellfun(@length, {pops_task.fprime_stimulus_mean}));
all_fprimes = zeros(1,n_all_fprimes);
all_ctdms   = zeros(1,n_all_fprimes);
i = 1;

hold on;
for pi=1:length(pops_task)
    pop = pops_task(pi);
    choice_triggered_delta_means = (nanmean(pop.spikeCounts_choiceA,2)-nanmean(pop.spikeCounts_choiceB,2))';
    scatter(pop.fprime_stimulus_mean, choice_triggered_delta_means, 5, colors(pi,:));
    
    next_i = i+length(pop.fprime_stimulus_mean);
    all_fprimes(i:next_i-1) = pop.fprime_stimulus_mean;
    all_ctdms(i:next_i-1) = choice_triggered_delta_means;
    i = next_i;
end
hold off;

valid_correlations = ~isnan(all_ctdms);
[R, P] = corrcoef(all_fprimes(valid_correlations), all_ctdms(valid_correlations));

axis square;
title(sprintf('Moment %d :: Correlation = %.4f :: p=%.2e', 1, R(2), P(2)));
xlabel('tuning curve f'' statistics');
ylabel('zero-stimulus noise statistics');

%% Compare higher-order moments
for moment = 2:params.moment
    Util.subplotsquare(params.moment, moment);

    % concatenation of all fprimes and CT?Ms for getting correlations
    n_all_fprimes = sum(cellfun(@(fp) sum(Util.ndtriu(length(fp) * ones(1,moment))), {pops_task.fprime_stimulus_mean}));
    all_fprimes = zeros(1,n_all_fprimes);
    all_ctdms   = zeros(1,n_all_fprimes);
    i = 1;

    if params.verbose, fprintf('Calculating moment %d\n', moment); end

    hold on;
    for pi=1:length(pops_task)
        pop = pops_task(pi);
        if params.verbose, fprintf('\tPopulation %d of %d (%d neurons)\n', pi, length(pops_task), length(pop.cellnos)); end;
        % get f'f'f'... up to moment times
        [stimulus_moments, ~, indices] = Util.nancomoment(pop.fprime_stimulus_mean, moment, true, false);
        choice_triggered_delta_means = Util.nancomoment(pop.spikeCounts_stim0', moment, true, true, true, params.min_pairs, params.min_rates);
        
        scatter(stimulus_moments(indices), choice_triggered_delta_means(indices), 5, colors(pi,:));
        
        next_i = i + length(indices);
        all_fprimes(i:next_i-1) = stimulus_moments(indices);
        all_ctdms(i:next_i-1)   = choice_triggered_delta_means(indices);
        i = next_i;
    end
    hold off;

    valid_correlations = ~isnan(all_ctdms);
    [R, P] = corrcoef(all_fprimes(valid_correlations), all_ctdms(valid_correlations));

    axis square;
    title(sprintf('Moment %d :: Correlation = %.4f :: p=%.2e', moment, R(2), P(2)));
    xlabel('tuning curve f'' statistics');
    ylabel('zero-stimulus noise statistics');
end
    
end