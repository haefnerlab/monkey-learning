% top-level script for monkey data analysis

%% Load and preprocess
if ~exist('monkey', 'var')
    monkey = 'lem';
end
if ~exist('populations', 'var')
    fprintf('loading data... ');
    populations = Load_Task_Data(monkey);
    fprintf('done\n');
end
populations = Split_Conditions( populations );
populations = Compute_fPrime_stimulus_means( populations );

verbose = true;
nmoments = 2;
min_pairs = 25;   % min # n-tuples of not-NaN trials for the moment to be considered valid
min_rates = 10; % min avg # spikes in the n-tuple for a trial to be counted

%% compare first moment (diff in means)

figure();
colors = hsv(length(populations));
subplotsquare(nmoments, 1);

% concatenation of all fprimes and CT?Ms for getting correlations
n_all_fprimes = sum(cellfun(@length, {populations.fprime_stimulus_mean}));
all_fprimes = zeros(1,n_all_fprimes);
all_ctdms   = zeros(1,n_all_fprimes);
i = 1;

hold on;
for pi=1:length(populations)
    pop = populations(pi);
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
for moment = 2:nmoments
    subplotsquare(nmoments, moment);

    % concatenation of all fprimes and CT?Ms for getting correlations
    n_all_fprimes = sum(cellfun(@(fp) sum(ndtriu(length(fp) * ones(1,moment))), {populations.fprime_stimulus_mean}));
    all_fprimes = zeros(1,n_all_fprimes);
    all_ctdms   = zeros(1,n_all_fprimes);
    i = 1;

    if verbose, fprintf('Calculating moment %d\n', moment); end

    hold on;
    for pi=1:length(populations)
        pop = populations(pi);
        if verbose, fprintf('\tPopulation %d of %d (%d neurons)\n', pi, length(populations), length(pop.cellnos)); end;
        % get f'f'f'... up to moment times
        [stimulus_moments, ~, indices] = nancomoment(pop.fprime_stimulus_mean, moment, true);
        choice_triggered_delta_means = nancomoment(pop.spikeCounts_stim0', moment, true, min_pairs, min_rates);
        
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
    
clearvars -except populations verbose;