function analyze_noisecorrelations_bothmonkeys(params)

if nargin < 1, params = New_Parameters(); end

%% Get noisecorrelation images for each monkey (180px X 180px)

image_total = zeros(180);
counts_total = zeros(180);

monkeys = {'lem', 'jbe'};
stim_conditions = {'spikeRates_stim0', 'spikeRates_stim0A', 'spikeRates_stim0B'};

figure();

for monk_idx = 1:length(monkeys);
    monkey = monkeys{monk_idx};
    [pops_task, pops_fix] = Load_Preprocess(params);

    % get a plot of noise correlations by stimulus condition and monkey
    for stim_idx = 1:length(stim_conditions)
        subplot(length(monkeys), length(stim_conditions), (monk_idx-1)*length(stim_conditions)+stim_idx);
        [stim_mean,~] = Vis.plot_noisecorrelations(pops_task, params, false, stim_conditions(stim_idx));
        Util.imagescnan(stim_mean); colorbar;
        title(sprintf('%s %s', monkey, stim_conditions{stim_idx}));
    end
    % get a total plot of averages
    [monkey_mean, monkey_count] = Vis.plot_noisecorrelations(pops_task, params, false, stim_conditions);
    countable = monkey_count ~= 0;
    image_total(countable) = image_total(countable) + monkey_mean(countable) .* monkey_count(countable);
    counts_total = counts_total + monkey_count;
end

%% plot totals

image_total(counts_total == 0) = NaN;
image_total = image_total ./ counts_total;

figure();
Util.imagescnan(image_total); colorbar;
title('All noise correlations, both monkeys included (task=0,90)');

end