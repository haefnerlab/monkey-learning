monkeys = {'lem', 'jbe', 'both'};
moments = 1:3;

for monk_idx = 1:length(monkeys)
    for moment = moments
        memo_file = fullfile('data', monkeys{monk_idx}, sprintf('analyze_fprime_moment%d.mat', moment));
        analyze_fprime_tuning_curves(...
            New_Parameters('monkey', monkeys{monk_idx}, 'bootstrap', 500), ...
            memo_file);
    end
end
