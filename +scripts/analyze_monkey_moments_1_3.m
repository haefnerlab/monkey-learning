monkeys = {'both'};
moments = 1:3;

for monk_idx = 1:length(monkeys)
    for moment = moments
        memo_file = fullfile('data', monkeys{monk_idx}, sprintf('analyze_fprime_moment%d.mat', moment));
        scripts.analyze_task_offset(...
            New_Parameters('monkey', monkeys{monk_idx}, 'bootstrap', 500, 'moment', moment), ...
            memo_file);
    end
end
