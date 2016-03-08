function [ boot_pops_task ] = Bootstrap_SpikeCounts(pops_task)
%BOOTSTRAP_SPIKECOUNTS given our standard pops_task array of recordings,
%returns a nearly identical array where the spikeCounts and spikeRates
%trials have been resampled with replacement.

boot_pops_task = pops_task;
for p_idx=1:length(pops_task)
    pop = pops_task(p_idx);
    
    % separately sample spike counts for each stimulus condition + choice
    % combination
    conditions = unique(pop.condVec);
    for s_idx=1:length(conditions)
        s = conditions(s_idx);
        for choice=[-1,1]
            trials_this_condition = find(pop.condVec == s & pop.realChoice == choice);
            n_trials = length(trials_this_condition);
            % resample trials at this condition (with replacement)
            resampled_trials = randi(n_trials, 1, n_trials);
            pop.spikeCounts(:, trials_this_condition) = ...
                pop.spikeCounts(:, trials_this_condition(resampled_trials));
        end
    end
end

boot_pops_task = Split_Conditions(boot_pops_task);
boot_pops_task = Compute_fPrime_stimulus_means( boot_pops_task );

end