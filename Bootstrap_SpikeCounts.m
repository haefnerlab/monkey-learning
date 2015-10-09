function [ boot_pops_task ] = Bootstrap_SpikeCounts(pops_task)
%BOOTSTRAP_SPIKECOUNTS given our standard pops_task array of recordings,
%returns a nearly identical array where the spikeCounts trials have been
%resampled with replacement.

boot_pops_task = pops_task;

for p_idx=1:length(pops_task)
    n_trials = size(pops_task(p_idx).spikeCounts, 2);
    rand_trials_with_replacement = randi(n_trials, 1, n_trials);
    boot_pops_task(p_idx).spikeCounts = pops_task(p_idx).spikeCounts(:,rand_trials_with_replacement);
end

end