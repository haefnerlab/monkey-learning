function [ boot_pops_task ] = Bootstrap_SpikeCounts(pops_task)
%BOOTSTRAP_SPIKECOUNTS given our standard pops_task array of recordings,
%returns a nearly identical array where the spikeCounts trials have been
%resampled with replacement.

spike_types = {'spikeCounts_stimA', 'spikeCounts_stimB', ...
    'spikeCounts_stim0', 'spikeCounts_stim0A', 'spikeCounts_stim0B', ...
    'spikeCounts_choiceA', 'spikeCounts_choiceB'};

boot_pops_task = pops_task;
for p_idx=1:length(pops_task)
    % separately sample spike counts for each stimulus condition
    for typ=1:length(spike_types)
        name = spike_types{typ};
        originals = pops_task(p_idx).(spike_types{typ});
        n_trials = size(originals, 2);
        if n_trials > 0
            rand_trials_with_replacement = randi(n_trials, 1, n_trials);
            boot_pops_task(p_idx).(name) = originals(:,rand_trials_with_replacement);
        end
    end

    % note that the above process is an approximation to what real
    % bootstrapping should do since some spikeCounts are subsets of others but
    % we resample each separately, which will likely break the subset
    % relationship. This *does not* affect analyze_fprim_tuningcurves, and
    % likely does not affect others
    
    boot_pops_task(p_idx).spikeCounts = horzcat(...
        boot_pops_task(p_idx).spikeCounts_stimA,...
        boot_pops_task(p_idx).spikeCounts_stimB,...
        boot_pops_task(p_idx).spikeCounts_stim0);
end

boot_pops_task = Compute_fPrime_stimulus_means( boot_pops_task );

end