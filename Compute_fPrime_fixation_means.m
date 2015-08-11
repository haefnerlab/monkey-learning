function [ pops_task ] = Compute_fPrime_fixation_means(pops_task, pops_fix)
%Compute_fPrime_fixation_means compute f' as from fixatation data using mean
% values from that task (or simple linear interpolation between them)

for p_idx=1:length(pops_task)
    pt = pops_task(p_idx);
    pf = pops_fix(p_idx);
    task_orientation = pt.Orientation;
    task_ortho = 1 + mod(task_orientation + 89, 180);
    
    orient_idx = strcmp('orientation', pf.condVecLabel);
    fixation_orientations = pf.condVec(:, orient_idx);
    
    usable_trials = ~isinf(fixation_orientations) & ~isnan(fixation_orientations);
    fixation_orientations = fixation_orientations(usable_trials);
    counts = pf.spikeCounts(:,usable_trials);
    
    unique_orientations = unique(fixation_orientations);
    
    for n_idx=length(pt.cellnos):-1:1
        spikeCount_means = arrayfun(@(o) nanmean(counts(n_idx, fixation_orientations == o)), unique_orientations);
        interpolated_spikecount = interp1(unique_orientations, spikeCount_means, task_orientation);
        interpolated_spikecount_offset = interp1(unique_orientations, spikeCount_means, task_ortho);
        pops_task(p_idx).fprime_fixation_means(n_idx) = interpolated_spikecount - interpolated_spikecount_offset;
    end
end

end