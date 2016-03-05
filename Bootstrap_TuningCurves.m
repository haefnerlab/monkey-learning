function [boot_pops_task, boot_pops_fix] = Bootstrap_TuningCurves( pops_task, pops_fix )
%BOOTSTRAP_TUNINGCURVES do bootstrap resampling of pops_fix spike counts at
%each stimulus orientation and recompute tuning curves (results of which
%are saved into boot_pops_task as usual)

boot_pops_task = pops_task;
boot_pops_fix = pops_fix;

%%% resample pops_fix spike counts
for p_idx = 1:length(pops_fix)
    pop = boot_pops_fix(p_idx);
    
    orientation_condition = strcmpi(pop.condVecLabel, 'orientation');
    orientations = pop.condVec(:, orientation_condition);
    unq_orients = unique(pop.condVec(:,orientation_condition));
    
    for o_idx = 1:length(unq_orients)
        o = unq_orients(o_idx);
        indices_this_stim = find(orientations == o);
        n_trials = length(indices_this_stim);
        rand_trials_with_replacement = randi(n_trials, 1, n_trials);
        pop.spikeCounts(:, indices_this_stim) = pop.spikeCounts(:, indices_this_stim(rand_trials_with_replacement));
    end
end
    
%%% recompute maximum likelihood tuning curves
[boot_pops_task, boot_pops_fix] = Split_Conditions(boot_pops_task, boot_pops_fix);
boot_pops_task = Compute_fPrime_bestfit( boot_pops_task, boot_pops_fix, true );
boot_pops_task = Compute_fPrime_fixation_means( boot_pops_task, boot_pops_fix, true);

end