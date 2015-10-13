% with 'otune' in the workspace..

deltas = [0, 0.005, 0.02, 0.08];
for di=1:length(deltas)
    for v=1:3
        savefile = fullfile('data', sprintf('SIM_d%.3f_v%d', deltas(di), v), 'preprocessed.mat');
        load(savefile);
        pops_fix = Convert_Sampling_Tuning(otune{di}, pops_task, 30);
        
        % do further preprocessing
        [pops_task, pops_fix] = Match_Corresponding_Populations(pops_task, pops_fix); % should keep everything...
        

        pops_task = Split_Conditions( pops_task );
        pops_task = Compute_fPrime_stimulus_means( pops_task );
        pops_task = Compute_fPrime_bestfit( pops_task, pops_fix );
        pops_task = Compute_fPrime_fixation_means( pops_task, pops_fix );
        
        save(savefile, 'pops_task', 'pops_fix');
    end
end