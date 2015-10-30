% with 'e' and 'otune' in the workspace..
% where 'e' is result of Sampling Model/S_Run_Experiments
% and 'otune' is the result of Sampling Model/Orientation_Tuning_Like(e,..)

e_types = [e{:,1}];
e_projections = [e_types.Projection];
deltas = [e_projections.delta];

for di=1:length(deltas)
    for v=1:3
        savefile = fullfile('data', sprintf('SIM_d%.3f_v%d', deltas(di), v), 'preprocessed.mat');
        disp(savefile);
        
        pops_task = Convert_Sampling_Output(e(di,:), 20, 20, 1000, 5);
        pops_fix = Convert_Sampling_Tuning(otune{di}, pops_task, 30);
        
        % do further preprocessing
        [pops_task, pops_fix] = Match_Corresponding_Populations(pops_task, pops_fix); % should keep everything...
        
        disp('fitting');
        pops_task = Split_Conditions( pops_task );
        pops_task = Compute_fPrime_stimulus_means( pops_task );
        pops_task = Compute_fPrime_bestfit( pops_task, pops_fix );
        pops_task = Compute_fPrime_fixation_means( pops_task, pops_fix );
        
        disp('saving');
        save(savefile, 'pops_task', 'pops_fix');
    end
end