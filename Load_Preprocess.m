function [pops_task, pops_fix, full_pops_task, full_pops_fix] = Load_Preprocess(params)
%LOAD_PREPROCESS load data for params.monkey and do all usual preprocessing

savefile = fullfile('data', params.monkey, 'preprocessed.mat');

% fitting tuning curves is time-consuming; precomputed results are put in
% a 'preprocessed.mat' file
if ~exist(savefile, 'file') || params.recompute_tuning
    fprintf('loading data (%s)... ', params.monkey);
    full_pops_task = Load_Task_Data(params.monkey);
    full_pops_fix = Load_Fixation_Data(params.monkey);
    [pops_task, pops_fix] = Match_Corresponding_Populations( full_pops_task, full_pops_fix );
    fprintf('done\n');
else
    fprintf('loading preprocessed data (%s)... ', params.monkey);
    savedata = load(savefile);
    full_pops_task = savedata.full_pops_task;
    full_pops_fix = savedata.full_pops_fix;
    pops_task = savedata.pops_task;
    pops_fix = savedata.pops_fix;
    fprintf('done\n');
end

[pops_task, pops_fix] = Split_Conditions( pops_task, pops_fix );
pops_task = Compute_fPrime_stimulus_means( pops_task, params.recompute_tuning );
pops_task = Compute_fPrime_bestfit( pops_task, pops_fix, params.recompute_tuning );
pops_task = Compute_fPrime_fixation_means( pops_task, pops_fix, params.recompute_tuning);

% For the 'full' data, we can only get fprime from stimulus means
% (the Match_Corresponding_Populations() function filtered them down into
% only the neurons shared across fixation and task data)
[full_pops_task, full_pops_fix] = Split_Conditions( full_pops_task, full_pops_fix );
full_pops_task = Compute_fPrime_stimulus_means( full_pops_task, params.recompute_tuning );

save(savefile, 'pops_task', 'pops_fix', 'full_pops_task', 'full_pops_fix');

end