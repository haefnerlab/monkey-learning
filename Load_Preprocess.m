function [pops_task, pops_fix] = Load_Preprocess(params)
%LOAD_PREPROCESS load data for params.monkey and do all usual preprocessing

savefile = fullfile('data', params.monkey, 'preprocessed.mat');

% fitting tuning curves is time-consuming; precomputed results are put in
% a 'preprocessed.mat' file
if ~exist(savefile, 'file') || params.recompute_tuning
    fprintf('loading data... ');
    pops_task = Load_Task_Data(params.monkey);
    pops_fix = Load_Fixation_Data(params.monkey);
    [pops_task, pops_fix] = Match_Corresponding_Populations( pops_task, pops_fix );
    fprintf('done\n');
else
    fprintf('loading preprocessed data... ');
    savedata = load(savefile);
    pops_task = savedata.pops_task;
    pops_fix = savedata.pops_fix;
    fprintf('done\n');
end

pops_task = Split_Conditions( pops_task );
pops_task = Compute_fPrime_stimulus_means( pops_task, params.recompute_tuning );
pops_task = Compute_fPrime_bestfit( pops_task, pops_fix, params.recompute_tuning );
pops_task = Compute_fPrime_fixation_means( pops_task, pops_fix, params.recompute_tuning);

save(savefile, 'pops_task', 'pops_fix');

end