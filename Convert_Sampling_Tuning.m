function [ pops_fix ] = Convert_Sampling_Tuning( tuning, sister_pops_task, keep_trials )
%CONVERT_SAMPLING_TUNING the sister function to Convert_Sampling_Output.
%This function takes as input a pops_task that came from the sampling model
%as well as the orientation tuning struct array results from the sampling
%model, and returns an analogous pops_fix struct array that has each
%orientation tested keep_trials times
%(keep_trials must be no more than the number of trials in tuning)

n_pops = length(sister_pops_task);
[~, n_samples, n_trials_tested] = size(tuning(1).X);
orientations = [tuning.orientation];
n_orientations = length(orientations);
trial_milliseconds = 425;

assert(keep_trials <= n_trials_tested);

% scale down the timescale assuming 0.2 seconds per sample and a total
% trial length of trial_dur milliseconds
sim_seconds = n_samples * 0.2;
trial_seconds = trial_milliseconds / 1000.0;
scale_by = trial_seconds / sim_seconds;

for p_idx = n_pops:-1:1
    
    cellnos = sister_pops_task(p_idx).cellnos;
    
    rand_trial_indices = randperm(n_orientations*n_trials_tested, n_orientations*keep_trials);
    orientations_by_trial = zeros(n_orientations*keep_trials,1);
    
    model_spikecounts = zeros(length(cellnos), keep_trials);
    for i=1:length(rand_trial_indices)
        t = rand_trial_indices(i);
        o_idx = floor((t-1) / n_trials_tested) + 1;
        t_idx = mod((t-1), n_trials_tested) + 1;
        
        orientations_by_trial(i) = orientations(o_idx);
        model_spikecounts(:,i) = round(scale_by * sum(tuning(o_idx).X(cellnos, :, t_idx),2));
    end
    
    pops_fix(p_idx) = struct(...
        'Header', struct('monkey', 'SIM', 'description', 'sampling model orientation tuning recast as pops_fix', 'SessionNumber', sister_pops_task(p_idx).Header.SessionNumber),...
        'condVec', orientations_by_trial, ...
        'condVecLabel', {'orientation'}, ...
        'cellnos', cellnos, ...
        'trialDur', trial_milliseconds, ...
        'spikeCounts', model_spikecounts,...
        'spikeTimes', {cell(length(cellnos),1)});
end

end

