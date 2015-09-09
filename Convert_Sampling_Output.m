function [ pops_task ] = Convert_Sampling_Output(e, n_pops, n_neurons, n_trials, n_contrasts_task)
%CONVERT_SAMPLING_OUTPUT convert from the struct output from the sampling
%model to a populations struct that looks like pops_task and pops_fix as
%used elsewhere in this repository
%
% [ pops_task ] = CONVERT_SAMPLING_OUTPUT(e, ...) where e is the
% output of the sampling model, makes structs pops_task and pops_fix that
% can be analyzed like the monkey data, but are missing some obvious things
% like date of recording, type of probe, etc...
%
% CONVERT_SAMPLING_OUTPUT(e, n_pops, n_neurons, n_trials, n_contrasts_task)
% creates n_pops populations, each with n_neurons drawn randomly from the X
% simulated neurons. Each population will contain n_trials trials that draw
% stimuli from +/- n_contrasts_task/2 contrasts.
%
% We can't get a struct like pops_fix easily since we would need multiple
% orientations, then we would need to use the same populations to do the
% 2AFC task.

n_contrasts_simulated = length(e);
assert(n_contrasts_simulated >= n_contrasts_task);
assert(mod(n_contrasts_task, 2) == 1, 'need an odd number of contrasts so zero will be included');

[n_trials_sim, n_neurons_sim, n_samples] = size(e{1}.X);
assert(n_neurons <= n_neurons_sim);
assert(n_trials <= n_trials_sim);

% scale down the timescale assuming 0.2 seconds per sample and a total
% trial length of 2 seconds
sim_seconds = n_samples * 0.2;
scale_by = 2 / sim_seconds;

h_contrasts = cellfun(@(expt) expt.InputImage.stimulus_contrast(1), e);
v_contrasts = cellfun(@(expt) expt.InputImage.stimulus_contrast(2), e);

% let's say contrasts -12:3:12 were tested, but we only want 3 of them
% (n_contrasts_task = 3). That means we want [-3,0,3]. This is gotten by
% 'trimming' [-12,-9,-6] and [6,9,12] from each side.
contrasts_tested = unique([-h_contrasts v_contrasts]);
n_trim = (length(contrasts_tested) - n_contrasts_task)/2;
contrasts_sample = contrasts_tested(n_trim+1:end-n_trim);

for p_idx = n_pops:-1:1
    random_neurons = randperm(n_neurons_sim, n_neurons);
    contrast_each_trial = randsample(contrasts_sample, n_trials, true);
    n_zerostim = sum(contrast_each_trial == 0);
    trial_within_contrast = randperm(n_trials_sim, n_trials);
    correct_choices = sign(contrast_each_trial);
    correct_choices(contrast_each_trial == 0) = round(rand(n_zerostim,1))*2-1;

    model_choices = arrayfun(...
        @(t) sampler_choice(...
            e{contrast_each_trial(t) == contrasts_tested}, ... % get experiment where this trial's contrast was used
            trial_within_contrast(t)), ... % ...and use this trial within that experiment
        1:n_trials);

    model_spikecounts = zeros(n_neurons, n_trials);
    for ti = 1:n_trials
        e_this_trial = e{contrast_each_trial(ti) == contrasts_tested};
        for ni = 1:n_neurons
            model_spikecounts(ni,ti) = get_spikecount(e_this_trial, random_neurons(ni), trial_within_contrast(ti), scale_by);
        end
    end
    
    pops_task(p_idx) = struct(...
        'Header', struct('monkey', 'SIM', 'description', 'sampling model recast as monkey data'), ...
        'condVec', contrast_each_trial, ...
        'condVecLabel', {'Signal Strength'}, ...
        'Orientation', 0, ... % TODO explicitly break into multiple tasks
        'trialDur', 2000, ...
        'correctChoice', correct_choices, ...
        'realChoice', model_choices, ...
        'reward', [], ...
        'trialSeed', [], ...
        'spikeTimes', [], ...
        'spikeCounts', model_spikecounts, ...
        'cellnos', random_neurons,...
        'tuning', e{1}.Projection.phi_x(random_neurons) * 180 / pi);
end

end

function choice = sampler_choice(e, trial)
    sampled_orientations = e.O(trial, 1, e.Sampling.number_burn_in+1:end);

    % assume perfect evidence accumulation
    if sum(sampled_orientations == 1) > sum(sampled_orientations == 2)
        choice = -1;
    else
        choice = 1;
    end
end

function count = get_spikecount(e, neuron, trial, scale)
    count = scale * sum(squeeze(e.X(trial,neuron,:)));
end