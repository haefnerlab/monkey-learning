% top-level script for monkey data analysis

%% Load and preprocess
% monkey = 'jbe';
monkey = 'lem';

savefile = fullfile('data', monkey, 'preprocessed.mat');

if ~exist('pops_task', 'var') || ~exist('pops_fix', 'var')
    if ~exist(savefile, 'file')
        fprintf('loading data... ');
        pops_task = Load_Task_Data(monkey);
        pops_fix = Load_Fixation_Data(monkey);
        [pops_task, pops_fix] = Match_Corresponding_Populations( pops_task, pops_fix );
        fprintf('done\n');
    else
        savedata = load(savefile);
        pops_task = savedata.pops_task;
        pops_fix = savedata.pops_fix;
    end
end
pops_task = Split_Conditions( pops_task );
pops_task = Compute_fPrime_stimulus_means( pops_task );
pops_task = Compute_fPrime_bestfit( pops_task, pops_fix );

save(savefile, 'pops_task', 'pops_fix');

n_pops = length(pops_task);
n_neurons = sum(arrayfun(@(p) length(p.cellnos), pops_task));

verbose = true;
min_pairs = 25; % min # n-tuples of not-NaN trials for the moment to be considered valid
min_rates = 10; % min avg # spikes in the n-tuple for a trial to be counted
popcolors = hsv(n_pops);

close all;

%% First moment: correlation of f' with change in choice-triggered-means.
%
% 1x3 subplots:
%  Left: scatter and correlation when f' aligned with task
%  Middle: scatter and correlation when f' 45 degrees off task
%  Right: smoothly varying the two; correlation as a function of distance
%         off task

TuningCurve_fPrime_At = @(curve, angle) (curve(angle) - curve(angle+90));

f1 = figure();
% LEFT subplot
all_fprimes = zeros(n_neurons,1);
all_ctdms = zeros(n_neurons,1);
start_idx = 1;

subplot(1,3,1);
hold on;
for p_idx=1:n_pops
    pop = pops_task(p_idx);
    popsize = length(pop.cellnos);
    choice_triggered_delta_means = (nanmean(pop.spikeCounts_choiceA,2)-nanmean(pop.spikeCounts_choiceB,2))';
%     fprime_at_task = arrayfun(@(n_idx) TuningCurve_fPrime_At(pop.tuning_vm_curves{n_idx}, pop.Orientation), 1:popsize);
    fprime_at_task = pop.fprime_stimulus_mean;
    scatter(fprime_at_task, choice_triggered_delta_means, 5, popcolors(p_idx,:));
    
    end_idx = start_idx + popsize - 1;
    all_fprimes(start_idx:end_idx) = fprime_at_task;
    all_ctdms(start_idx:end_idx) = choice_triggered_delta_means;
    start_idx = end_idx + 1;
end
hold off;
[R,P] = corrcoef(all_fprimes, all_ctdms);
title(sprintf('Corr. choice-triggered means vs task f''\nr=%.3f p=%.3e', R(2), P(2)));
xlabel('f'' aligned to task')
ylabel('choice-triggered diff means')

% MIDDLE subplot
all_fprimes = zeros(n_neurons,1);
all_ctdms = zeros(n_neurons,1);
start_idx = 1;

subplot(1,3,2);
hold on;
for p_idx=1:n_pops
    pop = pops_task(p_idx);
    popsize = length(pop.cellnos);
    choice_triggered_delta_means = (nanmean(pop.spikeCounts_choiceA,2)-nanmean(pop.spikeCounts_choiceB,2))';
    fprime_off_task = arrayfun(@(n_idx) TuningCurve_fPrime_At(pop.tuning_vm_curves{n_idx}, pop.Orientation + 45), 1:popsize);
    scatter(fprime_off_task, choice_triggered_delta_means, 5, popcolors(p_idx,:));
    
    end_idx = start_idx + popsize - 1;
    all_fprimes(start_idx:end_idx) = fprime_off_task;
    all_ctdms(start_idx:end_idx) = choice_triggered_delta_means;
    start_idx = end_idx + 1;
end
hold off;
[R,P] = corrcoef(all_fprimes, all_ctdms);
title(sprintf('Corr. choice-triggered means vs off-task f''\nr=%.3f p=%.3e', R(2), P(2)));
xlabel('f'' aligned 45 degrees off task')
ylabel('choice-triggered diff means')

% RIGHT subplot

offsets = -90:5:90;

all_correlations = zeros(popsize, length(offsets));
for o_idx=1:length(offsets)
    offset = offsets(o_idx);
    for p_idx=1:n_pops
        pop = pops_task(p_idx);
        popsize = length(pop.cellnos);
        choice_triggered_delta_means = (nanmean(pop.spikeCounts_choiceA,2)-nanmean(pop.spikeCounts_choiceB,2))';
        fprime_at_offset = arrayfun(@(n_idx) TuningCurve_fPrime_At(pop.tuning_vm_curves{n_idx}, pop.Orientation + offset), 1:popsize);

        R = corrcoef(fprime_at_offset, choice_triggered_delta_means);
        all_correlations(p_idx, o_idx) = R(2);
    end
end

means = mean(all_correlations, 1);
stds = std(all_correlations, 1);

subplot(1,3,3);
errorbar(offsets, means, stds);
title(sprintf('Corr. f'' of tuning curve with choice-triggered means\nas a function of distance-from-trial-center'));
xlabel('offset from trial alignment');
ylabel('correlation of tuning curve f'' at \theta+offset with choice-triggered means');