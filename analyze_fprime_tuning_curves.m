% top-level script for monkey data analysis

%% Load and preprocess
if ~exist('monkey', 'var')
    monkey = 'lem';
end

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
pops_task = Compute_fPrime_fixation_means( pops_task, pops_fix );

save(savefile, 'pops_task', 'pops_fix');

n_pops = length(pops_task);
n_neurons = sum(arrayfun(@(p) length(p.cellnos), pops_task));

verbose = true;
fprime_curve = 'tuning_vm_curves'; % f' defined by von Mises best-fit curve
% fprime_curve = 'tuning_pw_curves'; % f' defined by piecewise linear curve
min_pairs = 25; % (TODO?) min # n-tuples of not-NaN trials for the moment to be considered valid
min_rates = 10; % (TODO?) min avg # spikes in the n-tuple for a trial to be counted
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
if verbose, fprintf('Left plot: f'' task vs CT moment\n'); end
all_fprimes = zeros(n_neurons,1);
all_ctdms = zeros(n_neurons,1);
start_idx = 1;

subplot(1,3,1);
hold on;
for p_idx=1:n_pops
    pop = pops_task(p_idx);
    popsize = length(pop.cellnos);
    choice_triggered_delta_means = (nanmean(pop.spikeCounts_choiceA,2)-nanmean(pop.spikeCounts_choiceB,2))';
    fprime_at_task = arrayfun(@(n_idx) TuningCurve_fPrime_At(pop.(fprime_curve){n_idx}, pop.Orientation), 1:popsize);
    
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
if verbose, fprintf('Middle plot: f'' ortho vs CT moment\n'); end
all_fprimes = zeros(n_neurons,1);
all_ctdms = zeros(n_neurons,1);
start_idx = 1;

subplot(1,3,2);
hold on;
for p_idx=1:n_pops
    pop = pops_task(p_idx);
    popsize = length(pop.cellnos);
    choice_triggered_delta_means = (nanmean(pop.spikeCounts_choiceA,2)-nanmean(pop.spikeCounts_choiceB,2))';
    fprime_off_task = arrayfun(@(n_idx) TuningCurve_fPrime_At(pop.(fprime_curve){n_idx}, pop.Orientation + 45), 1:popsize);
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
if verbose, fprintf('Right plot: correlation as fn of offset\n'); end

offsets = -90:5:90;

all_correlations = zeros(1,length(offsets));
confidence_intervals = zeros(2, length(offsets));
for o_idx=1:length(offsets)
    if verbose, fprintf('\t%d/%d\n', o_idx, length(offsets)); end
    offset = offsets(o_idx);
    all_fprimes = zeros(n_neurons,1);
    all_ctdms = zeros(n_neurons,1);
    start_idx = 1;
    for p_idx=1:n_pops
        pop = pops_task(p_idx);
        popsize = length(pop.cellnos);
        choice_triggered_delta_means = (nanmean(pop.spikeCounts_choiceA,2)-nanmean(pop.spikeCounts_choiceB,2))';
        fprime_at_offset = arrayfun(@(n_idx) TuningCurve_fPrime_At(pop.(fprime_curve){n_idx}, pop.Orientation + offset), 1:popsize);
    
        end_idx = start_idx + popsize - 1;
        all_fprimes(start_idx:end_idx) = fprime_at_offset;
        all_ctdms(start_idx:end_idx) = choice_triggered_delta_means;
        start_idx = end_idx + 1;
    end
    % bootstrap confidence intervals, N=1000
    nboot=1000;
    
    % we occasionally get NaN in only some portion of the tuning curve
    %  when neurons aren't well-behaved; drop those values.
    valid_indices = ~isnan(all_fprimes) & ~isnan(all_ctdms);
    all_fprimes = all_fprimes(valid_indices);
    all_ctdms = all_ctdms(valid_indices);
        
    all_correlations(o_idx) = corr(all_fprimes, all_ctdms);
    confidence_intervals(:,o_idx) = bootci(nboot, {@corr, all_fprimes, all_ctdms}, ...
        'Options', statset('UseParallel', true));
end

subplot(1,3,3);
errorbar(offsets, all_correlations, confidence_intervals(1,:), confidence_intervals(2,:));
title(sprintf('Corr. f'' of tuning curve with choice-triggered means\nas a function of distance-from-trial-center'));
xlabel('offset from trial alignment');
ylabel('correlation of tuning curve f'' at \theta+offset with choice-triggered means');