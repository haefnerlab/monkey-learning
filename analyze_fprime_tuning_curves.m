function analyze_fprime_tuning_curves( params )
%ANALYZE_FPRIME_TUNING_CURVES complementary to analyze_scatter_moments,
% this function compares a "choice-triggered" (zero-stimulus) moment with
% f' for different notions of f' to see if the statistical moments of 
% choice-triggered responses are aligned to the task
%
% params.params.fprime_curve may be 'tuning_vm_curves' to use von Mises fits or
%   'tuning_pw_curves' to use piecewise-linear interpolations from the
%   fixation task as the tuning curve

%% Load and preprocess

savefile = fullfile('data', params.monkey, 'preprocessed.mat');

% fitting tuning curves is time-consuming; precomputed results are put in
% a 'preprocessed.mat' file
if ~exist(savefile, 'file')
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
pops_task = Compute_fPrime_stimulus_means( pops_task );
pops_task = Compute_fPrime_bestfit( pops_task, pops_fix );
pops_task = Compute_fPrime_fixation_means( pops_task, pops_fix );

save(savefile, 'pops_task', 'pops_fix');

n_pops = length(pops_task);
n_neurons = sum(arrayfun(@(p) length(p.cellnos), pops_task));
popcolors = hsv(n_pops);

% anonymous function to get f' from a tuning curve at a given stimulus
TuningCurve_fPrime_At = @(curve, angle) (curve(angle) - curve(angle+90));

% prepare scatter plot colors (neurons colored by population)
neuroncolors = zeros(n_neurons,3);
start_idx = 1;
for p_idx=1:n_pops
    pop = pops_task(p_idx);
    popsize = length(pop.cellnos); 
    
    end_idx = start_idx + popsize - 1;
    neuroncolors(start_idx:end_idx,:) = repmat(popcolors(p_idx,:),popsize,1);
    start_idx = end_idx + 1;
end

%% Compute "choice-triggered" moment

if params.verbose, fprintf('computing ctdms..\n'); end

all_ctdms = zeros(n_neurons,1);

% loop over populations to fill in (flattened) all_ctdms array, where now
% each entry corresponds to a neuron
start_idx = 1;
for p_idx=1:n_pops
    pop = pops_task(p_idx);
    popsize = length(pop.cellnos);
    
    choice_triggered_delta_means = (nanmean(pop.spikeCounts_choiceA,2)-nanmean(pop.spikeCounts_choiceB,2))';
    
    end_idx = start_idx + popsize - 1;
    all_ctdms(start_idx:end_idx) = choice_triggered_delta_means;
    start_idx = end_idx + 1;
end

%% Compute f' moment at each of a range of offsets from the task

if params.verbose, fprintf('computing fprimes..\n'); end

% TODO: make # offsets controllable in params
offsets = unique(horzcat(-90:5:90, [0,45]));

all_fprimes = zeros(n_neurons, length(offsets)); % each column corresponds to one task-offset
for o_idx=1:length(offsets)
    if params.verbose, fprintf('\t%d/%d\n', o_idx, length(offsets)); end
    offset = offsets(o_idx);
    
    start_idx = 1;
    for p_idx=1:n_pops
        pop = pops_task(p_idx);
        popsize = length(pop.cellnos);

        fprime_at_offset = arrayfun(@(n_idx) TuningCurve_fPrime_At(pop.(params.fprime_curve){n_idx}, pop.Orientation + offset), 1:popsize);
    
        end_idx = start_idx + popsize - 1;
        all_fprimes(start_idx:end_idx, o_idx) = fprime_at_offset;
        start_idx = end_idx + 1;
    end
end

%% loop to get bootstrapped confidence intervals at each offset

if params.verbose, fprintf('computing correlations..\n'); end

all_correlations = zeros(1,length(offsets));
confidence_intervals = zeros(2, length(offsets));
for o_idx=1:length(offsets)
    if params.verbose, fprintf('\t%d/%d\n', o_idx, length(offsets)); end
    
    fprimes_this_offset = all_fprimes(:, o_idx);
    valid_indices = ~isnan(fprimes_this_offset) & ~isnan(all_ctdms);
        
    all_correlations(o_idx) = corr(fprimes_this_offset(valid_indices), all_ctdms(valid_indices));
    confidence_intervals(:,o_idx) = bootci(params.bootstrap, ...
        {@corr, fprimes_this_offset(valid_indices), all_ctdms(valid_indices)}, ...
        'Options', statset('UseParallel', true));
end

%% PLOT I: scatter and correlation when f' aligned with task
if params.verbose, fprintf('Left plot: f'' task vs CT moment\n'); end

o_idx = offsets==0;

figure();
scatter(all_fprimes(:, o_idx), all_ctdms, 15, neuroncolors);
title(sprintf('Corr. choice-triggered means vs task f''\nr=%.3f +/- %.3e', all_correlations(o_idx), diff(confidence_intervals(:,o_idx))));
xlabel('f'' aligned to task')
ylabel('choice-triggered diff means')

%% PLOT II: scatter and correlation when f' 45 degrees off task
if params.verbose, fprintf('Middle plot: f'' ortho vs CT moment\n'); end

o_idx = offsets==45;

figure();
scatter(all_fprimes(:, o_idx), all_ctdms, 15, neuroncolors);
title(sprintf('Corr. choice-triggered means vs off-task f''\nr=%.3f +/- %.3e', all_correlations(o_idx), diff(confidence_intervals(:,o_idx))));
xlabel('f'' aligned 45 degrees off task')
ylabel('choice-triggered diff means')

%% PLOT III: interpolate (I) and (II): correlation as a function of distance off task
if params.verbose, fprintf('Right plot: correlation as fn of offset\n'); end

figure();
errorbar(offsets, all_correlations, confidence_intervals(1,:), confidence_intervals(2,:));
title(sprintf('Corr. f'' of tuning curve with choice-triggered means\nas a function of distance-from-trial-center'));
xlabel('offset from trial alignment');
ylabel('correlation of tuning curve f'' at \theta+offset with choice-triggered means');

end