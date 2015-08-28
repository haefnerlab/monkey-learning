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
popcolors = hsv(n_pops);

% lambda function: get f' from a tuning curve at a given stimulus
TuningCurve_fPrime_At = @(curve, angle) (curve(angle) - curve(angle+90));

% lambda function: number of non-duplicate moment data
Num_Unique_Moments = @(variables, moment) length(find(Util.ndtriu(variables * ones(1,moment))));

% calculate total number of data points we will get
n_momentdata = sum(arrayfun(@(p) Num_Unique_Moments(length(p.cellnos), params.moment), pops_task));

% prepare scatter plot colors (neurons colored by population)
neuroncolors = zeros(n_momentdata,3);
start_idx = 1;
for p_idx=1:n_pops
    pop = pops_task(p_idx);
    popsize = length(pop.cellnos); 
    popmoments = Num_Unique_Moments(popsize, params.moment);
    
    end_idx = start_idx + popmoments - 1;
    neuroncolors(start_idx:end_idx,:) = repmat(popcolors(p_idx,:),popmoments,1);
    start_idx = end_idx + 1;
end

%% Compute "choice-triggered" moment

if params.verbose, fprintf('computing ctdms..\n'); end

all_ctdms = zeros(n_momentdata, 1);

% loop over populations to fill in (flattened) all_ctdms array, where now
% each entry corresponds to a neuron
start_idx = 1;
for p_idx=1:n_pops
    pop = pops_task(p_idx);
    popsize = length(pop.cellnos);
    
    % by 'first moment' we really mean (choiceA mean) - (choiceB mean)
    % whereas 2nd+higher moments look at entire _stim0 distribution
    % together
    if params.moment == 1
        all_indices = 1:popsize;
        choice_triggered_moment = (nanmean(pop.spikeCounts_choiceA,2)-nanmean(pop.spikeCounts_choiceB,2))';
    else
        [choice_triggered_moment, ~, ~, all_indices] = ...
            Util.nancomoment(pop.spikeCounts_stim0', params.moment, true, params.min_pairs, params.min_rates);
        % note that anywhere min_pairs isn't satisfied, ctm is NaN
    end
    
    end_idx = start_idx + length(all_indices) - 1;
    all_ctdms(start_idx:end_idx) = choice_triggered_moment(all_indices);
    start_idx = end_idx + 1;
end

%% Compute f' moment at each of a range of offsets from the task

if params.verbose, fprintf('computing fprimes..\n'); end

% TODO: make # offsets controllable in params
% no matter what, make sure 0 and 45 are included (used in plots I and II)
offsets = unique(horzcat(-90:5:90, [0,45]));

all_fprimes = zeros(n_momentdata, length(offsets)); % each column corresponds to one task-offset
for o_idx=1:length(offsets)
    if params.verbose, fprintf('\t%d/%d\n', o_idx, length(offsets)); end
    offset = offsets(o_idx);
    
    start_idx = 1;
    for p_idx=1:n_pops
        pop = pops_task(p_idx);
        popsize = length(pop.cellnos);

        fprime_at_offset = arrayfun(@(n_idx) TuningCurve_fPrime_At(pop.(params.fprime_curve){n_idx}, pop.Orientation + offset), 1:popsize);
        
        % see comments in analogous part of previous section (~line 76)
        if params.moment == 1
            all_indices = 1:popsize;
            fprime_moment = fprime_at_offset;
        else
            [fprime_moment, ~, ~, all_indices] = Util.nancomoment(fprime_at_offset, params.moment, true);
        end
    
        end_idx = start_idx + length(all_indices) - 1;
        all_fprimes(start_idx:end_idx, o_idx) = fprime_moment(all_indices);
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
title(sprintf('Corr. choice-triggered moment vs task f''\nr=%.3f +/- %.3e', all_correlations(o_idx), diff(confidence_intervals(:,o_idx))));
xlabel('f'' aligned to task')
ylabel('choice-triggered diff means')

%% PLOT II: scatter and correlation when f' 45 degrees off task
if params.verbose, fprintf('Middle plot: f'' ortho vs CT moment\n'); end

o_idx = offsets==45;

figure();
scatter(all_fprimes(:, o_idx), all_ctdms, 15, neuroncolors);
title(sprintf('Corr. choice-triggered moment vs off-task f''\nr=%.3f +/- %.3e', all_correlations(o_idx), diff(confidence_intervals(:,o_idx))));
xlabel('f'' aligned 45 degrees off task')
ylabel('choice-triggered diff means')

%% PLOT III: interpolate (I) and (II): correlation as a function of distance off task
if params.verbose, fprintf('Right plot: correlation as fn of offset\n'); end

figure();
errorbar(offsets, all_correlations, confidence_intervals(1,:), confidence_intervals(2,:));
title(sprintf('Corr. f'' of tuning curve with choice-triggered moment\nas a function of distance-from-trial-center'));
xlabel('offset from trial alignment');
ylabel('correlation of tuning curve f'' at \theta+offset with choice-triggered moment');

end