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

TuningCurve_fPrime_At = @(curve, angle) (curve(angle) - curve(angle+90));

%% PLOT I: scatter and correlation when f' aligned with task
if params.verbose, fprintf('Left plot: f'' task vs CT moment\n'); end
all_fprimes = zeros(n_neurons,1);
all_ctdms = zeros(n_neurons,1);
start_idx = 1;

figure();
hold on;
for p_idx=1:n_pops
    pop = pops_task(p_idx);
    popsize = length(pop.cellnos);
    choice_triggered_delta_means = (nanmean(pop.spikeCounts_choiceA,2)-nanmean(pop.spikeCounts_choiceB,2))';
    fprime_at_task = arrayfun(@(n_idx) TuningCurve_fPrime_At(pop.(params.fprime_curve){n_idx}, pop.Orientation), 1:popsize);
    
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

%% PLOT II: scatter and correlation when f' 45 degrees off task
if params.verbose, fprintf('Middle plot: f'' ortho vs CT moment\n'); end
all_fprimes = zeros(n_neurons,1);
all_ctdms = zeros(n_neurons,1);
start_idx = 1;

figure();
hold on;
for p_idx=1:n_pops
    pop = pops_task(p_idx);
    popsize = length(pop.cellnos);
    choice_triggered_delta_means = (nanmean(pop.spikeCounts_choiceA,2)-nanmean(pop.spikeCounts_choiceB,2))';
    fprime_off_task = arrayfun(@(n_idx) TuningCurve_fPrime_At(pop.(params.fprime_curve){n_idx}, pop.Orientation + 45), 1:popsize);
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

%% PLOT III: interpolate (I) and (II): correlation as a function of distance off task
if params.verbose, fprintf('Right plot: correlation as fn of offset\n'); end

offsets = -90:5:90;

all_correlations = zeros(1,length(offsets));
confidence_intervals = zeros(2, length(offsets));
for o_idx=1:length(offsets)
    if params.verbose, fprintf('\t%d/%d\n', o_idx, length(offsets)); end
    offset = offsets(o_idx);
    all_fprimes = zeros(n_neurons,1);
    all_ctdms = zeros(n_neurons,1);
    start_idx = 1;
    for p_idx=1:n_pops
        pop = pops_task(p_idx);
        popsize = length(pop.cellnos);
        choice_triggered_delta_means = (nanmean(pop.spikeCounts_choiceA,2)-nanmean(pop.spikeCounts_choiceB,2))';
        fprime_at_offset = arrayfun(@(n_idx) TuningCurve_fPrime_At(pop.(params.fprime_curve){n_idx}, pop.Orientation + offset), 1:popsize);
    
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

figure();
errorbar(offsets, all_correlations, confidence_intervals(1,:), confidence_intervals(2,:));
title(sprintf('Corr. f'' of tuning curve with choice-triggered means\nas a function of distance-from-trial-center'));
xlabel('offset from trial alignment');
ylabel('correlation of tuning curve f'' at \theta+offset with choice-triggered means');

end