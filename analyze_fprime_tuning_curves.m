function [all_correlations, all_pvalues] = analyze_fprime_tuning_curves( params )
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

% calculate total number of data points we will get
n_momentdata = sum(arrayfun(@(p) length(Good_Pairs(p, params.moment, params.diagonal)), pops_task));

% prepare scatter plot colors (neurons colored by population)
neuroncolors = zeros(n_momentdata,3);
start_idx = 1;
for p_idx=1:n_pops
    pop = pops_task(p_idx);
    popmoments = length(Good_Pairs(pop, params.moment, params.diagonal));
    
    end_idx = start_idx + popmoments - 1;
    neuroncolors(start_idx:end_idx,:) = repmat(popcolors(p_idx,:),popmoments,1);
    start_idx = end_idx + 1;
end

%% Compute "choice-triggered" or "zero-stimulus" moment with bootstrapped estimates

if params.verbose, fprintf('computing noise correlations..\n'); end

all_0stim = cell(params.bootstrap,1);

parfor boot=1:params.bootstrap
    if params.verbose && mod(boot,10)==0, fprintf('\tbootstrapped spike counts %d/%d\n', boot, params.bootstrap); end
    bootpop = Bootstrap_SpikeCounts(pops_task);
    
    all_0stim{boot} = zeros(1,n_momentdata);
    
    % loop over populations to fill in (flattened) all_0stim array, where now
    % each entry corresponds to a neuron
    start_idx = 1;
    for p_idx=1:n_pops
        pop = bootpop(p_idx);
        pairs = Good_Pairs(pop, params.moment, params.diagonal);
        
        % get moment (divided by sqrt(variances); this means we're looking 
        % at correlations not covariances, etc)
        if params.moment == 1
            % use "choice-triggered" direction for first moment
            spikes_moment = (nanmean(pop.spikeCounts_choiceA,2)-nanmean(pop.spikeCounts_choiceB,2))';
            spikes_moment = spikes_moment ./ sqrt(nanvar(pop.spikeCounts_stim0', 1));
        else
            % for 2nd and higher moments, get stats of all _stim0 spikes
            [spikes_moment, ~, ~, ~] = ...
                Util.nancomoment(pop.spikeCounts_stim0', params.moment, true, true, true, params.min_pairs, params.min_rates);
        end
        
        end_idx = start_idx + length(pairs) - 1;
        all_0stim{boot}(start_idx:end_idx) = spikes_moment(pairs);
        start_idx = end_idx + 1;
    end
end

all_0stim = vertcat(all_0stim{:});

%% Compute f' moment at each of a range of offsets from the task

if params.verbose, fprintf('computing fprime moments..\n'); end

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
        pairs = Good_Pairs(pop, params.moment, params.diagonal);

        fprime_at_offset = arrayfun(@(n_idx) TuningCurve_fPrime_At(pop.(params.fprime_curve){n_idx}, pop.Orientation + offset), 1:popsize);
        
        [fprime_moment, ~, ~, ~] = Util.nancomoment(fprime_at_offset, params.moment, true, false);
    
        end_idx = start_idx + length(pairs) - 1;
        all_fprimes(start_idx:end_idx, o_idx) = fprime_moment(pairs);
        start_idx = end_idx + 1;
    end
end

%% get correlation of (noise_moment vs f' moment) at each 'task offset'
if params.verbose, fprintf('computing correlations..\n'); end

all_correlations = cell(params.bootstrap, 1);
all_pvalues = cell(params.bootstrap, 1);

for boot=1:params.bootstrap
    if params.verbose && mod(boot,10)==0, fprintf('\tboot %d/%d\n', boot, params.bootstrap); end
    all_correlations{boot} = zeros(1, length(offsets));
    boot_0stim = all_0stim(boot,:)';
    for o_idx=1:length(offsets)
        if params.verbose, fprintf('\t\t%d/%d\n', o_idx, length(offsets)); end
    
        fprimes_this_offset = all_fprimes(:, o_idx);
    
        valid_indices = ~isnan(fprimes_this_offset) & ~isnan(boot_0stim);
        [r,p] = corr(fprimes_this_offset(valid_indices), boot_0stim(valid_indices));
        all_correlations{boot}(o_idx) = r;
        all_pvalues{boot}(o_idx) = p;
    end
end

all_correlations = vertcat(all_correlations{:});
all_pvalues = vertcat(all_pvalues{:});

% get mean and confidence intervals on these correlations
[mean_corr, lo_corr, hi_corr] = Util.meanci(all_correlations, params.confidence);
plus_corr = hi_corr - mean_corr;
minus_corr = mean_corr - lo_corr;

%% PLOT I: scatter and correlation when f' aligned with task
if params.verbose, fprintf('First plot: f'' task vs CT moment\n'); end

o_idx = offsets==0;

figure();
scatter(all_fprimes(:, o_idx), nanmean(all_0stim), 15, neuroncolors);
title(sprintf('Corr. noise moment vs task f''\nr=%.3f +/- (%.3e/%.3e)', ...
    mean_corr(o_idx), plus_corr(o_idx), minus_corr(o_idx)));
xlabel('f'' aligned to task')
ylabel('choice-triggered diff means')

figpath = fullfile('figures', params.monkey, sprintf('moment%d', params.moment));

savefig(fullfile(figpath, 'scatter_aligned.fig'));

%% PLOT II: scatter and correlation when f' 45 degrees off task
if params.verbose, fprintf('Second plot: f'' ortho vs CT moment\n'); end

o_idx = offsets==45;

figure();
scatter(all_fprimes(:, o_idx), nanmean(all_0stim), 15, neuroncolors);
title(sprintf('Corr. noise moment vs off-task f''\nr=%.3f +/- (%.3e/%.3e)', ...
    mean_corr(o_idx), plus_corr(o_idx), minus_corr(o_idx)));
xlabel('f'' aligned 45 degrees off task')
ylabel('choice-triggered diff means')

savefig(fullfile(figpath, 'scatter_offset45.fig'));

%% PLOT III: interpolate (I) and (II): correlation as a function of distance off task
if params.verbose, fprintf('Third plot: correlation as fn of offset\n'); end

figure();
Vis.boundedline(offsets, mean_corr, [minus_corr', plus_corr'], 'alpha');
title(sprintf('Correlations (f''f'' ~ noise correlation) as a function of task-offset'));
xlabel('offset from trial alignment');
ylabel('corr(f''f'',NC)');

savefig(fullfile(figpath, 'corr_vs_offset.fig'));

%% PLOT IV: significance of 3rd plot's correlations as function of distance off task
if params.verbose, fprintf('Fourth plot: significance\n'); end

figure();
mean_pvalue = mean(all_pvalues);
sign_r = sign(mean_corr);
plot(offsets, mean_pvalue.*sign_r, '-r', 'LineWidth', 2);
title(sprintf('p*sign(r) as function of task offset'));
xlabel('offset from trial alignment');
ylabel('significance');

savefig(fullfile(figpath, 'significance.fig'));

end

function pairs = Good_Pairs(pop, moment, diagonal)
%GOOD_PAIRS finds indexes (i,j,k,..) where
% 1) they form a unique set (i.e. 1,2 and 2,1 not both included) and
% 2) all neurons i, j, and k have well-defined tuning.
%
% "pairs" is a bit of a misnomer; with moment=3 it's triples, etc.

if diagonal, k = 0;
else k = 1; end

% part 1: get unique pairs
n_neurons = length(pop.cellnos);
unq_pairs = find(Util.ndtriu(n_neurons*ones(1,moment), k));

% part 2: find where ~isnan(i) and ~isnan(j), etc.
well_tuned = ~isnan(pop.tuning)*1;

% now find where "all cells i,j,k,... had good tuning" by multiplying the
% logical array well_tuned by itself using a generalization of outer 
% product to n-dimensions (and there are 'moment' dimensions here)
if moment > 1
    % make a cell array, each entry containing a copy of well_tuned (which
    % is a row-vector)
    copies_of_well_tuned = mat2cell(repmat(well_tuned, moment, 1), ones(1,moment));
    % expand copies out into ndgrid (n-dimensional generalization of
    % meshgrid)
    [grids{1:moment}] = ndgrid(copies_of_well_tuned{:});
    
    grid = grids{1};
    for i=2:moment
        grid = grid .* grids{i};
    end
    well_tuned_idxs = find(grid);
else
    well_tuned_idxs = find(~isnan(well_tuned));
end

% end result is the intersection of the two sets of indexes
pairs = intersect(unq_pairs, well_tuned_idxs);

end