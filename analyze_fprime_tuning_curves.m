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

% calculate total number of data points we will get
n_momentdata = sum(arrayfun(@(p) length(Good_Pairs(p, params.diagonal)), pops_task));

% prepare scatter plot colors (neurons colored by population)
neuroncolors = zeros(n_momentdata,3);
start_idx = 1;
for p_idx=1:n_pops
    pop = pops_task(p_idx);
    popmoments = length(Good_Pairs(pop, params.diagonal));
    
    end_idx = start_idx + popmoments - 1;
    neuroncolors(start_idx:end_idx,:) = repmat(popcolors(p_idx,:),popmoments,1);
    start_idx = end_idx + 1;
end

%% Compute "choice-triggered" moment with bootstrapped estimates

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
        pairs = Good_Pairs(pop, params.diagonal);
        
        % get zero-stimulus noise correlations
        [noise_covariance, ~, ~, ~] = ...
            Util.nancomoment(pop.spikeCounts_stim0', params.moment, true, params.min_pairs, params.min_rates);
        % get *correlation* instead of covariance
        variances = diag(noise_covariance);
        noise_correlation = noise_covariance ./ sqrt(variances * variances');
        
        end_idx = start_idx + length(pairs) - 1;
        all_0stim{boot}(start_idx:end_idx) = noise_correlation(pairs);
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
        pairs = Good_Pairs(pop, params.diagonal);

        fprime_at_offset = arrayfun(@(n_idx) TuningCurve_fPrime_At(pop.(params.fprime_curve){n_idx}, pop.Orientation + offset), 1:popsize);
        
        [fprime_fprime, ~, ~, ~] = Util.nancomoment(fprime_at_offset, params.moment, true);
    
        end_idx = start_idx + length(pairs) - 1;
        all_fprimes(start_idx:end_idx, o_idx) = fprime_fprime(pairs);
        start_idx = end_idx + 1;
    end
end

%% get noise correlation vs f'f' at each 'task offset'
if params.verbose, fprintf('computing correlations..\n'); end

all_correlations = cell(params.bootstrap, 1);

for boot=1:params.bootstrap
    if params.verbose && mod(boot,10)==0, fprintf('\tboot %d/%d\n', boot, params.bootstrap); end
    all_correlations{boot} = zeros(1, length(offsets));
    boot_0stim = all_0stim(boot,:)';
    for o_idx=1:length(offsets)
        if params.verbose, fprintf('\t\t%d/%d\n', o_idx, length(offsets)); end
    
        fprimes_this_offset = all_fprimes(:, o_idx);
    
        valid_indices = ~isnan(fprimes_this_offset) & ~isnan(boot_0stim);
        all_correlations{boot}(o_idx) = corr(fprimes_this_offset(valid_indices), boot_0stim(valid_indices));
    end
end

all_correlations = vertcat(all_correlations{:});

% get mean and confidence intervals on these correlations
[mean_corr, lo_corr, hi_corr] = Util.meanci(all_correlations, params.confidence);
plus_corr = hi_corr - mean_corr;
minus_corr = mean_corr - lo_corr;

%% PLOT I: scatter and correlation when f' aligned with task
if params.verbose, fprintf('First plot: f'' task vs CT moment\n'); end

o_idx = offsets==0;

figure();
scatter(all_fprimes(:, o_idx), nanmean(all_0stim), 15, neuroncolors);
title(sprintf('Corr. choice-triggered moment vs task f''\nr=%.3f +/- (%.3e/%.3e)', ...
    mean_corr(o_idx), plus_corr(o_idx), minus_corr(o_idx)));
xlabel('f'' aligned to task')
ylabel('choice-triggered diff means')

%% PLOT II: scatter and correlation when f' 45 degrees off task
if params.verbose, fprintf('Second plot: f'' ortho vs CT moment\n'); end

o_idx = offsets==45;

figure();
scatter(all_fprimes(:, o_idx), nanmean(all_0stim), 15, neuroncolors);
title(sprintf('Corr. choice-triggered moment vs off-task f''\nr=%.3f +/- (%.3e/%.3e)', ...
    mean_corr(o_idx), plus_corr(o_idx), minus_corr(o_idx)));
xlabel('f'' aligned 45 degrees off task')
ylabel('choice-triggered diff means')

%% PLOT III: interpolate (I) and (II): correlation as a function of distance off task
if params.verbose, fprintf('Third plot: correlation as fn of offset\n'); end

figure();
Vis.boundedline(offsets, mean_corr, [minus_corr', plus_corr'], 'alpha');
title(sprintf('Correlations (f''f'' ~ noise correlation) as a function of task-offset'));
xlabel('offset from trial alignment');
ylabel('corr(f''f'',NC)');

end

function pairs = Good_Pairs(pop, diagonal)

if diagonal, k = 0;
else k = 1; end

n_neurons = length(pop.cellnos);
pair_indexes = find(triu(ones(n_neurons), k));
well_tuned_cells = ~isnan(pop.tuning)*1;
well_tuned_pairs = find(well_tuned_cells' * well_tuned_cells);

pairs = intersect(pair_indexes, well_tuned_pairs);

end