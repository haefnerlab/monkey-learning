function [all_correlations, all_pvalues, rot_sym_offsets] = analyze_fprime_tuning_curves( params, memo_file, recompute )
%ANALYZE_FPRIME_TUNING_CURVES complementary to analyze_scatter_moments,
% this function compares a "choice-triggered" (zero-stimulus) moment with
% f' for different notions of f' to see if the statistical moments of 
% choice-triggered responses are aligned to the task
%
% params.params.fprime_curve may be 'tuning_vm_curves' to use von Mises fits or
%   'tuning_pw_curves' to use piecewise-linear interpolations from the
%   fixation task as the tuning curve

%% just load precomputed results and remake plots if possible

need_computation = true;
if nargin > 1
    if nargin < 3, recompute = false; end
    if ~recompute && exist(memo_file, 'file')
        fprintf('already got results. loading %s\n', memo_file);
        load(memo_file);
        need_computation = false;
        % define things that would otherwise be skipped by need_comp=false
        if params.collapse_offsets
            sym_offset = @(o) min(abs(o), 90-abs(o));
        else
            sym_offset = @(o) o;
        end
    end
end
    
%% Load and preprocess
[pops_task, pops_fix] = Load_Preprocess(params);

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

if need_computation
    %% set up task offsets, collapsing rotationally symmetric data together (only need 0:45)
    
    % no matter what, make sure 0 and 45 are included (used in plots I and II)
    offsets = unique(horzcat(linspace(-45, 45, params.num_offsets), [0,45]));
    n_offsets = length(offsets);
    if params.collapse_offsets
        % [0,45] stays [0,45]
        % [0,-45] flips and becomes [0,45]
        % [45,90] is really getting "closer" to the task, so it's 90-[90:45]
        sym_offset = @(o) min(abs(o), 90-abs(o));
    else
        sym_offset = @(o) o;
    end
    
    rot_sym_offsets = unique(sym_offset(offsets));
    n_roffsets = length(rot_sym_offsets);
    
    % keep track of how many fprimes we are binning into a single
    % rotationally-symmetric one
    n_redundant_offsets = zeros(1,n_roffsets);
    
    for oi=1:n_roffsets
        n_redundant_offsets(oi) = sum(sym_offset(offsets) == rot_sym_offsets(oi));
    end
    
    %% Loop over bootstrapped spike data to get 0stim and fprime moments
    
    all_0stim = cell(params.bootstrap,1);
    all_fprimes = cell(params.bootstrap, 1);
    
    parfor boot=1:params.bootstrap
        if params.verbose && mod(boot,10)==0, fprintf('\tbootstrap loop %d/%d\n', boot, params.bootstrap); end
        
        % load precomputed results if they exist,
        % otherwise compute from scratch and save
        boot_file = fullfile('data', params.monkey, sprintf('boot_pop_%04d.mat', boot));
        if exist(boot_file, 'file')
            fprintf('loading %s\n', boot_file);
            precomputed = load(boot_file);
            bootpop = precomputed.pops_task;
        else
            [bootpop, bootfix] = Bootstrap_TuningCurves(pops_task, pops_fix);
            [bootpop] = Bootstrap_SpikeCounts(bootpop);
            save_pops(boot_file, bootpop, bootfix); 
        end
        
        %% Compute "choice-triggered" or "zero-stimulus" moment with bootstrapped estimates
        
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
                % remove data where min_rates not satisfied
                below_rate_threshold = nanmean(pop.spikeCounts_stim0,2)' <= params.min_rates;
                spikes_moment(below_rate_threshold) = NaN;
            elseif mod(params.moment,2) == 1
                % use outer product of choice-triggered direction for odd
                % moments
                % TODO - use min_pairs and min_rates here
                spikes_diff = (nanmean(pop.spikeCounts_choiceA,2)-nanmean(pop.spikeCounts_choiceB,2))';
                spikes_moment = Util.nancomoment(spikes_diff, params.moment, true, true);
            else
                % for 2nd and higher even moments, get stats of all _stim0 spikes
                [spikes_moment, ~, ~, ~] = ...
                    Util.nancomoment(pop.spikeCounts_stim0', params.moment, true, true, true, params.min_pairs, params.min_rates);
            end
            
            end_idx = start_idx + length(pairs) - 1;
            all_0stim{boot}(start_idx:end_idx) = spikes_moment(pairs);
            start_idx = end_idx + 1;
        end
        
        %% Compute f' moment at each of a range of offsets from the task
        all_fprimes{boot} = zeros(n_momentdata, n_offsets);
        for o_idx=1:n_offsets
            
            if params.verbose, fprintf('\t\t%d/%d\n', o_idx, n_offsets); end
            offset = offsets(o_idx);
            
            start_idx = 1;
            for p_idx=1:n_pops
                pop = bootpop(p_idx);
                popsize = length(pop.cellnos);
                pairs = Good_Pairs(pop, params.moment, params.diagonal);
                
                fprime_at_offset = arrayfun(@(n_idx) TuningCurve_fPrime_At(pop.(params.fprime_curve){n_idx}, pop.Orientation + offset), 1:popsize);
                
                if params.moment == 1
                    fprime_moment = fprime_at_offset;
                else
                    [fprime_moment, ~, ~, ~] = Util.nancomoment(fprime_at_offset, params.moment, true, false);
                end
                
                neuron_variances = nanvar(pop.spikeCounts_stim0',1);
                fprime_moment = fprime_moment ./ reshape(Util.ndouter(sqrt(neuron_variances), params.moment), size(fprime_moment));
                
                end_idx = start_idx + length(pairs) - 1;
                all_fprimes{boot}(start_idx:end_idx, o_idx) = fprime_moment(pairs);
                start_idx = end_idx + 1;
            end
        end
    end
    
    % collapse together data across bootstrapped trials
    all_0stim = vertcat(all_0stim{:});
    all_fprimes = cat(3, all_fprimes{:});
    
    %% get correlation of (noise_moment vs f' moment) at each 'task offset'
    if params.verbose, fprintf('computing correlations..\n'); end
    
    % both all_correlations and all_pvalues will turn into (boot x rot_sym_offsets)
    % arrays after concatenation.
    all_correlations = cell(params.bootstrap, 1);
    all_pvalues = cell(params.bootstrap, 1);
    
    parfor boot=1:params.bootstrap
        if params.verbose && mod(boot,10)==0, fprintf('\tboot %d/%d\n', boot, params.bootstrap); end
        all_correlations{boot} = zeros(1, n_roffsets);
        boot_0stim = all_0stim(boot,:)';
        boot_fprimes = squeeze(all_fprimes(:,:,boot));
        for o_idx=1:n_roffsets
            if params.verbose, fprintf('\t\t%d/%d\n', o_idx, n_roffsets); end
            
            % collapse together rotationally symmetric offsets
            fprimes_this_offset = boot_fprimes(:, sym_offset(offsets)==rot_sym_offsets(o_idx));
            % ... into a single column vector
            fprimes_this_offset = fprimes_this_offset(:);
            
            % replicate the 0stim data as many times as there are redundant
            % offsets at o_idx (so that fprimes_this_offset and
            % extended_0stim have the same number of elements and can be
            % correlated)
            extended_0stim = repmat(boot_0stim, n_redundant_offsets(o_idx), 1);
            
            valid_indices = ~isnan(fprimes_this_offset) & ~isnan(extended_0stim);
            [r,p] = corr(fprimes_this_offset(valid_indices), extended_0stim(valid_indices), 'type', params.corr_type);
            all_correlations{boot}(o_idx) = r;
            all_pvalues{boot}(o_idx) = p;
        end
    end
    
    all_correlations = vertcat(all_correlations{:});
    all_pvalues = vertcat(all_pvalues{:});
    
    %% SAVE RESULTS if specified
    if nargin > 1
        fprintf('saving to %s\n', memo_file);
        save(memo_file, 'params', 'all_fprimes', 'all_0stim', 'all_correlations', 'all_pvalues', 'rot_sym_offsets', 'offsets', 'n_redundant_offsets', 'n_offsets', 'n_roffsets');
    end

end
    
% compute variance of 0stim moment across bootstrapping
var_0stim = nanvar(all_0stim, 1, 1)';
    
% in correlations, fprime values will be weighted by their inverse
% variance across bootstrapping
var_fprimes= squeeze(nanvar(all_fprimes, 1, 3));
    
% get mean and confidence intervals on correlations
[mean_corr, lo_corr, hi_corr] = Util.meanci(all_correlations, params.confidence);
plus_corr = hi_corr - mean_corr;
minus_corr = mean_corr - lo_corr;

%% PLOT I: scatter and correlation when f' aligned with task
if params.verbose, fprintf('First plot: f'' task vs CT moment\n'); end

o_idx = rot_sym_offsets==0;
            
% collapse together rotationally symmetric offsets
fprimes_this_offset = nanmean(all_fprimes(:, sym_offset(offsets)==rot_sym_offsets(o_idx),:),3);
var_fprimes_this_offset = var_fprimes(:, sym_offset(offsets)==rot_sym_offsets(o_idx));
% ... into a single column vector
fprimes_this_offset = fprimes_this_offset(:);
var_fprimes_this_offset = var_fprimes_this_offset(:);

extended_0stim = repmat(nanmean(all_0stim), 1, n_redundant_offsets(o_idx));
extended_var0stim = repmat(var_0stim, n_redundant_offsets(o_idx), 1);
extended_colors = repmat(neuroncolors, n_redundant_offsets(o_idx), 1);

figure();
Vis.scatterr(fprimes_this_offset, extended_0stim, ...
    var_fprimes_this_offset, extended_var0stim, 15, extended_colors);
title(sprintf('Corr. noise moment vs task f''\nr=%.3f +/- (%.3e/%.3e)', ...
    mean_corr(o_idx), plus_corr(o_idx), minus_corr(o_idx)));
xlabel('f'' aligned to task')
ylabel('choice-triggered diff means')

figpath = fullfile('figures', params.monkey, sprintf('moment%d', params.moment));

savefig(fullfile(figpath, 'scatter_aligned.fig'));

%% PLOT II: scatter and correlation when f' 45 degrees off task
if params.verbose, fprintf('Second plot: f'' ortho vs CT moment\n'); end

o_idx = rot_sym_offsets==45;

% collapse together rotationally symmetric offsets
fprimes_this_offset = nanmean(all_fprimes(:, sym_offset(offsets)==rot_sym_offsets(o_idx),:),3);
var_fprimes_this_offset = var_fprimes(:, sym_offset(offsets)==rot_sym_offsets(o_idx));
% ... into a single column vector
fprimes_this_offset = fprimes_this_offset(:);
var_fprimes_this_offset = var_fprimes_this_offset(:);

extended_0stim = repmat(nanmean(all_0stim), 1, n_redundant_offsets(o_idx));
extended_var0stim = repmat(var_0stim, n_redundant_offsets(o_idx), 1);
extended_colors = repmat(neuroncolors, n_redundant_offsets(o_idx), 1);

figure();
Vis.scatterr(fprimes_this_offset, extended_0stim, ...
    var_fprimes_this_offset, extended_var0stim, 15, extended_colors);
title(sprintf('Corr. noise moment vs off-task f''\nr=%.3f +/- (%.3e/%.3e)', ...
    mean_corr(o_idx), plus_corr(o_idx), minus_corr(o_idx)));
xlabel('f'' aligned 45 degrees off task')
ylabel('choice-triggered diff means')

savefig(fullfile(figpath, 'scatter_offset45.fig'));

%% PLOT III: interpolate (I) and (II): correlation as a function of distance off task
if params.verbose, fprintf('Third plot: correlation as fn of offset\n'); end

figure();
Vis.boundedline(rot_sym_offsets, mean_corr, [minus_corr', plus_corr'], 'alpha');
title(sprintf('Correlations (f''f'' ~ noise correlation) as a function of task-offset'));
xlabel('offset from trial alignment');
ylabel('corr(f''f'',NC)');

savefig(fullfile(figpath, 'corr_vs_offset.fig'));

%% PLOT IV: significance of 3rd plot's correlations as function of distance off task
if params.verbose, fprintf('Fourth plot: significance\n'); end

figure();
mean_pvalue = mean(all_pvalues);
semilogy(rot_sym_offsets, mean_pvalue, '-r', 'LineWidth', 2);
title(sprintf('p-value as function of task offset'));
xlabel('offset from trial alignment');
ylabel('significance');

savefig(fullfile(figpath, 'significance.fig'));

end

function save_pops(file, pops_task, pops_fix)
save(file, 'pops_task', 'pops_fix');
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