function [all_correlations, all_pvalues, rot_sym_offsets] = analyze_task_offset( params, memo_file, recompute, figpath )
%ANALYZE_TASK_OFFSET complementary to analyze_scatter_moments,
% this function compares a "choice-triggered" (zero-stimulus) moment with
% f' for different notions of f' to see if the statistical moments of 
% choice-triggered responses are aligned to the task
%
% params.params.fprime_curve may be 'tuning_vm_curves' to use von Mises fits or
%   'tuning_pw_curves' to use piecewise-linear interpolations from the
%   fixation task as the tuning curve
%
% see New_Parameters() for other inputs that may be set

%% just load precomputed results, filter, and redo plots if possible

need_computation = true;
if nargin > 1
    if nargin < 3, recompute = false; end
    if ~recompute && exist(memo_file, 'file')
        fprintf('already got results. loading %s\n', memo_file);
        load(memo_file);
        need_computation = false;
    end
end
    
%% Load and preprocess
[pops_task, pops_fix] = Load_Preprocess(params);
    
n_pops = length(pops_task);
n_momentdata = sum(arrayfun(@(p) length(p.cellnos)^params.moment, pops_task));
% make array of offsets to be tested, ensuring that 0 and 45 are included
offsets = unique(horzcat(linspace(-45, 45, params.num_offsets), [0,45]));

% lambda function: get f' from a tuning curve at a given stimulus
TuningCurve_fPrime_At = @(curve, angle) (curve(angle) - curve(angle+90));

%% Loop over bootstrapped spike data to get 0stim and fprime moments as cell arrays
if need_computation
    all_0stim = cell(params.bootstrap,1); % each entry is (n_momentdata x 1) array
    all_fprimes = cell(params.bootstrap, 1); % each entry is (n_momentdata x n_offsets) matrix
    
    parfor boot=1:params.bootstrap
        if params.verbose && mod(boot,10)==0, fprintf('\tbootstrap loop %d/%d\n', boot, params.bootstrap); end
        
        % load precomputed bootstrapped data if they exist,
        % otherwise compute from scratch and save
        boot_file = fullfile('data', params.monkey, sprintf('boot_pop_%04d.mat', boot));
        if exist(boot_file, 'file')
            fprintf('loading %s\n', boot_file);
            precomputed = load(boot_file);
            bootpop = precomputed.pops_task;
            % TODO - remove this bit once all recomputed and saved
            if ~isfield(bootpop, 'anova')
                [bootpop, bootfix] = Split_Conditions(bootpop, precomputed.pops_fix);
                bootpop = Compute_Sensitivity_Anova(bootpop, bootfix);
            end
            if ~isfield(bootpop, 'fprime_pvalue')
                bootpop = Compute_fPrime_stimulus_means(bootpop,true);
            end
        else
            [bootpop, bootfix] = Bootstrap_TuningCurves(pops_task, pops_fix);
            [bootpop] = Bootstrap_SpikeCounts(bootpop);
            save_pops(boot_file, bootpop, bootfix); 
        end
        
        %% Compute "choice-triggered" or "zero-stimulus" moment with bootstrapped estimates
        
        boot_0stim = cell(1,n_pops);
        
        % loop over populations to fill in all_0stim array
        for p_idx=1:n_pops
            pop = bootpop(p_idx);
            
            % get moment (divided by sqrt(variances); this means we're looking
            % at correlations not covariances, etc)
            if params.moment == 1
                % use "choice-triggered" direction for first moment
                spikes_moment = (nanmean(pop.spikeRates_choiceA,2)-nanmean(pop.spikeRates_choiceB,2))';
                % ---
                % bugfix: what is called 'A' or 'B' is not consistently the positive
                % direction of the stimulus axis 's'. This causes a large fraction of
                % the choice_triggered_delta_means to have the wrong sign, which
                % completely ruins the correlation.. We correct for it here.
                signal_trials = pop.condVec > 0;
                % dot product between 'correctChoice' and 'sign(condVec)' will be
                % positive if they align and negative if they don't.
                choice_sign = sign(pop.correctChoice(signal_trials)' * sign(pop.condVec(signal_trials)));
                % ---
                spikes_moment = choice_sign * spikes_moment;
                spikes_moment = spikes_moment ./ sqrt(nanvar(pop.spikeRates_stim0', 1));
                % remove data where min_rates not satisfied
                below_rate_threshold = nanmean(pop.spikeRates_stim0,2)' <= params.min_rates;
                spikes_moment(below_rate_threshold) = NaN;
            elseif mod(params.moment,2) == 1
                % use outer product of choice-triggered direction for odd
                % moments
                % TODO - use min_pairs and min_rates here, do choice_sign
                spikes_diff = (nanmean(pop.spikeRates_choiceA,2)-nanmean(pop.spikeRates_choiceB,2))';
                spikes_moment = Util.nancomoment(spikes_diff, params.moment, true, true);
            else
                % for 2nd and higher even moments, get stats of all _stim0 spikes
                [spikes_moment, ~, ~, ~] = ...
                    Util.nancomoment(pop.spikeRates_stim0', params.moment, true, true, true, params.min_pairs);
            end
            
            boot_0stim{p_idx} = spikes_moment(:);
        end
        
        all_0stim{boot} = vertcat(boot_0stim{:});
        
        %% Compute f' moment at each of a range of offsets from the task
        all_fprimes{boot} = cell(n_pops,1);
        for p_idx=1:n_pops
            pop = bootpop(p_idx);
            % create n_moments_this_pop x n_offsets array (to be 'vertcat'd later)
            fprime_moment_each_offset = zeros(length(pop.cellnos)^params.moment, params.num_offsets);
            for o_idx=1:params.num_offsets
                
                if params.verbose, fprintf('\t\t[OR] %d/%d\n', o_idx, params.num_offsets); end
                offset = offsets(o_idx);
                
                fprime_at_offset = cellfun(@(curve) TuningCurve_fPrime_At(curve, pop.Orientation+offset), pop.(params.fprime_curve));
                
                fprime_moment = Util.ndouter(fprime_at_offset, params.moment);
                
                neuron_variances = nanvar(pop.spikeRates_stim0',1);
                fprime_moment = fprime_moment ./ reshape(Util.ndouter(sqrt(neuron_variances), params.moment), size(fprime_moment));
                
                fprime_moment_each_offset(:,o_idx) = fprime_moment(:);
            end
            all_fprimes{boot}{p_idx} = fprime_moment_each_offset;
        end
        all_fprimes{boot} = vertcat(all_fprimes{boot}{:});
    end
    
end

%% SAVE RESULTS if specified
if nargin > 1
    fprintf('saving to %s\n', memo_file);
    save(memo_file, 'all_fprimes', 'all_0stim');
end

%% set up task offsets, collapsing rotationally symmetric data together (only need 0:45)
    
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

%% Apply 'filters' to see which neurons/pairs should be kept
%but note that min_pairs and min_rates were already taken care of and can't
%easily be moved to this part of the analysis as post-hoc filters

filtered_indices = cell(n_pops,1);
for p_idx=1:n_pops
    pop = pops_task(p_idx);
    % the Good_Pairs function returns the indices that pass the
    % 'significance' and 'uniqueness' tests as specified in params
    idxs = Good_Pairs(pop, params);
    filtered_indices{p_idx} = false(length(pop.cellnos)^params.moment,1);
    filtered_indices{p_idx}(idxs) = true;
end
filtered_indices = vertcat(filtered_indices{:});

if params.verbose, fprintf('filter is keeping %d/%d\n', sum(filtered_indices), length(filtered_indices)); end

%% get correlation of (noise_moment vs f' moment) at each 'task offset'
% and taking into account which 'filters' should be applied

if params.verbose, fprintf('computing correlations..\n'); end

% both all_correlations and all_pvalues will turn into (boot x rot_sym_offsets)
% arrays after concatenation.
all_correlations = cell(params.bootstrap, 1);
all_pvalues = cell(params.bootstrap, 1);

parfor boot=1:params.bootstrap
    if params.verbose && mod(boot,10)==0, fprintf('\t[CORR] boot %d/%d\n', boot, params.bootstrap); end
    all_correlations{boot} = zeros(1, n_roffsets);
    boot_0stim = all_0stim{boot}; % (n_momentdata x 1) array
    boot_fprimes = all_fprimes{boot}; % (n_momentdata x n_offsets) matrix
    for o_idx=1:n_roffsets
        if params.verbose, fprintf('\t\t[CORR:OR] %d/%d\n', o_idx, n_roffsets); end
        
        % collapse together rotationally symmetric offsets
        collapse_indices = sym_offset(offsets)==rot_sym_offsets(o_idx);
        fprimes_this_offset = boot_fprimes(filtered_indices, collapse_indices);
        % ... into a single column vector
        fprimes_this_offset = fprimes_this_offset(:);
        
        % replicate the 0stim data as many times as there are redundant
        % offsets at o_idx (so that fprimes_this_offset and
        % extended_0stim have the same number of elements and can be
        % correlated)
        extended_0stim = repmat(boot_0stim(filtered_indices), n_redundant_offsets(o_idx), 1);
        
        % further reduce data - don't correlate NaN values
        valid_indices = ~isnan(fprimes_this_offset) & ~isnan(extended_0stim);
        [r,p] = corr(fprimes_this_offset(valid_indices), extended_0stim(valid_indices), 'type', params.corr_type);
        all_correlations{boot}(o_idx) = r;
        all_pvalues{boot}(o_idx) = p;
    end
end

all_correlations = vertcat(all_correlations{:});
all_pvalues = vertcat(all_pvalues{:});
    
% compute variance of 0stim moment across bootstrapping
var_0stim = nanvar(horzcat(all_0stim{:}), 1, 2);
    
% likewise get variance across bootstrapping of fprime moments
var_fprimes= squeeze(nanvar(cat(3,all_fprimes{:}), 1, 3));
    
% get mean and confidence intervals on correlations
[mean_corr, lo_corr, hi_corr] = Util.meanci(all_correlations, params.confidence);
plus_corr = hi_corr - mean_corr;
minus_corr = mean_corr - lo_corr;

% get mean and confidence intervals on p-values
[mean_pv, lo_pv, hi_pv] = Util.meanci(all_pvalues, params.confidence);
plus_pv = hi_pv - mean_pv;
minus_pv = mean_pv - lo_pv;

%% SCATTER PLOTS - commented out b/c redundant with analyze_scatter_moments
%%% PLOT I: scatter and correlation when f' aligned with task
% if params.verbose, fprintf('First plot: f'' task vs CT moment\n'); end
% 
% o_idx = rot_sym_offsets==0;
%             
% % collapse together rotationally symmetric offsets
% fprimes_this_offset = nanmean(all_fprimes(:, sym_offset(offsets)==rot_sym_offsets(o_idx),:),3);
% var_fprimes_this_offset = var_fprimes(:, sym_offset(offsets)==rot_sym_offsets(o_idx));
% % ... into a single column vector
% fprimes_this_offset = fprimes_this_offset(:);
% var_fprimes_this_offset = var_fprimes_this_offset(:);
% 
% extended_0stim = repmat(nanmean(all_0stim), 1, n_redundant_offsets(o_idx));
% extended_var0stim = repmat(var_0stim, n_redundant_offsets(o_idx), 1);
% extended_colors = repmat(neuroncolors, n_redundant_offsets(o_idx), 1);
% 
% figure();
% Vis.scatterr(fprimes_this_offset, extended_0stim, ...
%     var_fprimes_this_offset, extended_var0stim, 15, extended_colors);
% title(sprintf('Corr. noise moment vs task f''\nr=%.3f +/- (%.3e/%.3e)', ...
%     mean_corr(o_idx), plus_corr(o_idx), minus_corr(o_idx)));
% xlabel('f'' aligned to task')
% ylabel('choice-triggered diff means')
% 
% 
% if nargin >= 4
%     savefig(fullfile(figpath, 'scatter_aligned.fig'));
% end
% 
% %% PLOT II: scatter and correlation when f' 45 degrees off task
% if params.verbose, fprintf('Second plot: f'' ortho vs CT moment\n'); end
% 
% o_idx = rot_sym_offsets==45;
% 
% % collapse together rotationally symmetric offsets
% fprimes_this_offset = nanmean(all_fprimes(:, sym_offset(offsets)==rot_sym_offsets(o_idx),:),3);
% var_fprimes_this_offset = var_fprimes(:, sym_offset(offsets)==rot_sym_offsets(o_idx));
% % ... into a single column vector
% fprimes_this_offset = fprimes_this_offset(:);
% var_fprimes_this_offset = var_fprimes_this_offset(:);
% 
% extended_0stim = repmat(nanmean(all_0stim), 1, n_redundant_offsets(o_idx));
% extended_var0stim = repmat(var_0stim, n_redundant_offsets(o_idx), 1);
% extended_colors = repmat(neuroncolors, n_redundant_offsets(o_idx), 1);
% 
% figure();
% Vis.scatterr(fprimes_this_offset, extended_0stim, ...
%     var_fprimes_this_offset, extended_var0stim, 15, extended_colors);
% title(sprintf('Corr. noise moment vs off-task f''\nr=%.3f +/- (%.3e/%.3e)', ...
%     mean_corr(o_idx), plus_corr(o_idx), minus_corr(o_idx)));
% xlabel('f'' aligned 45 degrees off task')
% ylabel('choice-triggered diff means')
% 
% if nargin >= 4
%     savefig(fullfile(figpath, 'scatter_offset45.fig'));
% end
% 
%% PLOT III: correlation as a function of distance off task
if params.verbose, fprintf('Third plot: correlation as fn of offset\n'); end

figure();
Vis.boundedline(rot_sym_offsets, mean_corr, [minus_corr', plus_corr'], 'b', 'alpha');
title(sprintf('Correlations (f''f'' ~ noise correlation) as a function of task-offset'));
xlabel('offset from trial alignment');
ylabel('corr(f''f'',NC)');

if nargin >= 4
    savefig(fullfile(figpath, sprintf('[%s][%s]corr_vs_offset_m%d.fig', params.fprime_curve, params.corr_type, params.moment)));
end

%% PLOT IV: significance of 3rd plot's correlations as function of distance off task
if params.verbose, fprintf('Fourth plot: significance\n'); end

figure();
Vis.boundedline(rot_sym_offsets, mean_pv, [zeros(size(plus_pv')), plus_pv'], 'r', 'alpha');
set(gca, 'YScale', 'log');
title(sprintf('p-value of correlation as function of task offset'));
xlabel('offset from trial alignment');
ylabel('significance');

fprintf('correlation significance (bootstrap 95%%) = %f +%f -%f\n', mean_pv(1), plus_pv(1), minus_pv(1));

if nargin >= 4
    savefig(fullfile(figpath, sprintf('[%s][%s]significance_vs_offset_m%d.fig', params.fprime_curve, params.corr_type, params.moment)));
end

end

function save_pops(file, pops_task, pops_fix)
save(file, 'pops_task', 'pops_fix');
end