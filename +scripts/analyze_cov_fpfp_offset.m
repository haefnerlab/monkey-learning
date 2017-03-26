function [all_slopes, all_rsquared, rot_sym_offsets] = analyze_cov_fpfp_offset( params, memo_file, recompute, figpath )
%ANALYZE_COV_FPFP_OFFSET complementary to analyze_scatter_moments,
% this function regresses noise covariance to neurons' tuning derivative
% outer product. Based on analyze_task_offset
%
% params.fprime_curve may be 'tuning_vm_curves' to use von Mises fits or
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
% Get number of unique pairs in each population.
n_pairs = arrayfun(@(p) length(p.cellnos)^2, pops_task);
n_tot_pairs = sum(n_pairs);
% Make array of offsets to be tested, ensuring that 0 and 45 are included
offsets = unique(horzcat(linspace(-45, 45, params.num_offsets), [0,45]));

% lambda function: get f' from a tuning curve at a given stimulus
TuningCurve_fPrime_At = @(curve, angle) (curve(angle) - curve(angle+90));

%% Loop over bootstrapped spike data to get 0stim and fprime moments as cell arrays
if need_computation
    all_covs = cell(params.bootstrap,1); % each entry is (n_tot_pairs x 1) array
    all_fprimes = cell(params.bootstrap, 1); % each entry is (n_tot_pairs x n_offsets) matrix
    
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
        
        %% Compute noise covariance for this bootstrapped population.
        
        boot_covs = cell(1,n_pops);
        
        % loop over populations to fill in all_covs array
        for p_idx=1:n_pops
            pop = bootpop(p_idx);
            pop_covs = Util.nancomoment(pop.spikeRates_stim0', 2, true, true, false, params.min_pairs);
            boot_covs{p_idx} = pop_covs(:);
        end
        
        all_covs{boot} = vertcat(boot_covs{:});
        
        %% Compute f' moment at each of a range of offsets from the task
        all_fprimes{boot} = cell(n_pops,1);
        for p_idx=1:n_pops
            pop = bootpop(p_idx);
            % create n_pairs_this_pop x n_offsets array (to be 'vertcat'd later)
            fprime_moment_each_offset = zeros(n_pairs(p_idx), params.num_offsets);
            for o_idx=1:params.num_offsets
                
                if params.verbose, fprintf('\t\t[OR] %d/%d\n', o_idx, params.num_offsets); end
                offset = offsets(o_idx);
                
                fprime_at_offset = cellfun(@(curve) TuningCurve_fPrime_At(curve, pop.Orientation+offset), pop.(params.fprime_curve));
                
                fprime_outer = fprime_at_offset(:) * fprime_at_offset(:)';
                
                fprime_moment_each_offset(:,o_idx) = fprime_outer(:);
            end
            all_fprimes{boot}{p_idx} = fprime_moment_each_offset;
        end
        all_fprimes{boot} = vertcat(all_fprimes{boot}{:});
    end
end

%% SAVE RESULTS if specified
if nargin > 1
    fprintf('saving to %s\n', memo_file);
    save(memo_file, 'all_fprimes', 'all_covs');
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

if params.verbose, fprintf('filter is keeping %d/%d datapoints\n', sum(filtered_indices), length(filtered_indices)); end

%% get slope of (covariance vs f' moment) relation at each 'task offset'
% and taking into account which 'filters' should be applied

if params.verbose, fprintf('computing regressions..\n'); end

% both all_correlations and all_pvalues will turn into (boot x rot_sym_offsets)
% arrays after concatenation.
all_slopes = cell(params.bootstrap, 1);
all_rsquared = cell(params.bootstrap, 1);

parfor boot=1:params.bootstrap
    if params.verbose && mod(boot,10)==0, fprintf('\t[REG] boot %d/%d\n', boot, params.bootstrap); end
    all_slopes{boot} = zeros(1, n_roffsets);
    boot_covs = all_covs{boot}; % (n_pairs x 1) array
    boot_fprimes = all_fprimes{boot}; % (n_pairs x n_offsets) matrix
    for o_idx=1:n_roffsets
        if params.verbose, fprintf('\t\t[REG:OR] %d/%d\n', o_idx, n_roffsets); end
        
        % collapse together rotationally symmetric offsets
        collapse_indices = sym_offset(offsets)==rot_sym_offsets(o_idx);
        fprimes_this_offset = boot_fprimes(filtered_indices, collapse_indices);
        % ... into a single column vector
        fprimes_this_offset = fprimes_this_offset(:);
        
        % replicate the cov data as many times as there are redundant
        % offsets at o_idx (so that fprimes_this_offset and
        % extended_covs have the same number of elements and can be
        % correlated)
        extended_covs = repmat(boot_covs(filtered_indices), n_redundant_offsets(o_idx), 1);
        
        % further reduce data - don't regress NaN values
        valid_indices = ~isnan(fprimes_this_offset) & ~isnan(extended_covs);
        [s, r2] = linearfit(fprimes_this_offset(valid_indices), extended_covs(valid_indices));
        all_slopes{boot}(o_idx) = s;
        all_rsquared{boot}(o_idx) = r2;
    end
end

all_slopes = vertcat(all_slopes{:});
all_rsquared = vertcat(all_rsquared{:});

% get mean and confidence intervals on regression slopes
[mean_slope, lo_slope, hi_slope] = Util.meanci(all_slopes, params.confidence);
plus_slope = hi_slope - mean_slope;
minus_slope = mean_slope - lo_slope;

% get mean and confidence intervals on r2 values
[mean_r2, lo_r2, hi_r2] = Util.meanci(all_rsquared, params.confidence);
plus_r2 = hi_r2 - mean_r2;
minus_r2 = mean_r2 - lo_r2;

%% PLOT: regression as a function of task offset
if params.verbose, fprintf('Plot: correlation as fn of offset\n'); end

figure();
Vis.boundedline(rot_sym_offsets, mean_slope, [minus_slope', plus_slope'], 'b', 'alpha');
title(sprintf('Regression slopes (f''f'' ~ noise covariance) as a function of task-offset'));
xlabel('offset from trial alignment');
ylabel('regression slope');

if nargin >= 4
    savefig(fullfile(figpath, sprintf('[%s]slope_vs_offset.fig', params.fprime_curve)));
end

%% PLOT: R2 of previous plot's slopes
if params.verbose, fprintf('Fourth plot: significance\n'); end

figure();
Vis.boundedline(rot_sym_offsets, mean_r2, [zeros(size(plus_r2')), plus_r2'], 'r', 'alpha');
set(gca, 'YScale', 'log');
title(sprintf('R^2 of regression as function of task offset'));
xlabel('offset from trial alignment');
ylabel('R^2');

if nargin >= 4
    savefig(fullfile(figpath, sprintf('[%s]r2_vs_offset.fig', params.fprime_curve)));
end

end

function save_pops(file, pops_task, pops_fix)
save(file, 'pops_task', 'pops_fix');
end

function [slope, R2] = linearfit(X, Y)
slope = X\y;

SStotal = sum((Y - mean(Y)).^2);
SSresid = sum((X*slope - Y).^2);

R2 = 1 - SSresid / SStotal;
end