%% Load data
clearvars -except pops_task pops_fix params;

if ~exist('params', 'var')
    params = New_Parameters();
end

if ~exist('pops_task', 'var') || ~exist('pops_fix', 'var') || ~strcmp(pops_task(1).Header.monkey, params.monkey)
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
end

%% get image that sums noise correlations plotted with preferred orientations

fig_each = figure();

image_total = zeros(180);
counts_total = zeros(180);

for p_idx = 1:length(pops_task)
    pop = pops_task(p_idx);
    n_neurons = length(pop.cellnos);
    [noise_correlations,~,indices] = Util.nancomoment(pop.spikeCounts_stim0', 2, false, params.min_pairs, params.min_rates);

    % 3 flattened arrays of correlations for all valid pairs
    orientations_1 = zeros(length(indices),1);
    orientations_2 = zeros(length(indices),1);
    correlations = noise_correlations(indices);

    for pair_idx = 1:length(indices)
        [n1, n2] = ind2sub([n_neurons, n_neurons], indices(pair_idx));
        orientations_1(pair_idx) = pop.tuning(n1);
        orientations_2(pair_idx) = pop.tuning(n2);
    end
    
    % take out neurons that don't have well-defined tuning
    invalid = isnan(orientations_1) | isnan(orientations_2);
    orientations_1 = orientations_1(~invalid);
    orientations_2 = orientations_2(~invalid);
    correlations = correlations(~invalid);
    
    % align to task so that a preferred direction of 0 means 'task-aligned'
    orientations_1 = orientations_1 - pop.Orientation;
    orientations_2 = orientations_2 - pop.Orientation;

    [im_this_pop, counts_this_pop] = Vis.image_pref_orientation(orientations_1, orientations_2, correlations);
    counts_total = counts_total + counts_this_pop;
    countable = counts_this_pop > 0;
    image_total(countable) = image_total(countable) + im_this_pop(countable);

    Util.subplotsquare(length(pops_task), p_idx);
    Util.imagescnan(im_this_pop); colorbar;
    title(sprintf('Noise correlations pop %d', p_idx));
end

%% plot totals

image_total(counts_total == 0) = NaN;
image_total = image_total ./ counts_total;

fig_tot = figure();
Util.imagescnan(image_total); colorbar;
title('All noise correlations across populations (task=0,90)');