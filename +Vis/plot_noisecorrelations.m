function [image_total, counts_total] = plot_noisecorrelations(pops_task, params, plot)

if nargin < 3, plot = true; end

%% get image that sums noise correlations plotted with preferred orientations

if plot, figure(); end

image_total = zeros(180);
counts_total = zeros(180);

for p_idx = 1:length(pops_task)
    pop = pops_task(p_idx);
    n_neurons = length(pop.cellnos);
    [noise_covariances,~,indices] = Util.nancomoment(pop.spikeCounts_stim0', 2, false, params.min_pairs, params.min_rates);

    % 3 flattened arrays of correlations for all valid pairs
    orientations_1 = zeros(length(indices),1);
    orientations_2 = zeros(length(indices),1);
    covariances = noise_covariances(indices);
    variances = diag(noise_covariances);
    denominator = sqrt(variances * variances');
    correlations = covariances ./ denominator(indices);

    for pair_idx = 1:length(indices)
        [n1, n2] = ind2sub([n_neurons, n_neurons], indices(pair_idx));
        orientations_1(pair_idx) = pop.(params.nc_tuning_method)(n1);
        orientations_2(pair_idx) = pop.(params.nc_tuning_method)(n2);
    end
    
    % take out neurons that don't have well-defined tuning
    invalid = isnan(orientations_1) | isnan(orientations_2);
    orientations_1 = orientations_1(~invalid);
    orientations_2 = orientations_2(~invalid);
    correlations = correlations(~invalid);
    
    % align to task so that a preferred direction of 0 means 'task-aligned'
    orientations_1 = orientations_1 - pop.Orientation;
    orientations_2 = orientations_2 - pop.Orientation;

    [im_this_pop, counts_this_pop] = Vis.image_pref_orientation(orientations_1, orientations_2, correlations, 180, params.discsize);
    counts_total = counts_total + counts_this_pop;
    countable = counts_this_pop > 0;
    image_total(countable) = image_total(countable) + im_this_pop(countable) .* counts_this_pop(countable);

    if plot
        Util.subplotsquare(length(pops_task), p_idx);
        Util.imagescnan(im_this_pop, [1 .6 .8]); colorbar;
        title(sprintf('Noise correlations pop %d', p_idx));
    end
end

%% plot totals

image_total(counts_total == 0) = NaN;
image_total = image_total ./ counts_total;

if plot
    figure();
    Util.imagescnan(image_total, [1 .6 .8]); colorbar;
    title('All noise correlations across populations (task=0,90)');
end

end