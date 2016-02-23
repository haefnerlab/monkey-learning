function [] = plot_sampling_eigenvectors(e, standout, subsample)
%PLOT_SAMPLING_EIGENVECTORS recreate plots of eigenvectors from e.X, 


if nargin < 2, standout = 5; end
if nargin < 3, subsample = 1; end

[n_trials, n_neurons, n_samples] = size(e.X);

% scale down the timescale assuming 0.2 seconds per sample and a total
% trial length of 2 seconds
sim_seconds = n_samples * 0.2;
scale_by = 2 / sim_seconds;

[pref_orientations, reorder] = sort(e.Projection.phi_x);

% turn into trial x neurons matrix of spike count totals
spike_counts = scale_by * sum(e.X, 3);
spike_counts = spike_counts(:, reorder);

n_keep = n_neurons;
if subsample < 1
    n_keep = round(n_neurons * subsample);
    random_idxs = sort(randperm(n_neurons, n_keep));
    spike_counts = spike_counts(:, random_idxs);
    pref_orientations = pref_orientations(random_idxs);
end

[evecs, ~, evals] = pca(spike_counts);

figure();
plot(evals, 'o');
hold on;
colors = hsv(standout);
for i=1:standout
    plot(i, evals(i), 'o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
end
xlabel('rank');
ylabel('eigenvalue');
savefig(sprintf('sampling_eigenvalues_t%d_n%d.fig', n_trials, n_keep));

figure();
plot(pref_orientations, zeros(size(pref_orientations)), '--k');
hold on;
for i=1:standout
    plot(pref_orientations, evecs(:,i), '-', 'Color', colors(i,:));
end
axis tight;
set(gca, 'YLim', [-.1,.1]);
xlabel('preferred direction');
savefig(sprintf('sampling_eigenvectors_t%d_n%d.fig', n_trials, n_keep));
end

