function test_fprime_methods(params, figdir)
% f-prime refers to the derivative of the tuning function of a neuron with
% respect to the change in stimulus. We have two ways of computing it:
%
% in Compute_fPrime_stimulus_means, we use the 'standard' method of fitting
%   a linear function to the responses in the 'task' dataset
%
% in Compute_fPrime_bestfit, we use the fixation dataset to fit a vonMises
%   tuning curve, then we use the difference in the fit tuning curve at the
%   two stimulus values.
%
% in Compute_fPrime_fixation_means, we create a piecewise linear function
% from the 'fixation' dataset to approximate the full tuning curve
%
% This script compares the different methods.

[pops_task, ~] = Load_Preprocess(params);

colors = hsv(length(pops_task));

all_fprimes_stimulus_means = arrayfun(@(p) p.fprime_stimulus_mean, pops_task, 'UniformOutput', false);
all_fprimes_stimulus_means = horzcat(all_fprimes_stimulus_means{:});
all_fprimes_bestfit = arrayfun(@(p) p.fprime_bestfit, pops_task, 'UniformOutput', false);
all_fprimes_bestfit = horzcat(all_fprimes_bestfit{:});
all_fprimes_fixation_means = arrayfun(@(p) p.fprime_fixation_means, pops_task, 'UniformOutput', false);
all_fprimes_fixation_means = horzcat(all_fprimes_fixation_means{:});

[~,Ax,BigAx] = plotmatrix([all_fprimes_stimulus_means', all_fprimes_bestfit', all_fprimes_fixation_means']);
set(get(BigAx,'title'), 'string', 'scatter comparison of three f'' methods');
set(get(Ax(1,1), 'ylabel'), 'string', 'Stimulus-driven');
set(get(Ax(2,1), 'ylabel'), 'string', '(high contrast) von Mises fit');
set(get(Ax(3,1), 'ylabel'), 'string', '(high contrast) means');
set(get(Ax(3,1), 'xlabel'), 'string', 'Stimulus-driven');
set(get(Ax(3,2), 'xlabel'), 'string', '(high contrast) von Mises fit');
set(get(Ax(3,3), 'xlabel'), 'string', '(high contrast) means');

[r,~] = corrcoef(all_fprimes_stimulus_means, all_fprimes_bestfit);
fprintf('stimulus means -- von mises fit:\tr=%f\n', r(2));

[r,~] = corrcoef(all_fprimes_stimulus_means, all_fprimes_fixation_means);
fprintf('stimulus means -- linear fit:\tr=%f\n', r(2));

[r,~] = corrcoef(all_fprimes_fixation_means, all_fprimes_bestfit);
fprintf('linear fit -- von mises fit:\tr=%f\n', r(2));

if nargin >= 2, savefig(fullfile(figdir, 'fprime_plotmatrix.fig')); end

%% now make a bunch of plots showing how well the 'stimulus means' method fits the data

figure();

total_neurons = numel([pops_task.cellnos]);
rmse = zeros(1,total_neurons);

i = 1;
for p_idx=1:length(pops_task)
    pop = pops_task(p_idx);
    stims = unique(pop.condVec);
    for n_idx=1:length(pop.cellnos)
        Util.subplotsquare(total_neurons,i);
        scatter(pop.condVec, pop.spikeRates(n_idx,:)');
        hold on;
        % redo linear fit (same as in stimulus-means function, but this
        % gets us the offset too)
        valid = ~isnan(pop.spikeRates(n_idx, :));
        stats = regstats(pop.spikeRates(n_idx,valid)', pop.condVec(valid), ...
            'linear', 'tstat');
        % plot linear fit
        fit_means = stats.tstat.beta(1) + stats.tstat.beta(2) * stims;
        plot(stims, fit_means);
        
        % compute RMSE
        means = zeros(size(stims));
        for s_idx=1:length(stims)
            means(s_idx) = nanmean(pop.spikeRates(n_idx, pop.condVec == stims(s_idx)));
        end
        rmse(i) = sqrt(sum((means - fit_means).^2));
        
        hold off;
        i = i+1;
    end
end

if nargin >= 2, savefig(fullfile(figdir, 'fprime_fits.fig')); end

figure();
hist(rmse, 50);
xlabel('rmse of f'' fit');
ylabel('count');

if nargin >= 2, savefig(fullfile(figdir, 'fprime_rmse.fig')); end

end