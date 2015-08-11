% f-prime refers to the derivative of the tuning function of a neuron with
% respect to the change in stimulus. We have two ways of computing it:
%
% in Compute_fPrime, we use the 'standard' method of taking the difference
%   in mean firing rates of a neuron under the 2 conditions. This requires
%   only the 'task' dataset.
%
% in Compute_fPrime_bestfit, we use the fixation dataset to fit a vonMises
%   tuning curve, then we use the difference in the fit tuning curve at the
%   two stimulus values.
%
% This script compares the two.

clearvars -except pops_task pops_fix; close all;

% monkey = 'jbe';
monkey = 'lem';

if ~exist('pops_task', 'var') || ~exist('pops_fix', 'var')
    pops_fix = Load_Fixation_Data(monkey);
    pops_task = Load_Task_Data(monkey);
    [pops_task, pops_fix] = Match_Corresponding_Populations(pops_task, pops_fix);
end

pops_task = Split_Conditions( pops_task );
pops_task = Compute_fPrime_stimulus_means( pops_task );
pops_task = Compute_fPrime_bestfit( pops_task, pops_fix );
pops_task = Compute_fPrime_fixation_means( pops_task, pops_fix );

colors = hsv(length(pops_task));

all_fprimes_stimulus_means = arrayfun(@(p) p.fprime_stimulus_mean, pops_task, 'UniformOutput', false);
all_fprimes_stimulus_means = horzcat(all_fprimes_stimulus_means{:});
all_fprimes_bestfit = arrayfun(@(p) p.fprime_bestfit, pops_task, 'UniformOutput', false);
all_fprimes_bestfit = horzcat(all_fprimes_bestfit{:});
all_fprimes_fixation_means = arrayfun(@(p) p.fprime_fixation_means, pops_task, 'UniformOutput', false);
all_fprimes_fixation_means = horzcat(all_fprimes_fixation_means{:});

[~,~,BigAx] = plotmatrix([all_fprimes_stimulus_means', all_fprimes_bestfit', all_fprimes_fixation_means']);

labels = {'Stimulus-means', 'best-fit von Mises', 'Fixation (linear interpolated)'};
ylabels = arrayfun(@(i) strjust(sprintf('%40s', labels{i}), 'center'), 3:-1:1, 'UniformOutput', false);
xlabels = arrayfun(@(i) strjust(sprintf('%40s', labels{i}), 'center'), 1:3,    'UniformOutput', false);
set(get(BigAx,'ylabel'),'string', [ylabels{:}]);
set(get(BigAx,'xlabel'),'string', [xlabels{:}]);
set(get(BigAx,'title'), 'string', 'scatter comparison of three f'' methods');  

% compare the three f' methods
% f = figure();
% y=x line
% plot([-40,40],[-40,40], '--', 'Color', 0.8*[1,1,1]);
% for p_idx=1:length(pops_task)
%     hold on;
%     scatter(pops_task(p_idx).fprime_stimulus_mean, pops_task(p_idx).fprime_bestfit, 8, colors(p_idx, :));
%     hold off;
% end
% xlabel('f'' (from stimulus response means)')
% ylabel('f'' (from fit tuning curve in fixation task)');
% title('Sanity-Check: f'' measured 2 different ways');
% axis('equal');
