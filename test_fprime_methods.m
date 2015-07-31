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
pops_task = Compute_fPrime( pops_task );
pops_task = Compute_fPrime_bestfit( pops_task, pops_fix );

colors = hsv(length(pops_task));

% scatter fprime-according-to-means and fprime-according-to-fit-tuning
f = figure();
% y=x line
plot([-40,40],[-40,40], '--', 'Color', 0.8*[1,1,1]);
for p_idx=1:length(pops_task)
    hold on;
    scatter(pops_task(p_idx).fprime, pops_task(p_idx).fprime_bestfit, 8, colors(p_idx, :));
    hold off;
end
xlabel('f'' (from stimulus response means)')
ylabel('f'' (from fit tuning curve in fixation task)');
title('Sanity-Check: f'' measured 2 different ways');
axis('equal');