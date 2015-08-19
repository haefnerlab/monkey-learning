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

clearvars -except pops_task pops_fix params; close all;

if ~exist('params', 'var')
    params = New_Parameters('monkey', 'lem');
end

if ~exist('pops_task', 'var') || ~exist('pops_fix', 'var')
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

[~,Ax,BigAx] = plotmatrix([all_fprimes_stimulus_means', all_fprimes_bestfit', all_fprimes_fixation_means']);
set(get(BigAx,'title'), 'string', 'scatter comparison of three f'' methods');
set(get(Ax(1,1), 'ylabel'), 'string', 'Stimulus-driven');
set(get(Ax(2,1), 'ylabel'), 'string', '(high contrast) von Mises fit');
set(get(Ax(3,1), 'ylabel'), 'string', '(high contrast) means');
set(get(Ax(3,1), 'xlabel'), 'string', 'Stimulus-driven');
set(get(Ax(3,2), 'xlabel'), 'string', '(high contrast) von Mises fit');
set(get(Ax(3,3), 'xlabel'), 'string', '(high contrast) means');

