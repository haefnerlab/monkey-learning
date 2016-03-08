%% Use ML and user-feedback to estimate what priors on model parameters should look like
clearvars; close all;

nvotes = 100;

pops = Load_Fixation_Data('jbe');

tot_n_neurons = numel(horzcat(pops.cellnos));

% for uniform draws of pop,neuron pairs, keep one Nx2 list where
% [pop,neuron] = pop_and_neuron_idxs(idx,:)
pop_and_neuron_idxs = zeros(tot_n_neurons, 2);
first = 1;
for p_idx=1:length(pops)
    n_neurons_in_pop = length(pops(p_idx).cellnos);
    last = first + n_neurons_in_pop - 1;
    pop_and_neuron_idxs(first:last, 1) = p_idx;
    pop_and_neuron_idxs(first:last, 2) = 1:n_neurons_in_pop;
    first = last+1;
end

votes = zeros(nvotes,1);
fits = zeros(nvotes,4);

% Do nvotes random fits on random neurons
f = figure();
os = linspace(0,180,201);

for i=1:nvotes
    idxs = pop_and_neuron_idxs(randi(tot_n_neurons),:);
    pop = pops(idxs(1));

    % condVecLabel is a cell array of strings identifying contents of each
    %   condVec column. Find which one is 'orientation'.
    orientation_condition = find(strcmp('orientation', pop.condVecLabel));
    orientations = pop.condVec(:,orientation_condition);
    rates = pop.spikeRates(idxs(2),:);

    % Get a single randomly-initialized ML Fit (we want overfitting here
    % occasionally)
    [best, curve, best_map, worst, worst_map] = TuningCurves.fitVonMises(orientations, rates, false, 1);

    % scatter plot with tuning curve overlayed
    scatter(orientations, rates);
    hold on;
    plot(os, TuningCurves.vonMises(os, best), 'LineWidth', 2);
    axis([0,180,0,max(rates)+10]);
    hold off;
    title(sprintf('Population %d Neuron %d', idxs(1), idxs(2)));
    
    % keep track of fits
    fits(i,:) = best;
    choice = input('(bad fit, good fit, bad data, quit) ([0],1,2,3): ', 's');
    if strcmp(choice,'1')
        votes(i) = 1;
    elseif strcmp(choice, '2')
        votes(i) = 2;
    elseif strcmp(choice,'3')
        i = i-1;
        break;
    else
        votes(i) = 0;
    end
end

%% now we have i votes for things and can compare parameter values for
% 'good' vs 'bad' fits
votes = votes(1:i);
fits = fits(1:i,:);

binmethod = 'auto';

subplot(2,3,1);
[n,e] = histcounts(fits(votes==1,1), 'BinMethod', binmethod);
bar(e(1:end-1) + diff(e)/2,n,'b')
title('distribution of r_0 GOOD');
subplot(2,3,4);
[n,e] = histcounts(fits(votes==0,1), 'BinMethod', binmethod);
bar(e(1:end-1) + diff(e)/2,n,'r');
title('distribution of r_0 BAD');

subplot(2,3,2);
[n,e] = histcounts(fits(votes==1,2), 'BinMethod', binmethod);
bar(e(1:end-1) + diff(e)/2, n, 'b');
title('distribution of r_m_a_x GOOD');
subplot(2,3,5);
[n,e] = histcounts(fits(votes==0,2), 'BinMethod', binmethod);
bar(e(1:end-1) + diff(e)/2, n, 'r');
title('distribution of r_m_a_x BAD');

subplot(2,3,3);
[n,e] = histcounts(fits(votes==1,3), 'BinMethod', binmethod);
bar(e(1:end-1) + diff(e)/2, n, 'b');
title('distribution of k GOOD');
subplot(2,3,6);
[n,e] = histcounts(fits(votes==0,3), 'BinMethod', binmethod);
bar(e(1:end-1) + diff(e)/2, n, 'r');
title('distribution of k BAD');