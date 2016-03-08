function [ pops_task ] = Compute_Choice_Probabilities( pops_task, zscore )
%COMPUTE_CHOICE_PROBABILITY compute a 'cp' field for each population: an
%array with one element per neuron, the choice probability of that neuron
%
%choice probability indirectly measures how correlated a single neuron is
%with the subject's *choice* at the end of the trial. It is measured during
%zero-stimulus conditions. Specifically, it is the area under the roc curve
%for spike rates to predict the choice

if nargin < 2, zscore = false; end

pops_task = arrayfun(@(p) Compute_Choice_Probabilities_Single_Pop(p, zscore), pops_task);

end

function [ pop ] = Compute_Choice_Probabilities_Single_Pop( pop, zscore )

n_neurons = length(pop.cellnos);

stim0 = pop.condVec == 0;
choices = pop.realChoice(stim0);
rates = pop.spikeRates(:, stim0);
pos = max(choices);

% % TODO
% if zscore
%     conditions = unique(pop.condVec);
%     stim0A = min(conditions(conditions > 0));
%     stim0B = max(conditions(conditions < 0));
%     
%     if stim0A
%         choices0A = pop.realChioce(pop.CondVec == stim0A);
%         spikes0A = 
%         
%     end
% end

pop.cp = zeros(1,n_neurons);
for i=1:n_neurons
    % manually handle spikeRates NaN
    valid = ~isnan(rates(i,:));
    if length(unique(choices(valid))) < 2
        pop.cp(i) = NaN;
    else
        [~,~,~,auroc] = perfcurve(choices(valid), rates(i,valid), pos);
        pop.cp(i) = auroc;
    end
end

end