function [ populations_task ] = Fit_Population_TuningCurves( populations_task, populations_fixation )
%FIT_POPULATION_TUNINGCURVES use populations_fixation to fit a tuning
%curve, which will be saved into the populations_task.tuning_vm_curves cell
%array (of function handles) and tuning_vm preferred orientation values
%
% note: populations_task already comes with a field 'tuning'.. which is
% sometimes NaN. We don't discriminate well- from poorly-tuned neurons yet

for p_idx=1:length(populations_task)
    pfix = populations_fixation(p_idx);
    col_idx = strcmp('orientation', pfix.condVecLabel);
    orientations = pfix.condVec(:, col_idx);
    
    n_neurons = length(pfix.cellnos);
    tuning_curves = cell(1,n_neurons);
    preferred_directions = zeros(1,n_neurons);
    
    for n_idx=1:length(pfix.cellnos)
        fprintf('pop %d/%d neuron %d/%d\n', p_idx, length(populations_task), n_idx, length(pfix.cellnos));
        spikeCounts = pfix.spikeCounts(n_idx,:);
        [model, curve] = TuningCurves.fitVonMises(orientations, spikeCounts, false);
        tuning_curves{n_idx} = curve;
        preferred_directions(n_idx) = model(4);
    end
        
    populations_task(p_idx).tuning_vm_curves = tuning_curves;
    populations_task(p_idx).tuning_vm = preferred_directions;
end

end

