function [ pops_task ] = Get_PiecewiseLinear_TuningCurves( pops_task, pops_fix )
%Get_PiecewiseLinear_TuningCurves use mean firing rates in pops_fix and linearly
% interpolate between them to make a piecewise tuning curve

for p_idx=1:length(pops_task)
    pt = pops_task(p_idx);
    pf = pops_fix(p_idx);
    
    orient_idx = strcmp('orientation', pf.condVecLabel);
    fixation_orientations = pf.condVec(:, orient_idx);
    
    usable_trials = ~isinf(fixation_orientations) & ~isnan(fixation_orientations);
    fixation_orientations = fixation_orientations(usable_trials);
    fixation_orientations = mod(fixation_orientations, 180);
    counts = pf.spikeCounts(:,usable_trials);
    
    % we make the tuning curve periodic by having 0 and 180 the same and by
    % using mod(o,180) in the anonymous function evaluation
    unique_orientations = vertcat(unique(fixation_orientations), 180);
        
    n_neurons = length(pt.cellnos);
    curves = cell(1,n_neurons);
    
    for n_idx=n_neurons:-1:1
        spikeCount_means = arrayfun(@(o) nanmean(counts(n_idx, fixation_orientations == mod(o,180))), unique_orientations);
        % use linear 'interp1' to interpolate
        curve = @(o) interp1(unique_orientations, spikeCount_means, mod(o,180), 'linear');
        curves{n_idx} = curve;
    end
    
    pops_task(p_idx).tuning_pw_curves = curves;
end
end
