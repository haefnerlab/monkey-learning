function [ populations_task ] = Compute_fPrime_bestfit(populations_task, populations_fixation)
%Compute_fPrime_bestfit compute f' as change-in-tuning-curve-values at the
%  two stimulus values.
%
% note: Fit_Population_TuningCurves takes a long time. Use parpool if
% possible

if ~isfield(populations_task, 'tuning_vm_curves')
    populations_task = Fit_Population_TuningCurves(populations_task, populations_fixation);
end

populations_task = arrayfun(@Compute_fPrime_bestfit_Single_Pop, populations_task);

end

function [pop] = Compute_fPrime_bestfit_Single_Pop( pop )

stimA = pop.Orientation;
stimB = stimA + 90; % always orthogonal gratings

pop.fprime_bestfit = cellfun(@(curve) curve(stimA) - curve(stimB), pop.tuning_vm_curves);

end

