function [ populations_task ] = Compute_fPrime_bestfit(populations_task, populations_fixation, recompute)
%Compute_fPrime_bestfit compute f' as change-in-tuning-curve-values at the
%  two stimulus values.
%
% note: Fit_Population_TuningCurves takes a long time. Use parpool if
% possible

if nargin < 3, recompute = false; end

if ~isfield(populations_task, 'tuning_vm_curves') || recompute
    populations_task = TuningCurves.Fit_Population_TuningCurves(populations_task, populations_fixation);
end

populations_task = arrayfun(@Compute_fPrime_bestfit_Single_Pop, populations_task);

end

function [pop] = Compute_fPrime_bestfit_Single_Pop( pop )

pop.fprime_bestfit = cellfun(@(curve) curve(pop.stimA) - curve(pop.stimB), pop.tuning_vm_curves);

end

