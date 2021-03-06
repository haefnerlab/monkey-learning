function [ pops_task ] = Compute_fPrime_fixation_means(pops_task, pops_fix, recompute)
%Compute_fPrime_fixation_means compute f' as from fixatation data using mean
% values from that task (or simple linear interpolation between them)
%
% (analogous to the other Compute_fprime_... functions)

if nargin < 3, recompute = false; end


if ~isfield(pops_task, 'tuning_pw_curves') || recompute
    pops_task = TuningCurves.Get_PiecewiseLinear_TuningCurves(pops_task, pops_fix);
end

for p_idx = 1:length(pops_task)
    n_neurons = length(pops_task(p_idx).cellnos);
    
    for n_idx=n_neurons:-1:1
        pw_curve = pops_task(p_idx).tuning_pw_curves{n_idx};
        pops_task(p_idx).fprime_fixation_means(n_idx) = pw_curve(pops_task(p_idx).stimA) - pw_curve(pops_task(p_idx).stimB);
    end
end

end