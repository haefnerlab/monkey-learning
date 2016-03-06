function [ pops_task, pops_fix ] = Compute_Sensitivity_Anova( pops_task, pops_fix )

for p_idx = 1:length(pops_task)
	pt = pops_task(p_idx);
	pf = pops_fix(p_idx);
	n_neurons = length(pt.cellnos);
	anova_pvalues = zeros(1,n_neurons);
	o_idx = strcmpi(pops_fix(p_idx).condVecLabel, 'orientation');
	for n_idx=1:n_neurons
		anova_pvalues(n_idx) = anova1(pf.spikeRates(n_idx,:)', pf.condVec(:,o_idx), 'off');
	end
	pops_task(p_idx).anova = anova_pvalues;
end

end