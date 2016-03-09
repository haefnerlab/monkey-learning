function preboot(fro,count)

lem_params = New_Parameters('monkey', 'lem');
[orig_lem_pops_task, orig_lem_pops_fix] = Load_Preprocess(lem_params);
jbe_params = New_Parameters('monkey', 'jbe');
[orig_jbe_pops_task, orig_jbe_pops_fix] = Load_Preprocess(jbe_params);

nboot = 500;
if nargin < 1, fro=1; end
if nargin < 2, count=nboot-fro+1; end

parfor loop_idx=1:count
	b = fro+loop_idx-1;
	fprintf('BOOT %d/%d\n', b, nboot);
	if ~exist(sprintf('/scratch/rlange/data/monkey_data/both/boot_pop_%04d.mat', b), 'file')
		% bootstrap LEM data
		lem_pt = Bootstrap_SpikeCounts(orig_lem_pops_task);
		[lem_pt, lem_pf] = Bootstrap_TuningCurves(lem_pt, orig_lem_pops_fix);
		save_boot(b, 'lem', lem_pt, lem_pf);
		% bootstrap JBE data
		jbe_pt = Bootstrap_SpikeCounts(orig_jbe_pops_task);
		[jbe_pt, jbe_pf] = Bootstrap_TuningCurves(jbe_pt, orig_jbe_pops_fix);
		save_boot(b, 'jbe', jbe_pt, jbe_pf);
		% concatenate the two into BOTH data
		pops_task = horzcat(lem_pt, jbe_pt);
		pops_fix = horzcat(lem_pf, jbe_pf);
		save_boot(b, 'both', pops_task, pops_fix);
	else
		fprintf('~~skipping %04d~~\n', b);
	end
end

end
