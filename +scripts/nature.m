function nature(m_idx)
% start by creating directories where figures will be saved

monkeys = {'both', 'lem', 'jbe'};
conditions = {'reject anova', 'reject fprime', 'keep all'};

savedir = 'paper figures';
if ~exist(savedir, 'dir'), mkdir(savedir); end

for m_dir = 1:length(monkeys)
    for c_idx = 1:length(conditions)
        newdir = fullfile(savedir, monkeys{m_dir}, conditions{c_idx});
        if ~exist(newdir, 'dir'), mkdir(newdir); end
    end
end

if nargin < 1, m_idx = 1; end

monkey = monkeys{m_idx};

params = New_Parameters(...
    'monkey', monkey, ...
    'moment', 1, ...
    'bootstrap', 500, ...
    'min_pairs', 25, ...
    'min_rates', 7, ...
    'exclusion_rule', 'anova', ...
    'exclusion_threshold', 0.05, ...
    'diagonal', false, ... % exclude all (i,i) pairs
    'num_offsets', 37, ... % every 2.5 degrees
    'collapse_offsets', false);

% scatter moments - ANOVA rejection rule
scripts.analyze_scatter_moments(params, fullfile(savedir, monkey, 'reject anova'));
% scatter moments - fprime rejection rule
params.exclusion_rule = 'fprime_pvalue';
scripts.analyze_scatter_moments(params, fullfile(savedir, monkey, 'reject fprime'));
% scatter moments - keep all
params.exclusion_threshold = inf;
scripts.analyze_scatter_moments(params, fullfile(savedir, monkey, 'keep all'));

close all;

% rotated task - 1st moment - all orientations - ANOVA rejection rule
params.exclusion_rule = 'anova';
params.exclusion_threshold = 0.05;
scripts.analyze_task_offset(params, fullfile('data', monkey, 'analyze_task_offset_moment1.mat'), false, fullfile(savedir, monkey, 'reject anova'));
% rotated task - 1st moment - all orientations - fprime rejection rule
params.exclusion_rule = 'fprime_pvalue';
scripts.analyze_task_offset(params, fullfile('data', monkey, 'analyze_task_offset_moment1.mat'), false, fullfile(savedir, monkey, 'reject fprime'));
% rotated task - 1st moment - all orientations - keep all
params.exclusion_threshold = inf;
scripts.analyze_task_offset(params, fullfile('data', monkey, 'analyze_task_offset_moment1.mat'), false, fullfile(savedir, monkey, 'keep all'));

close all;

% rotated task - 2nd moment - all orientations - ANOVA rejection rule
params.moment = 2;
params.collapse_offsets = false;
params.exclusion_rule = 'anova';
params.exclusion_threshold = 0.05;
scripts.analyze_task_offset(params, fullfile('data', monkey, 'analyze_task_offset_moment2.mat'), false, fullfile(savedir, monkey, 'reject anova'));
% rotated task - 2nd moment - all orientations - fprime rejection rule
params.exclusion_rule = 'fprime_pvalue';
scripts.analyze_task_offset(params, fullfile('data', monkey, 'analyze_task_offset_moment2.mat'), false, fullfile(savedir, monkey, 'reject fprime'));
% rotated task - 2nd moment - all orientations - keep all
params.exclusion_threshold = inf;
scripts.analyze_task_offset(params, fullfile('data', monkey, 'analyze_task_offset_moment2.mat'), false, fullfile(savedir, monkey, 'keep all'));

% rotated task - 2nd moment - collapsed orientations - ANOVA rejection rule
params.collapse_offsets = true;
params.exclusion_rule = 'anova';
params.exclusion_threshold = 0.05;
scripts.analyze_task_offset(params, fullfile('data', monkey, 'analyze_task_offset_moment2.mat'), false, fullfile(savedir, monkey, 'reject fprime'));
% rotated task - 2nd moment - collapsed orientations - fprime rejection rule
params.exclusion_rule = 'fprime_pvalue';
scripts.analyze_task_offset(params, fullfile('data', monkey, 'analyze_task_offset_moment2.mat'), false, fullfile(savedir, monkey, 'reject fprime'));
% rotated task - 2nd moment - collapsed orientations - keep all
params.exclusion_threshold = inf;
scripts.analyze_task_offset(params, fullfile('data', monkey, 'analyze_task_offset_moment2.mat'), false, fullfile(savedir, monkey, 'keep all'));

close all;

% f' linearity
test_fprime_methods;
fprintf('save %s fprime plots manually. (todo: automate)\n', monkey);
pause;

close all;

end
