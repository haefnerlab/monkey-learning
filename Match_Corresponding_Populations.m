function [ ptasks, pfixs ] = Match_Corresponding_Populations( ptasks, pfixs )
%MATCH_CORRESPONDING_POPULATIONS align 'task' and 'fixation' populations
%
% I.E. we keep only sessions and cells shared across both fixation and task
% datasets
%
% NOTE none of the 'additional' fields are filtered here like
% 'tuning_vm' or 'fprime_bestfit'. 
% RUN THIS FUNCTION BEFORE ANY OF THE OTHER PREPROCESSING

sessions_task = arrayfun(@(p) p.Header.SessionNumber, ptasks, ...
    'UniformOutput', false);
sessions_fix = arrayfun(@(p) p.Header.SessionNumber, pfixs, ...
    'UniformOutput', false);

%% find session names common to both
common_sessions = intersect(sessions_task, sessions_fix);
% filter each array down to those with common sessions, in same order
ptask_idxs = intersect_indexes(sessions_task, common_sessions);
pfix_idxs = intersect_indexes(sessions_fix, common_sessions);

ptasks = ptasks(ptask_idxs);
pfixs  =  pfixs(pfix_idxs);

%% find cellnos common to both
ptask_cells = arrayfun(@(p) p.cellnos, ptasks, 'UniformOutput', false);
pfix_cells = arrayfun(@(p) p.cellnos, pfixs, 'UniformOutput', false);

common_cells = arrayfun(@(p_idx) intersect(ptask_cells{p_idx}, pfix_cells{p_idx}), 1:length(ptasks), ...
    'UniformOutput', false);

ptasks_copy  = ptasks;
pfixs_copy = pfixs;

for p_idx=1:length(ptasks)
    % find corresponding indexes in cellnos arrays that match common_cells
    task_cell_idxs = intersect_indexes(ptasks_copy(p_idx).cellnos, common_cells{p_idx});
    fix_cell_idxs = intersect_indexes(pfixs_copy(p_idx).cellnos, common_cells{p_idx});
    
    % keep only data corresponding to shared cells
    ptasks_copy(p_idx).cellnos = ptasks_copy(p_idx).cellnos(task_cell_idxs);
    ptasks_copy(p_idx).spikeTimes = ptasks_copy(p_idx).spikeTimes(task_cell_idxs);
    ptasks_copy(p_idx).spikeCounts = ptasks_copy(p_idx).spikeCounts(task_cell_idxs,:);
    ptasks_copy(p_idx).tuning = ptasks_copy(p_idx).tuning(task_cell_idxs);
    
    pfixs_copy(p_idx).cellnos = pfixs_copy(p_idx).cellnos(fix_cell_idxs);
    pfixs_copy(p_idx).spikeTimes = pfixs_copy(p_idx).spikeTimes(fix_cell_idxs);
    pfixs_copy(p_idx).spikeCounts = pfixs_copy(p_idx).spikeCounts(fix_cell_idxs,:);
end

ptasks = ptasks_copy;
pfixs = pfixs_copy;

end

function [indexes] = intersect_indexes(A, B)
[~,indexes] = intersect(A,B);
end