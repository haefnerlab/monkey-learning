function [eigenvectors] = analyze_eigenvectors(params)
%ANALYZE_EIGENVECTORS get eigenvectors of zero-signal stimulus and plot as
%a function of preferred orientation - task orientation

%% Load and preprocess
[pops_task, ~] = Load_Preprocess(params);

%% PCA on each population

for p_idx=1:length(pops_task)
    %TODO 
    %  *  either remove nan or resample
    %  *  z-score
    % (?) deal with multiplicative term; z-scoring neurons may kill some
    % useful variation?
end
end