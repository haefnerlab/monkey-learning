function [ parameters ] = New_Parameters( varargin )
%NEW_PARAMETERS Creates a 'parameters' struct that is used to control
% analyze_* functions.

% begin with default parameters
parameters = struct(...
    ... MONKEY ANALYSIS
    'monkey', 'lem', ...
    'recompute_tuning', false, ...
    'verbose', true, ...
    'fprime_curve', 'tuning_vm_curves', ... % or 'tuning_pw_curves'
    'min_pairs', 25, ... % minimum # shared trials between neurons, otherwise we throw out their correlations
    'min_rates', 10, ... % minimum average spike rate of a set of neurons for their moment to be calculated
    'moment', 2, ...
    'bootstrap', 1000,...
    'confidence', 0.95,... % confidence interval for plotting bootstrapped results
    'diagonal', false,... % whether to include self-correlations in analysis
    'num_offsets', 37, ... % for analyze_fprime_tuning_curves, number of 'hypothetical tasks' to test
    'collapse_offsets', true, ... % whether to collapse together rotationally symmetric offsets (which is mildly broken for odd moments)
    'corr_type', 'Pearson', ... % may also be 'Spearman' (i.e. the 'type' argument to the corr() function)
    ... SAMPLING MODEL
    'nc_tuning_method', 'tuning',...
    'discsize', 8, ...
    'sampling_npops', 20, ...
    'sampling_nneurons', 10, ...
    'sampling_delta', 0.08);

% add in any name, value pair from varargin
for i=1:2:length(varargin)
    if ~isfield(parameters, varargin{i})
        warning('Unkown parameter: %s', varargin{i});
    end
    parameters.(varargin{i}) = varargin{i+1};
end

end

