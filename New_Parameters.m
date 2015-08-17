function [ parameters ] = New_Parameters( varargin )
%NEW_PARAMETERS Creates a 'parameters' struct that is used to control
% analyze_* functions.

% begin with default parameters
parameters = struct(...
    'monkey', 'lem', ...
    'verbose', true, ...
    'fprime_curve', 'tuning_vm_curves', ... % or 'tuning_pw_curves'
    'min_pairs', 25, ...
    'min_rates', 10, ...
    'moment', 3);

% add in any name, value pair from varargin
for i=1:2:length(varargin)
    parameters.(varargin{i}) = varargin{i+1};
end

end

