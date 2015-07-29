function [ populations ] = Load_Fixation_Data( monkey_name )
%LOAD_TASK_DATA Load all .mat files for given monkey+condition into single
%struct array

datapath = fullfile('data', monkey_name, 'fixation');
pattern = fullfile(datapath, '*.mat');

files = dir(pattern);
size = length(files);

vars = {'Header', 'condVec', 'condVecLabel', 'trialDur', ...
    'trialStart', 'trialSeed', 'spikeTimes', 'spikeCounts', 'cellnos'};

for i=size:-1:1
    matdata = load(fullfile(datapath, files(i).name));
    for vi=1:length(vars)
        populations(i).(vars{vi}) = matdata.S.(vars{vi});
    end
end

end

