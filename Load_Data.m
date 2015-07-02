function [ populations ] = Load_Data( monkey_name )
%LOAD_DATA Load all .mat files for given monkey+condition into single
%struct array

datapath = fullfile('data', monkey_name, 'task');
pattern = fullfile(datapath, '*.mat');

files = dir(pattern);
size = length(files);

vars = {'Header', 'condVec', 'condVecLabel', 'Orientation', 'trialDur', ...
    'trialStart', 'correctChoice', 'realChoice', 'reward', 'trialSeed', ...
    'spikeTimes', 'spikeCounts', 'tuning', 'cellnos'};

for i=size:-1:1
    matdata = load(fullfile(datapath, files(i).name));
    for vi=1:length(vars)
        populations(i).(vars{vi}) = matdata.S.(vars{vi});
    end
end

end

