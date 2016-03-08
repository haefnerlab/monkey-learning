function [ pops_task, pops_fix ] = Split_Conditions( pops_task, pops_fix )
%SPLIT_CONDITIONS split each population's spike Rates into seven
%conditions:
% - pop.spikeCounts_stimA = spike counts on trials with one sign stim
% - pop.spikeCounts_stimB = spike counts on trials with other sign stim
% - pop.spikeCounts_stim0A = spike counts on smallest nonzero signal-to-noise A
% - pop.spikeCounts_stim0B = spike counts on smallest nonzero signal-to-noise B
% - pop.spikeCounts_stim0 = spike counts on zero-stimulus trials
% - pop.spikeCounts_choiceA = zero-stimulus trials, but chose A
% - pop.spikeCounts_choiceB = zero-stimulus trials, but chose B
%
% ..likewise for each 'spikeRates' entry

pops_task = arrayfun(@Split_Conditions_Single_Pop, pops_task);
if nargin >= 2
    pops_fix = arrayfun(@Process_Fixation_Conditions, pops_fix);
end

end

function [ pop ] = Split_Conditions_Single_Pop( pop )

stimA = pop.condVec > 0;
stimB = pop.condVec < 0;
stim0 = pop.condVec == 0;
stim0A = pop.condVec == min(pop.condVec(stimA)); % smallest positive stim
stim0B = pop.condVec == max(pop.condVec(stimB)); % largest negative stim

pop.stimA = mod(pop.Orientation+90, 180);
pop.stimB = pop.Orientation;

choiceA = stim0 & pop.realChoice > 0;
choiceB = stim0 & pop.realChoice < 0;

pop.spikeCounts_stimA = pop.spikeCounts(:, stimA);
pop.spikeCounts_stimB = pop.spikeCounts(:, stimB);
pop.spikeCounts_stim0 = pop.spikeCounts(:, stim0);
pop.spikeCounts_stim0A = pop.spikeCounts(:, stim0A);
pop.spikeCounts_stim0B = pop.spikeCounts(:, stim0B);
pop.spikeCounts_choiceA = pop.spikeCounts(:, choiceA);
pop.spikeCounts_choiceB = pop.spikeCounts(:, choiceB);

% counts start after 50ms in and trialDur measured in ms
pop.spikeRates = pop.spikeCounts / (pop.trialDur/1000 - 0.05);
pop.spikeRates_stimA = pop.spikeCounts(:, stimA) / (pop.trialDur/1000 - 0.05);
pop.spikeRates_stimB = pop.spikeCounts(:, stimB) / (pop.trialDur/1000 - 0.05);
pop.spikeRates_stim0 = pop.spikeCounts(:, stim0) / (pop.trialDur/1000 - 0.05);
pop.spikeRates_stim0A = pop.spikeCounts(:, stim0A) / (pop.trialDur/1000 - 0.05);
pop.spikeRates_stim0B = pop.spikeCounts(:, stim0B) / (pop.trialDur/1000 - 0.05);
pop.spikeRates_choiceA = pop.spikeCounts(:, choiceA) / (pop.trialDur/1000 - 0.05);
pop.spikeRates_choiceB = pop.spikeCounts(:, choiceB) / (pop.trialDur/1000 - 0.05);

end

function [ pop ] = Process_Fixation_Conditions( pop )

% counts start after 50ms in
pop.spikeRates = pop.spikeCounts / (pop.trialDur/1000 - 0.05);
end