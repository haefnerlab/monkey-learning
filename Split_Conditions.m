function [ populations ] = Split_Conditions( populations )
%SPLIT_CONDITIONS split each population's spike counts into seven
%conditions:
% - pop.spikeCounts_stimA = spike counts on trials with one sign stim
% - pop.spikeCounts_stimB = spike counts on trials with other sign stim
% - pop.spikeCounts_stim0A = spike counts on smallest nonzero signal-to-noise A
% - pop.spikeCounts_stim0B = spike counts on smallest nonzero signal-to-noise B
% - pop.spikeCounts_stim0 = spike counts on zero-stimulus trials
% - pop.spikeCounts_choiceA = zero-stimulus trials, but chose A
% - pop.spikeCounts_choiceB = zero-stimulus trials, but chose B

populations = arrayfun(@Split_Conditions_Single_Pop, populations);

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

end