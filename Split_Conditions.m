function [ populations ] = Split_Conditions( populations )
%SPLIT_CONDITIONS split each population's spike counts into five
%conditions:
% - pop.spikeCounts_stimA = spike counts on trials with one sign stim
% - pop.spikeCounts_stimB = spike counts on trials with other sign stim
% - pop.spikeCounts_stim0 = spike counts on zero-stimulus trials
% - pop.spikeCounts_choiceA = zero-stimulus trials, but chose A
% - pop.spikeCounts_choiceB = zero-stimulus trials, but chose B

populations = arrayfun(@Split_Conditions_Single_Pop, populations);

end

function [ pop ] = Split_Conditions_Single_Pop( pop )

stimA = pop.condVec > 0;
stimB = pop.condVec < 0;
stim0 = pop.condVec == 0;

% some cases, sign(condVec) is opposite of correctChoice, but consistent
% for all trials
if ~all(sign(pop.condVec(~stim0)) == pop.correctChoice(~stim0))
    pop.condVec = -pop.condVec;
    stimA = pop.condVec > 0;
    stimB = pop.condVec < 0;
end
% sanity-check that it worked
assert(all(sign(pop.condVec(~stim0)) == pop.correctChoice(~stim0)));

choiceA = stim0 & pop.realChoice > 0;
choiceB = stim0 & pop.realChoice < 0;

pop.spikeCounts_stimA = pop.spikeCounts(:, stimA);
pop.spikeCounts_stimB = pop.spikeCounts(:, stimB);
pop.spikeCounts_stim0 = pop.spikeCounts(:, stim0);
pop.spikeCounts_choiceA = pop.spikeCounts(:, choiceA);
pop.spikeCounts_choiceB = pop.spikeCounts(:, choiceB);

end