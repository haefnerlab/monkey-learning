function [ fig ] = PCA_2neuron( pop, fig )
%PCA_2neuron make 2D projection of neural responses showing graded
% stimulus levels with red->blue fade
%
% make sure Split_Conditions is run first so that condVec is sane

if nargin < 2, fig = figure(); end

stimuli = unique(pop.condVec);

nonzerosignal = pop.condVec ~= 0;
% zerosignal_choiceA = pop.condVec == 0 & pop.realChoice > 0;
% zerosignal_choiceB = pop.condVec == 0 & pop.realChoice < 0;

% show stimulus trials graded red->blue
blues = (sign(stimuli) + 1) / 2;
colors = [ones(length(stimuli),1)-blues, zeros(length(stimuli),1), blues];
sizes = abs(stimuli) * 20 / max(abs(stimuli));
classes = arrayfun(@(s) find(stimuli == s, 1), pop.condVec);
Vis.PCA_Projection(pop.spikeCounts(:,nonzerosignal)', classes(nonzerosignal), colors, sizes, fig);

% (TODO) show zero-stimulus trials

end