function [ fig ] = PCA_2neuron( pop, fig )
%PCA_2neuron make 2D projection of neural responses showing graded
% stimulus levels with red->blue fade
%
% make sure Split_Conditions is run first so that condVec is sane

if nargin < 2, fig = figure(); end

stimuli = unique(pop.condVec);
% zero-stim handled separately
stimuli(stimuli == 0) = [];
nstim = length(stimuli);

zerosignal_choiceB = pop.condVec == 0 & pop.realChoice < 0;

% show stimulus trials red vs blue
blues = (sign(stimuli) + 1) / 2;
colors = zeros(nstim + 2, 3);
colors(1:nstim,:) = [ones(nstim,1)-blues, zeros(nstim,1), blues];
% make point sizes proportional to signal strength
sizes = zeros(nstim+2, 1);
sizes(1:nstim) = abs(stimuli) * 20 / max(abs(stimuli));
% assign classes to stimulus cases and to zerostimulus choice A
stimuli(end+1) = 0;
classes = arrayfun(@(s) find(stimuli == s, 1), pop.condVec);

% make zero-stimulus case different (size = 15, 2 classes based on choice)
sizes(nstim+1:end) = 15;
classes(zerosignal_choiceB) = nstim+2;

colors(nstim+1,:) = [1, 0.8, 0];
colors(nstim+2,:) = [0, 0.8, 1];

% show zero-stimulus trials yellow and cyan (i.e. red/blue + some green)


[eigenvectors, mean_vec, covariance, fig, handles] = Vis.PCA_Projection(pop.spikeRates', classes, colors, sizes, fig);
set(handles(nstim+1:end), 'Marker', 'diamond');

disp(mean_vec);
disp(eigenvectors);
disp(covariance);
end