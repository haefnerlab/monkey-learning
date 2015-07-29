clearvars; close all;

pops = Load_Fixation_Data('jbe');

f = figure();

os = linspace(0,180,201);

for pi=2:length(pops)
    pop = pops(pi);
    
    % condVecLabel is a cell array of strings identifying contents of each
    %   condVec column. Find which one is 'orientation'.
    orientation_condition = find(strcmp('orientation', pop.condVecLabel));
    
    orientations = pop.condVec(:,orientation_condition);
    for ni=1:length(pop.cellnos)
        counts = pop.spikeCounts(ni,:);
        
        if numel(counts) == numel(orientations)
        
            [params, curve] = fitVonMises(orientations, counts);

            % scatter plot with tuning curve overlayed
            scatter(orientations, counts);
            hold on;
            plot(os, curve(os));
            hold off;
            title(sprintf('Population %d Neuron %d', pi, ni));

            pause; 
        end
    end
    fprintf('end of population %d\n', pi);
end