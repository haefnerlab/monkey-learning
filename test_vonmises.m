clearvars; close all;

pops = Load_Fixation_Data('jbe');

f = figure();

os = linspace(0,180,201);

for pi=1:length(pops)
    pop = pops(pi);
    
    % condVecLabel is a cell array of strings identifying contents of each
    %   condVec column. Find which one is 'orientation'.
    orientation_condition = find(strcmp('orientation', pop.condVecLabel));
    
    orientations = pop.condVec(:,orientation_condition);
    for ni=1:length(pop.cellnos)
        counts = pop.spikeCounts(ni,:);
        
        if numel(counts) == numel(orientations)
        
            [best, curve, best_map, worst, worst_map] = fitVonMises(orientations, counts);
            disp(best);
            disp(best_map);
            disp(worst);
            disp(worst_map);

            % scatter plot with tuning curve overlayed
            scatter(orientations, counts);
            hold on;
            plot(os, vonMises(os, best), 'LineWidth', 2);
            plot(os, vonMises(os, worst), 'LineStyle', '--');
            axis([0,180,0,max(counts)+10]);
            hold off;
            title(sprintf('Population %d Neuron %d', pi, ni));

            pause; 
        else
            warning('#counts was not same as #orientations');
        end
    end
    fprintf('end of population %d\n', pi);
end