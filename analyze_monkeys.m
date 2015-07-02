% top-level script for monkey data analysis

%% Load and preprocess
if ~exist('populations', 'var')
    populations = Load_Data('lem');
end
populations = Split_Conditions( populations );
populations = Compute_fPrime( populations );

verbose = true;

%% Scatter plots of moments and 
figure();
colors = hsv(length(populations));

nmoments = 4;

for moment = 1:nmoments
    subplotsquare(nmoments, moment);
    title(sprintf('Moment %d', moment));
    
    if verbose, fprintf('Calculating moment %d\n', moment); end
    
    hold on;
    for pi=1:length(populations)
        if verbose, fprintf('\tPopulation %d of %d (%d neurons)\n', pi, length(populations), length(pop.cellnos)); end;
        pop = populations(pi);
        % get f'f'f'... up to moment times
        stimulus_moments = nancomoment(pop.fprime, moment);
        noise_moments = nancomoment(pop.spikeCounts_stim0', moment);
        
        scatter(stimulus_moments(:), noise_moments(:), 5, colors(pi,:));
    end
    hold off;
    
    xlabel('tuning curve f'' statistics');
    ylabel('zero-stimulus noise statistics');
end

clearvars -except populations verbose;