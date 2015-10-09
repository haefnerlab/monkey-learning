%% TEST Deriv

params = [2.3071    6.5792    1.8156   22.4661];
vals = linspace(0,180,201);
vm = TuningCurves.vonMises(vals, params);
vmp = TuningCurves.vonMisesDeriv(vals, params);
vmd = diff(vm) / (vals(2)-vals(1));

plot(vals, vmp);
hold on;
plot(vals(1:end-1), vmd, '.');
hold off;
legend('analytic', 'diff-approximation');

%% TEST FIT
clearvars; close all;

pops_fix = Load_Fixation_Data('jbe');
pops_task = Load_Task_Data('jbe');

[pops_task, pops_fix] = Match_Corresponding_Populations(pops_task, pops_fix);

f = figure();

os = linspace(0,180,201);

n_neurons = sum(arrayfun(@(pop) length(pop.cellnos), pops_fix));
tun = zeros(1,n_neurons);
fve = zeros(1,n_neurons);
amp = zeros(1,n_neurons);
j = 1;

for pi=length(pops_fix):-1:1
    pop = pops_fix(pi);
    
    % condVecLabel is a cell array of strings identifying contents of each
    %   condVec column. Find which one is 'orientation'.
    orientation_condition = strcmp('orientation', pop.condVecLabel);
    
    orientations = pop.condVec(:,orientation_condition);
    for ni=1:length(pop.cellnos)
        fprintf('%d/%d\n', j, n_neurons);
        counts = pop.spikeCounts(ni,:);
        
        if numel(counts) == numel(orientations)
            
            [best, curve, best_map, worst, worst_map] = TuningCurves.fitVonMises(orientations, counts);
            disp(best);
%             disp(best_map);
%             disp(worst);
%             disp(worst_map);
            
            % scatter plot with tuning curve overlayed
            %             scatter(orientations, counts);
            %             hold on;
            %             plot(os, TuningCurves.vonMises(os, best), 'LineWidth', 2);
            %             plot(os, TuningCurves.vonMises(os, worst), 'LineStyle', '--');
            % plot where we thing preferred orientation is (black)
            %             plot([best(4) best(4)], [0, max(counts)], '--k')
            % neurons already have some estimated tuning from A.B. et al
            % (green)
            %             et = popst(pi).tuning(ni);
            %             if ~isnan(et)
            %                 plot([et et], [0, max(counts)], '--g')
            %             end
            %             axis([0,180,0,max(counts)+10]);
            %             hold off;
            %             title(sprintf('Population %d Neuron %d', pi, ni));
            
            amp(j) = best(2) * exp(best(3));
            fve(j) = Util.Variance_Explained(counts, TuningCurves.vonMises(orientations, best));
            tun(j) = pops_task(pi).tuning(ni);
            if isnan(fve(j))
                keyboard
            end
%             if (fve(j) < 0.2 && ~isnan(tun(j))) || isnan(tun(j))
%                 % look at ones without VE but with tuning well-defined?
%                 % scatter plot with tuning curve overlayed
%                 scatter(orientations, counts);
%                 hold on;
%                 plot(os, TuningCurves.vonMises(os, best), 'LineWidth', 2);
%                 plot(os, TuningCurves.vonMises(os, worst), 'LineStyle', '--');
%                 % plot where we thing preferred orientation is (black)
%                 plot([best(4) best(4)], [0, max(counts)], '--k')
%                 % neurons already have some estimated tuning from A.B. et al
%                 % (green)
%                 et = pops_task(pi).tuning(ni);
%                 if ~isnan(et)
%                     plot([et et], [0, max(counts)], '--g')
%                 else
%                     plot([best(4) best(4)], [0, max(counts)], '--r')
%                 end
%                 axis([0,180,0,max(counts)+10]);
%                 hold off;
%                 title(sprintf('Population %d Neuron %d', pi, ni));
%                 drawnow; pause;
%             end
            j = j+1;
        else
            warning('#counts was not same as #orientations');
        end
    end
    fprintf('end of population %d\n', pi);
end