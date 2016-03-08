%% TEST Deriv

curve = [2.3071    6.5792    1.8156   22.4661];
vals = linspace(0,180,201);
vm = TuningCurves.vonMises(vals, curve);
vmp = TuningCurves.vonMisesDeriv(vals, curve);
vmd = diff(vm) / (vals(2)-vals(1));

plot(vals, vmp);
hold on;
plot(vals(1:end-1), vmd, '.');
hold off;
legend('analytic', 'diff-approximation');

%% TEST FIT
clearvars -except params; close all;

if ~exist('params', 'var'), params = New_Parameters('monkey', 'jbe'); end

[pops_task, pops_fix] = Load_Preprocess(params);

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
        rates = pop.spikeRates(ni,:);
        
        if numel(rates) == numel(orientations)
            
            [best, curve, best_map, worst, worst_map] = TuningCurves.fitVonMises(orientations, rates);
            disp(best);
            disp(best_map);
            disp(worst);
            disp(worst_map);
            
            % scatter plot with tuning curve overlayed
            scatter(orientations, rates);
            hold on;
            plot(os, TuningCurves.vonMises(os, best), 'LineWidth', 2);
            plot(os, TuningCurves.vonMises(os, worst), 'LineStyle', '--');
            % plot where we thing preferred orientation is (black)
            plot([best(4) best(4)], [0, max(rates)], '--k')
            % neurons already have some estimated tuning from A.B. et al
            % (green)
            et = pops_task(pi).tuning(ni);
            if ~isnan(et)
                plot([et et], [0, max(rates)], '--g')
            end
            axis([0,180,0,max(rates)+10]);
            hold off;
            title(sprintf('Population %d Neuron %d', pi, ni));
            drawnow; pause;
            
            %             amp(j) = best(2) * exp(best(3));
            %             fve(j) = Util.Variance_Explained(rates, TuningCurves.vonMises(orientations, best));
            %             tun(j) = pops_task(pi).tuning(ni);
            %             if isnan(fve(j))
            %                 keyboard
            %             end
            %             if (fve(j) < 0.2 && ~isnan(tun(j))) || isnan(tun(j))
            %                 % look at ones without VE but with tuning well-defined?
            %                 % scatter plot with tuning curve overlayed
            %                 scatter(orientations, rates);
            %                 hold on;
            %                 plot(os, TuningCurves.vonMises(os, best), 'LineWidth', 2);
            %                 plot(os, TuningCurves.vonMises(os, worst), 'LineStyle', '--');
            %                 % plot where we thing preferred orientation is (black)
            %                 plot([best(4) best(4)], [0, max(rates)], '--k')
            %                 % neurons already have some estimated tuning from A.B. et al
            %                 % (green)
            %                 et = pops_task(pi).tuning(ni);
            %                 if ~isnan(et)
            %                     plot([et et], [0, max(rates)], '--g')
            %                 else
            %                     plot([best(4) best(4)], [0, max(rates)], '--r')
            %                 end
            %                 axis([0,180,0,max(rates)+10]);
            %                 hold off;
            %                 title(sprintf('Population %d Neuron %d', pi, ni));
            %                 drawnow; pause;
            %             end
            j = j+1;
        else
            warning('#rates was not same as #orientations');
        end
    end
    fprintf('end of population %d\n', pi);
end