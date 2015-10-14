deltas = [0, 0.005, 0.02, 0.08];

cs = cell(4,3);
ps = cell(4,3);

for di=1:length(deltas)
    sim = sprintf('SIM_d%.3f_v1', deltas(di));
    [c,p,o] = analyze_fprime_tuning_curves(New_Parameters('monkey', sim, 'bootstrap', 500), fullfile('data', sim, 'analyze_fprime_results.mat'));
    [cs{di,:}] = Util.meanci(c);
    [ps{di,:}] = Util.meanci(p);
end

figure();
hold on;

colors = 'cgrb';

for di=1:length(deltas)
    cm = cs{di,1};
    cminus = cm - cs{di,2};
    cplus = cs{di,3} - cm;
    
    Vis.boundedline(o, cm, [cminus;cplus]', colors(di), 'alpha');
end
