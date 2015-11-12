% this script makes visualizations to show that the Bootstrap_* functions
% work as expected

if ~exist('params', 'var'), params = New_Parameters('monkey', 'jbe'); end

[pops_task, pops_fix] = Load_Preprocess(params);

%%% test Bootstrap_SpikeCounts
% bootpop = Bootstrap_SpikeCounts(pops_task);

%%% test Bootstrap_TuningCurves

figure();
os = linspace(0,180,181);

subplot(4,3,1);
hold on;
plot(pops_fix(1).condVec(:,strcmpi(pops_fix(1).condVecLabel, 'orientation')), pops_fix(1).spikeCounts(1,:), 'o');
plot(os, pops_task(1).tuning_vm_curves{1}(os));

subplot(4,3,2);
hold on;
plot(pops_fix(1).condVec(:,strcmpi(pops_fix(1).condVecLabel, 'orientation')), pops_fix(1).spikeCounts(2,:), 'o');
plot(os, pops_task(1).tuning_vm_curves{2}(os));

subplot(4,3,3);
hold on;
plot(pops_fix(1).condVec(:,strcmpi(pops_fix(1).condVecLabel, 'orientation')), pops_fix(1).spikeCounts(3,:), 'o');
plot(os, pops_task(1).tuning_vm_curves{3}(os));

for resample=1:3
    [bpt, bpf] = Bootstrap_TuningCurves(pops_task(1), pops_fix(1));

    subplot(4,3,resample*3+1);
    hold on;
    plot(bpf(1).condVec(:,strcmpi(bpf(1).condVecLabel, 'orientation')), bpf(1).spikeCounts(1,:), 'o');
    plot(os, bpt(1).tuning_vm_curves{1}(os));
    
    subplot(4,3,resample*3+2);
    hold on;
    plot(bpf(1).condVec(:,strcmpi(bpf(1).condVecLabel, 'orientation')), bpf(1).spikeCounts(2,:), 'o');
    plot(os, bpt(1).tuning_vm_curves{2}(os));
    
    subplot(4,3,resample*3+3);
    hold on;
    plot(bpf(1).condVec(:,strcmpi(bpf(1).condVecLabel, 'orientation')), bpf(1).spikeCounts(3,:), 'o');
    plot(os, bpt(1).tuning_vm_curves{3}(os));
end