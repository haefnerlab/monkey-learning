
figure();
for p_idx = 1:length(pops_task)
    Util.subplotsquare(length(pops_task), p_idx);
    
    [orientations,idxs] = sort(pops_task(p_idx).tuning_vm);
    
    plot(orientations, pops_task(p_idx).fprime_stimulus_mean(idxs), 'b-o');
    hold on;
    plot(orientations, pops_task(p_idx).fprime_fixation_means(idxs), 'g-+');
    legend('f''=r(A)-r(B)', 'f''=pwl(A)-pwl(B)');
    task = pops_task(p_idx).Orientation;
    plot([task,task],[-30,30], '--r');
    task90 = mod(task + 90, 180);
    plot([task90,task90],[-30,30], '--r');
    hold off;
    xlim([0,180]);
    title(['population' num2str(p_idx)]);
end