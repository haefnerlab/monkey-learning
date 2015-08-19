% Created by Shuchen Wu
% 08/2015
%% sample distribution of choice probability from randam draw with replacement from the data
% draws n samples and form a distribution of CP assuming a gaussian
% distribution
% n is number of trails used to measure the spike rates
function [err_up, err_down] = sampling_Data_Error(pref, anti,n)
CP_draw = zeros(1,n);
for i = 1:n
  draw_pref = datasample(pref, length(pref));
  draw_anti = datasample(anti, length(anti));
  CP_draw(1,i) = ROC(draw_pref,draw_anti);
end 
err_up = prctile(CP_draw,16);
err_down = prctile(CP_draw,84);



