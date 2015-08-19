function p = ROC(x,y)

% function p = ROC(x,y)
%
% Computes "area under the curve" in ROC analyis
% p=1 implies max(x)<min(y)
counter = 0;
nx=length(x);
ny=length(y);

%for i=1:nx
  for j=1:ny
    %if x(i) ==0 || y(j) == 0;
    array = x<y(j);
    counter = counter + sum(array);
    array = x==y(j);
    counter = counter + .5 * sum(array);
    %if x(i)<y(j),      counter=counter+1;
    %elseif x(i)==y(j),     counter=counter+0.5;
    %end
  end
  p=counter/nx/ny;% choice probability


