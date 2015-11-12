function [ h,l ] = scatterr( x, y, errx, erry, s, c )
%SCATTERR scatter plot with error bars
%
% [h,l] = scatterr(X,Y,ErrX,ErrY,S,C) behaves like the built-in scatter
% function scatter(X,Y,S,C) but with lines X+/-ErrX and Y+/-Erry. returns a
% handle to the scatter plot (h) and a list of handles to line objects (l)

tmp_hold = ishold;

N = length(x);
assert(N == length(y), 'X and Y must have same number of elements');

% first draw the lines so they appear behind the points

l = zeros(1, 2*N);
for i=1:length(x)
    xs = [x(i)-errx(i), x(i), x(i)+errx(i)];
    ys = [y(i)-erry(i), y(i), y(i)+erry(i)];
    % draw horizontal error bar
    l(2*i-1) = line([xs(1) xs(3)], [ys(2) ys(2)], 'Color', c(i,:));
    % draw vertical error bar
    l(2*i)   = line([xs(2) xs(2)], [ys(1) ys(3)], 'Color', c(i,:));
end

% make the scatter plot
hold on;
h = scatter(x, y, s, c);

% return hold state to how it was entering this function
if tmp_hold, hold on;
else hold off;
end

end

