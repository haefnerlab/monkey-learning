function [ handle ] = subplotsquare( total, p )
%SUBPLOTSQUARE do p-th subplot of a total of m, where they are arranged as
%squarely as possible, with 1-fewer rows than cols if necessary

m = ceil(sqrt(total));
n = ceil(total / m);

handle = subplot(n, m, p);

end

