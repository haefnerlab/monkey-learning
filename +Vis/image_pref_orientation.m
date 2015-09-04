function [image, counts] = image_pref_orientation(orient_x, orient_y, xydata, imsize, circlesize)
%IMAGE_PREF_ORIENTATION 
%
%	I = IMAGE_PREF_ORIENTATION(orient_x, orient_y, xydata, imsize, circlesize) 
%	creates a square image where the axes correspond to 0:180. Each data point 
%	xydata(i) is plotted by adding a circle of diameter circlesize (degrees) to the image centered 
%	at orient_x(i) orient_y(i). The gaussian has width (in degrees) defined by 
%	standard deviation circlesize. The image is imsizeXimsize pixels.
%   Each pixel is averaged based on the number of data points that
%   contributed to it. threshold determines how much of each gaussian is
%   taken. The image is set to NaN where zero data are counted.
%	
%	imsize and circlesize are optional parameters

if nargin < 4, imsize = 180; end
if nargin < 5, circlesize = 10; end

image = zeros(imsize);
counts = zeros(imsize);

x = linspace(0,180,imsize+1); x = x(1:end-1);
y = linspace(0,180,imsize+1); y = y(1:end-1);
[xx, yy] = meshgrid(x,y);

ndata = length(orient_x);
assert(length(orient_y) == ndata);
assert(length(xydata) == ndata);

function circle = gaussian_centered_at(cx, cy)
    % center on cx,cy and wrap around edges (toroidal)
	shifted_xx = mod(xx - cx, 180); shifted_xx = min(shifted_xx, 180-shifted_xx);
	shifted_yy = mod(yy - cy, 180); shifted_yy = min(shifted_yy, 180-shifted_yy);
	circle = exp(-0.5*(shifted_xx.^2 + shifted_yy.^2) / (circlesize^2)) > 0.67;
end

for i = 1:ndata
    circ = gaussian_centered_at(orient_x(i), orient_y(i));
    image = image + xydata(i) * circ;
    counts = counts + circ; % count # values that contribute to each pixel
end

% average the image based on counts, setting zero counts to NaN
image(counts == 0) = NaN;
image = image ./ counts;

end