function [image, norm] = image_pref_orientation(orient_x, orient_y, xydata, varargin)
%IMAGE_PREF_ORIENTATION 
%
%	I = IMAGE_PREF_ORIENTATION(orient_x, orient_y, xydata, ...) 
%	creates a square image where the axes correspond to 0:180. Each data point 
%	xydata(i) is plotted with some width to "blend" into areas with no data.
%   This can be done either with discs (constant contribution to all points
%   within a radius) or with circular hills. The image is imsizeXimsize pixels.
%   Each pixel is a weighted average from discs/hills. 
%   The image is NaN where zero data are counted.
%	
%   I = IMAGE_PREF_ORIENTATION(..., 'imsize', w) makes a w x w image
%	
%   I = IMAGE_PREF_ORIENTATION(..., 'method', m) where m is 'disc' or
%   'smooth' selects which method is used to fill in missing data
%	
%   I = IMAGE_PREF_ORIENTATION(..., 'radius', r) sets the 'width' of the
%   disc or smooth-hill kernels. It is the radius of the disc or stdev of a
%   gaussian kernel
%	
%   I = IMAGE_PREF_ORIENTATION(..., 'diagonal', d) a boolean flag, defaults
%   to false. Determins whether to ignore data where orient_x == orient_y

imsize = 180;
method = 'disc';
radius = 10;
diagonal = false;

nvararg = length(varargin);
for i=1:2:nvararg
    if strcmpi(varargin{i}, 'imsize'),     imsize = varargin{i+1};
    elseif strcmpi(varargin{i}, 'method'), method = varargin{i+1};
    elseif strcmpi(varargin{i}, 'radius'), radius = varargin{i+1};
    elseif strcmpi(varargin{i}, 'diagonal'), diagonal = varargin{i+1};
    end
end

image = zeros(imsize);
norm = zeros(imsize); % sum of kernels; denominator for averaging

x = linspace(0,180,imsize+1); x = x(1:end-1);
y = linspace(0,180,imsize+1); y = y(1:end-1);
[xx, yy] = meshgrid(x,y);

ndata = length(orient_x);
assert(length(orient_y) == ndata);
assert(length(xydata) == ndata);

function kernel = kernel_centered_at(cx, cy)
    % center on cx,cy and wrap around edges (toroidal)
	shifted_xx = mod(xx - cx, 180); shifted_xx = min(shifted_xx, 180-shifted_xx);
	shifted_yy = mod(yy - cy, 180); shifted_yy = min(shifted_yy, 180-shifted_yy);
	kernel = exp(-0.5*(shifted_xx.^2 + shifted_yy.^2) / (radius^2));
    if strcmp(method, 'disc')
        % disc is just binarized gaussian at 1 std deviation
        kernel = kernel > 0.67;
    end
end

for i = 1:ndata
    if diagonal || (orient_x(i) ~= orient_y(i))
        if ~isnan(xydata(i))
            kernel = kernel_centered_at(orient_x(i), orient_y(i));
            image = image + xydata(i) * kernel;
            norm = norm + kernel; % count # values that contribute to each pixel
        end
    end
end

% average the image based on counts, setting zero counts to NaN
image(norm == 0) = NaN;
image = image ./ norm;

end