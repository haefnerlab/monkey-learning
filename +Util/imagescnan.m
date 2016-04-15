function [ handle ] = imagescnan( imdata, nancolor )
%IMAGESCNAN imagesc, but with NaN data white

if any(isnan(imdata))
    if nargin < 2
        nancolor = [1,1,1];
    end
    
    % thanks to user 'pipo' on this thread: http://www.mathworks.com/matlabcentral/newsreader/view_thread/140607
    maxval = max(imdata(:));
    imdata(isnan(imdata)) = maxval + maxval/10;
    handle = imagesc(imdata);
    
    colordata = colormap;
    colordata(end,:) = nancolor;
    colormap(colordata);
else
    handle = imagesc(imdata);
end

end

