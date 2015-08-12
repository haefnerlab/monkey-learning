function [ mov ] = NaN_resample_movie( pop, nsample )
%NAN_RESAMPLE_MOVIE repeatedly call PCA_2neuron to visualize error
% introduced by 'filling' NaN values with sampled values

mov = struct('cdata', [], 'colormap', []);

f = figure();

for frame=1:nsample
    Vis.PCA_2neuron( pop, f );
    if frame == 1
        ax = axis();
    else
        axis(ax);
    end
    mov(frame) = getframe;
end

end

