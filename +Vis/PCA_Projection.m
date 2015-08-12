function [ eigenvectors, mean_vec, covariance, fig, handles ] = PCA_Projection( data, classes, colors, sizes, fig )
%PCA_PROJECTION Create a 2D projection of data, colored by class labels.
%  'data' must have variables along columns and observations along rows
%  'classes' are Nx1 labels {1,2,...,C}, which together with 'colors'
%            determines how to color each point
%  'colors' are Cx3 RGB colors to use for each class
%  'sizes' are Cx1 scatter-plot sizes for each class
% 
% Creates a new figure if one not given

if nargin < 2
    classes = ones(size(data,1),1);
end

if nargin < 3
    colors = hsv(length(unique(classes)));
end

if nargin < 4
    sizes = 20*ones(length(unique(classes)), 1);
end
if nargin < 5
    fig = figure();
else
    figure(fig);
end

mean_vec = nanmean(data,1);

data_centered = data - repmat(mean_vec, size(data,1), 1);

% Here, we do a version of PCA that can handle NaN measurements. PCA simply
% finds the top-K eigenvectors of the covariance matrix. So we begin with a
% NaN-safe computation of covariance
covariance = nancov(data_centered, 'pairwise');

% NOTE that in rare cases there may be two neurons which share zero non-nan
%  trials. We set those covariances to 0 here, which is justified by the
%  fact that we use the covariance for predictions to fill in missing
%  values, and zero covariance corresponds to nil predictive power, which
%  is essentially what we have.
covariance(isnan(covariance)) = 0;

% Principle components (here we want 2 of them to project into 2D) are
% largest eigenvectors of the covariance matrix
[eigenvectors,~] = eigs(covariance, 2);

% for the sake of consistency across runs, our standard will be a positive
% first component of each
if eigenvectors(1,1) < 0
    eigenvectors(:,1) = eigenvectors(:,1) * -1;
end
if eigenvectors(1,2) < 0
    eigenvectors(:,2) = eigenvectors(:,2) * -1;
end

% Before projecting, we sample from sliced-covariances to fill in 'missing'
% (i.e. NaN) data.
% NOTE that by using covariance, we only preserve mean and covariance of
% the samples; no higher-order statistics are employed.
data_augmented = Vis.resample_mvn(data_centered, covariance);

% project data: NxP * Px2 --> Nx2
data_projected = data_augmented * eigenvectors;

hold on;
for c=size(colors,1):-1:1
    class_data = data_projected(classes == c, :);
    handles(c) = scatter(class_data(:,1), class_data(:,2), sizes(c), colors(c,:), 'filled');
end
hold off;

end