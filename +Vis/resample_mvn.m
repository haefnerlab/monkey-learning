function [ rows ] = resample_mvn(rows, covariance)
% resample_mvn fills in NaN columns sampled from the covariance sliced by
% the non-nan values in row
% (assuming data centered with mean zero already)

for r=1:size(rows,1)
    row = rows(r,:);
    if ~any(isnan(row)), continue; end

    observed = ~isnan(row);

    cov_11 = covariance(~observed, ~observed);
    cov_12 = covariance(~observed,  observed);
    cov_21 = covariance( observed, ~observed);
    cov_22 = covariance( observed,  observed);

    % see https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
    mu_slice = cov_12 * (cov_22 \ row(observed)');
    cov_slice = cov_11 - cov_12 * (cov_22 \ cov_21);

    % sample from sliced distribution
    rows(r,~observed) = mvnrnd(mu_slice, cov_slice);
    
end

end