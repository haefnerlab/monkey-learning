function [ fve ] = Variance_Explained( data_y, predict_y )
%VARIANCE_EXPLAINED calculates fraction variance explained as the fraction
%of variance in the model over variance in the data
%
%The fraction of variance unexplained (FVU) is 1-(fve)

fve = nanvar(predict_y) / nanvar(data_y);

end

