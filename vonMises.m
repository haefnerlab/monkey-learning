function [val] = vonMises(theta, r_0, r_max, k, theta_pref)
% return value of vonMises function, where theta is in degrees
val = r_0 + r_max .* exp(k .* cos((theta*pi/180 - theta_pref)).^2);

end