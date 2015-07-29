function [val] = vonMises(theta, varargin)
% return value of vonMises function, where theta is in degrees
%
% varargin is either [r_0, r_max, k, theta_pref], or that same array
% unpacked into separate args

if length(varargin) == 1
    r_0 = varargin{1}(1);
    r_max = varargin{1}(2);
    k = varargin{1}(3);
    theta_pref = varargin{1}(4);
else % multiple args in
    r_0 = varargin{1};
    r_max = varargin{2};
    k = varargin{3};
    theta_pref = varargin{4};
end

val = r_0 + r_max .* exp(k .* cos((theta*pi/180 - theta_pref)).^2);

val(val < eps) = eps;

end