function [val] = vonMisesDeriv(theta, varargin)
% return derivative of vonMises function (with respect to theta),
% where theta is in degrees
%
% inputs identical to vonMises() function

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

theta = theta*pi/180;
theta_pref = theta_pref*pi/180;

% note that I used the identity sin(x)cos(x) = 1/2 sin(2x)
val = - pi * r_max * k * exp(k*cos(theta - theta_pref).^2).*sin(2*(theta-theta_pref)) / 180;

end