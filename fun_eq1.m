% Created by Shuchen Wu
% 08/2015
% 
%% computes the negative log likelihood specified by the model
% input:
% - n: vector of numbers specifying the number of trials measured for each
% stimuli
% - theta: fminsearched optimized model parameter
% - model: string constant specifying 'FF' or 'FB'
% log likelihood =  xlog(p)+(n-x)log(1-p)
% x = cp * n
%
% return:
% - y: negative log likelihood
% k = (log(factorial(round(n(i)))) - log(factorial(round(cp(i)*n(i))))-log(factorial(round(n(i) - cp(i)*n(i)))))
%lgh = @(theta) (-1).*(log(theta)*sum(CP(1,i,:)*n) + log(1-theta)*(sum(n-CP(1,i,:).*n)));
% theta = fminsearch(lgh,0); %parameter correspond to maximum likelihood
function y = fun_eq1(model,theta,EXP,cp,n)
y = 0;

switch model
  case 'FB'   %Constant model
    for i = 1:length(cp)
      y =  y + (-1).*(log(theta).*(cp(i).*n(i)) + log(1-theta).*((n(i)-cp(i).*n(i))));
    end 
    
  case 'FF'   %Equation 1 model
    for i = 1:length(cp)
      y = y + (-1).*(log((theta-0.5).*(EXP(i))+0.5).*(cp(i).*n(i))+...
        (log(0.5-(theta-0.5).*EXP(i)).*(n(i)-cp(i).*n(i))));
      if (~isreal(y))
        y = inf;
      end
    end
end

    
      