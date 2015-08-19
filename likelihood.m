% return the likelihood of a fit function and its associated values\
% input:
% - cr: choice ratio
% - cp: choice probability
% - n:
% returns:
% - llh_eqn1: log likelihood fit for equation 1
% - llh_cst: log likelihood fit for constant model
% - theta_FF: parameter estimated for Feed Forward Model(equation 1)
% - theta_FB: parameter estimated for Feed Back Model(constant)
% - CP_model_eqn1: modelled CP for equation 1 
% - CP_model_const: modelled CP for constant 


function [llh_eqn1,llh_cst, theta_FF,theta_FB,CP_model_eqn1,CP_model_const]...
  = likelihood(cr, cp,n)
      %  the exponential part of equation 1, returns an array with length k
      EXP = exp((norminv(cr,0,1)).^2.*(-1/2))./(4.*cr.*(1-cr));% with length of contrast level
      
       %% calculate maximum likelihood and theta_FF for equation 1 (FF) model
      theta_FF = fminsearch(@(x) fun_eq1('FF',x,EXP,cp,n),0.5, optimset('MaxFunEvals', 4000, 'MaxIter', 2000)); %CP(0.5) parameter correspond to maximum likelihood
      llh_eqn1 = fun_eq1('FF',theta_FF,EXP,cp,n);
      % test(theta_FF, EXP);
      %alternative:
      %theta_star =( sum(CP(1,i,:)*n)*EXP + 1/2*n - 1/2*n*EXP)./(2*EXP*sum(CP(1,i,:).*n)-EXP*n);
      
      % constant (FB) model
      theta_FB = fminsearch(@(x) fun_eq1('FB',x,EXP,cp,n), 0.5, optimset('MaxFunEvals', 4000, 'MaxIter', 2000));
      % lof likelihood constant model
      llh_cst = fun_eq1('FB',theta_FB,EXP,cp,n);
     
      % model predicted value
      CP_model_eqn1 = 0.5 + (theta_FF-0.5).*EXP;
      %%CP_eq1_all(k) = CP_model_eqn1;
      CP_model_const = theta_FB;
      
      
%     log likelihood  
%       LR_neuron_eqn1(k) =  -2*llh_cst +2*llh_eqn1;
%       LR_neuron_const(k) =  -2*llh_eqn1 +2*llh_cst;
%       