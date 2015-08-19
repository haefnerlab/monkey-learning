%created by Shuchen Wu July 23/2015
% Created by Shuchen Wu 
% 08/2015


%returns the likehood of constant model and the model by arxiv paper and computes the likelihood ratio
%Equation Declaration for maximum likelihood fitting
%Likelihood_ratio = 2 log(P(y|theta)) - 2 log(P(y|H0))
%p = x/n; %parameter to fit the model x = sum(CPi*n)
% calculation of errorbar for math models: 
%sqrt(p*(n-p)/n)

%%the parameter fit by paper(the best choice ratio refer to maximized probability)
%%theta_star = sum(CPi*ni - ni)/(sum(ni).*exp(-1/2.*(norminv(CR,0,1)))^2)
% lgh_constant = @(theta) (-1).*(log(theta)*sum(CP(1,i,:)*n) + log(1-theta)*(sum(n-CP(1,i,:).*n)));
% theta = fminsearch(lgh_constant,0.1);
% loglikelihood = log(theta)*sum(CP(1,i,:)*n) + log(1-theta)*(sum(n-CP(1,i,:).*n));
      

% 
% inputs:
% - Size: the number of structs in the data
% - populations: a single struct containing all the information
% - displayPlot: 'ON' or 'OFF' depending on whether plot display is wanted
% 
% return values: 
% - const_likelihood: negative log likelihood for each neuron correspond to constant model
% - eqn1_likelihood: negative log likelihood correspond to equation 1 model
% - LR_neuron_eqn1: likelihood ratio for each neuron assuming null hypothesis of constant
% model
% - LR_neuron_const:likelihood ratio for each neuron assuming null
% hypothesis of equation 1 model
% - llh_cst_all: cumulative negative log likelihood for neurons fit by
% constant model
% - llh_eqn1_all:cumulative negative log likelihood for neurons fit by
% equation 1 model
% - number_of_trail: correspond to measurement of that particular neuron
% - least_square_sum: another method to measure goodness of fit


function [const_likelihood,eqn1_likelihood,LR_neuron_eqn1,LR_neuron_const,...
  llh_cst_all,llh_eqn1_all,number_of_trail,least_square_sum_eq1, ...
  least_square_sum_const,const_CP] = fitMonkeyData(Size,populations, displayPlot)
%% initialize variables
llh_eqn1_all = 0;
llh_cst_all = 0;
n_neuron = 0;

% LR_neuron_eqn1 = zeros(1,188);%likelihood ratio for one neuron
% LR_neuron_const = zeros(1,188);%likelihood ratio for one neuron
% eqn1_likelihood = zeros(1,188);
% const_likelihood = zeros(1,188);
% best_eqn1 = [153,173,152,176,187];
% best_const = [169,184,170,154,185];
% best_eqn1 = [29,26,30,25,28];
% best_const = [29,26,30,28,25];
best_lsf_eqn1 = [145,34,105,29,62];
best_lsf_const = [145,34,105,29,62];
best_eqn1 = [6,26,5,29,40];
best_const = [22,37,23,7,38];
% worst_const = [176,173,177,175,172,164,178,174,163,159];
worst_const = [91,99,97,94,103,93,101,92,96,102];
b_lsf_c = 1;
b_lsf_e = 1;
w_c = 1;
b_e = 1;
b_c = 1;

% best by least square for equation 1
% [145,34,105,29,62,18,60,20,91,13]
b_lsf_eq1 = [145,34,105,29,62,18,60,20,91,13];
% worst by least square for equation 1 
%  [172   173   170   187   178   176   184   174   160   162]
% correspond to
% [  25    26    23    40    31    29    37    27    13    15] in lem
w_lsf_eq1 = [  25    26    23    40    31    29    37    27    13    15];

% best by least square for constant
% [145,34,105,29,62,18,60,104,20,91]
b_lsf_const = [145,34,105,29,62,18,60,104,20,91];
% worst by least square for constant
% [172   170   184   173   178   174   162   160   187   176]
% correspond to 
% [25    23    37    26    31    27    15    13    40    29]
% in lem
w_lsf_const = [25    23    37    26    31    27    15    13    40    29];
% most trials with the most number of n and thus smaller error bar
most_trails = [90    91    92    93    94    95    96    97    98    99   ...
  100   101   102   103  104   105    24    25    26    27    28    29    ...
  30    31    32    33    34    35   36    37];
largest_CP = [17    29    26     6    16     7     5    37    12    41];%%lem


% %initializing value for scatter plot separating 2 monkeys
% %jbe has 147 neurons being measured
% jbe_constant = zeros(1,147);
% jbe_eqn1 = zeros(1,147);
% jbe_size = zeros(1,147);
% %jbe has 147 neurons being measured
% lem_constant = zeros(1,41);
% lem_eqn1 = zeros(1,41);
% lem_size = zeros(1,41);

%% loop through structs
for j = 1:Size
  %% declare variables
 RealChoice = populations(j).realChoice;
  condVec = populations(j).condVec;
  TPR_allSignal = find(RealChoice == -1);%index where monkey chooses 1 with unsorted signal strength
  FPR_allSignal = find(RealChoice ==1);%index where monkey chooses -1
  spikeCounts = populations(j).spikeCounts;
  spikeCounts(isnan(spikeCounts)) = -1;
  
  [neurons,trails] = size(spikeCounts);
  
  CP = zeros(1,neurons);%%cp of neurons across time
  err_up = zeros(1,neurons);% upward error for original data
  err_down = zeros(1,neurons);% downward error for original data
  choiceRatio = zeros(1,neurons);
  %minlength = min(length(pref),length(anti))
  
  %within each struct, find the CP according to different signal strength
  N = length(spikeCounts(1,:));
  n = zeros(1,7);
  num_trial = 0;
 %% collect CP data across different stimulus
  for i = 1: neurons %collect CP and choiceRatio for all neurons
    containsNAN = false;%initialize containsNAN   
    % calculate CP and CR
    for l = 1:7
       [CP(1,i,l), choiceRatio(1,i,l),n(1, l),err_up(1,i,l),err_down(1,i,l)] = calculateCPCR(l,condVec,TPR_allSignal, FPR_allSignal, i,spikeCounts);
    end
    
    %theta= sum(CP(1,i,:))*n/(n*7);
    cr = squeeze(choiceRatio(1,i,:));
    cp = squeeze(CP(1,i,:));
    for w=1:7%%check for nans in data
      if isnan(cr(w))
        containsNAN = true;
      end
    end
    
    if (containsNAN == true)%choosing the neurons without nans in measured values in any of the 7 stimulus
    else
      %count the number of neurons  and record number of trails measured
      %for that neuron
      n_neuron = n_neuron +1;
      number_of_trail(1,n_neuron) = length(spikeCounts);
      
      %% calculate log likelihood for each model 
      
      %  the exponential part of equation 1, returns an array with length n_neuron     
      EXP = exp((norminv(cr,0,1)).^2.*(-1/2))./(4.*cr.*(1-cr));
      
      % calculate maximum likelihood and theta_n for equation 1 (FF) model
      theta_n = fminsearch(@(x) fun_eq1('FF',x,EXP,cp,n),0.5, optimset('MaxFunEvals', 4000, 'MaxIter', 2000)); %CP(0.5) parameter correspond to maximum likelihood
      llh_eqn1 = fun_eq1('FF',theta_n,EXP,cp,n);
      llh_eqn1_all = llh_eqn1_all+llh_eqn1;
%       test(theta_n, EXP);
      %alternative:
      %theta_star =( sum(CP(1,i,:)*n)*EXP + 1/2*n - 1/2*n*EXP)./(2*EXP*sum(CP(1,i,:).*n)-EXP*n);
      
      
      if isnan(llh_eqn1)
        disp(i);
      end%in case something is wrong
            
      % constant (FB) model
      theta = fminsearch(@(x) fun_eq1('FB',x,EXP,cp,n), 0.5, optimset('MaxFunEvals', 4000, 'MaxIter', 2000));
      const_CP(1,n_neuron) = theta;
      % lof likelihood constant model
      llh_cst = fun_eq1('FB',theta,EXP,cp,n);
      llh_cst_all = llh_cst_all+llh_cst;
      
      LR_neuron_eqn1(n_neuron) =  -2*llh_cst +2*llh_eqn1;
      LR_neuron_const(n_neuron) =  -2*llh_eqn1 +2*llh_cst;
      
      %record likelihood for histogram plotting
      eqn1_likelihood(n_neuron) =  llh_eqn1;
      const_likelihood(n_neuron) =  llh_cst;
      
      % model predicted value
      CP_model_eqn1 = 0.5 + (theta_n-0.5).*EXP;
      %%CP_eq1_all(n_neuron) = CP_model_eqn1;
      CP_model_const(1:7) = theta;
%       MockData(cr,CP_model_const, EXP, n);
      % Least Square Sum
      least_square_sum_eq1(n_neuron) = sum((cp - CP_model_eqn1 ).^2);
      least_square_sum_const(n_neuron) = sum((cp - (CP_model_const') ).^2);
      
      %error bar, should it be not symmetric along the center?? CHECK
      %LATER!!
      error_bar = sqrt((cp).*(1-cp)./(n'));%calculate error bar
      x = linspace(0,1,100);
      EXP_smooth = exp((norminv(x,0,1)).^2.*(-1/2))./(4.*x.*(1-x));
      CP_model_eqn1_smooth = 0.5 + (theta_n-0.5).*EXP_smooth;
      mean_CP(1,n_neuron) = mean(cp);
      %% plot largest CPs
      if ((~isempty(find(largest_CP==n_neuron, 1)))&& (Size == 6))%%Size == 22
        label = rem(w_c,4);
        if (label == 0)
          label = 4;
        end
        if (label==1)
        figure% create new plot every 4 sunplots
        end
        subplot(2,2,label)
        Likelihood_Ratio = -2*llh_cst +2*llh_eqn1;

        h1 = errorbar(cr,cp,cp-squeeze(err_down(1,i,:)),squeeze(err_up(1,i,:))-cp,'b*');
        hold on
        plot(x,CP_model_eqn1_smooth,'r-');
        hold on;
        h = errorbar(cr,CP_model_eqn1,error_bar,'r-');
        set(h,'Linestyle','none')
        hold on
        plot(cr,CP_model_const,'m-')
        errorbar(cr,CP_model_const,error_bar,'m-')
        legend('measured CP',...
           'Eqn1 predicted CP','Constant Model predicted CP');
%         legend('Measured Choice Probability','Equation 1 Model Predicted Choice Probability')
        xlabel('CR')
        ylabel('CP')
%         str = {'neuron with the most number of trials',...
%           'Negative Log Likelihood of FB(Constant fit)', llh_cst};
        str = {'neuron with biggest deviation of CP', 'Likelihood Ratio: ',...
          Likelihood_Ratio,'constant model parameter: ', theta, ...
          'equation 1 model parameter',theta_n};
        title(str)
        axis([0 1 0 1]);
        w_c = w_c + 1;
      end
      %% plot the 5 best least square fit  both model
      if ((~isempty(find(w_lsf_const==n_neuron, 1)))&& (Size == 29))
        
        if (rem(w_c,4)==1)
          figure% create new plot every 4 sunplots
        end
        subplot(1,4,rem(w_c, 4))
        h2 = errorbar(cr,cp,cp-squeeze(err_down(1,i,:)),squeeze(err_up(1,i,:))-cp,'b*');
        
        hold on
        plot(cr,CP_model_const,'r')
        errorbar(cr,CP_model_const,error_bar,'r')
        legend('Measured Choice Probability','Constant Model Predicted Choice Probability')
        
        xlabel('CR')
        ylabel('CP')
        str = {'neuron with the biggest least square sum by constant (FB) model',...
          'Least Square Sum of FB(Constant fit)', least_square_sum_const(n_neuron)};
        title(str)
        axis([0 1 0.1 1]);
        w_c = w_c + 1;
      
      end 
      %% plot characteristic deviation from the constant model
      if ((~isempty(find(most_trails==n_neuron, 1)))&& (Size ~= 22))%%
        label = rem(w_c,4);
        if (label == 0)
          label = 4;
        end
        if (label==1)
        figure% create new plot every 4 sunplots
        end
        
        subplot(2,2,label)
        h1 = errorbar(cr,cp,cp-squeeze(err_down(1,i,:)),squeeze(err_up(1,i,:))-cp,'b*');
        hold on
        plot(x,CP_model_eqn1_smooth,'r');

%         plot(cr,CP_model_const,'r')
%         errorbar(cr,CP_model_const,error_bar,'r')
                   hold on;
        h = errorbar(cr,CP_model_eqn1,error_bar,'r');
        set(h,'Linestyle','none')

%         legend('Measured Choice Probability',...
%           'Constant Model Predicted Choice Probability')
        legend('Measured Choice Probability','Equation 1 Model Predicted Choice Probability')
        xlabel('CR')
        ylabel('CP')
%         str = {'neuron with the most number of trials',...
%           'Negative Log Likelihood of FB(Constant fit)', llh_cst};
        str = {'neuron with the most number of trials',...
          'Negative Log Likelihood of FF(Equation 1 fit)', llh_eqn1};
        title(str)
        axis([0 1 0.1 1]);
        w_c = w_c + 1;
      end
       %% If specified, display the best 5 fits for FF and FB model
       if (strcmp(displayPlot,'ON') && (Size == 22))
         
         % plot the best modeling cases for eqn1
         if (~isempty(find(best_eqn1==n_neuron, 1)))
           
           subplot(2,5,b_e)

           plot(x,CP_model_eqn1_smooth,'r');
           hold on 
           h1 = errorbar(cr,cp,cp-squeeze(err_down(1,i,:)),squeeze(err_up(1,i,:))-cp,'b*');
           hold on;
           h = errorbar(cr,CP_model_eqn1,error_bar,'r');

           set(h,'Linestyle','none')

           legend('Measured Choice Probability','Equation 1 Model Predicted Choice Probability')
           xlabel('CR')
           ylabel('CP')
           str = {'top neuron with the best likelihood ratio fit by equation 1 (FF) model',...
             'Negative Log Likelihood of FF(Equation 1 fit)',llh_eqn1,...
             'Negative Log Likelihood of FB(Constant fit)', llh_cst};
           title(str)
           axis([0 1 0.1 1]);
           b_e = b_e + 1;
         end    
         
         % plot the worst modeling cases for eqn1(best for constant model)
         if (~isempty(find(best_const==n_neuron, 1)))
           subplot(2,5, 5+ b_c);
           %plot(cr,,'b*')
           h1 = errorbar(cr,cp,cp-squeeze(err_down(1,i,:)),squeeze(err_up(1,i,:))-cp,'b*');

           hold on
           plot(cr,CP_model_const,'r')
           errorbar(cr,CP_model_const,error_bar,'r')
           legend('Measured Choice Probability','Constant Model Predicted Choice Probability')
           
           xlabel('CR')
           ylabel('CP')
           str = {'top neuron with the best likelihood ratio fit by constant (FB) model',...
             'Negative Log Likelihood of FF(Equation 1 fit)',llh_eqn1,...
             'Negative Log Likelihood of FB(Constant fit)', llh_cst};
           title(str)
           axis([0 1 0.1 1]);
           b_c = b_c + 1;
         end         
       end 
    end    
  end
end
% figure
% plot(number_of_trail,CP_eq1_all,'*');
% xlabel('number of trials')
% ylabel('predicted choice probability')





      
%%        Alternative Method
%       if (llh_eqn1==Inf)
%                 disp(i);
% 
%         min = 1000000;
%         min_theta = 0;
%             for q=1:100
%              theta_n_plot = q/100;
%              % logl = sum(log((theta_n_plot-0.5).*(EXP(:))+0.5).*(cp(:)*n))+sum(log(0.5-(theta_n_plot-0.5).*EXP(:)).*(n-cp(:).*n));
%              logl = fun_eq1(theta_n_plot,EXP,cp,n);
%              plot(q,logl);
%              hold on;
%              if (min > logl )
%                min = logl;
%                min_theta = theta_n_plot;
%              end 
%             end
%           theta_n =   min_theta;
%           llh_eqn1 =   min;        
%       end 
