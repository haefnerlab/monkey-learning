% Created by Shuchen Wu 
% 08/2015
% Main function test the goodness of fit(log likelihood) from the measured
% monkey data

function testMonkeyFit

%%load monkey's data separate
% Variables Explained
% num_trail: array of number of trials used to measure each neuron
[size populations] = Load_Task_Data( 'jbe');% load all data into one population struct
[const_ll_jbe,eqn1_ll_jbe,LR_neuron_eqn1_jbe,LR_neuron_const_jbe,...
llh_cst_all_jbe,llh_eqn1_all_jbe,number_of_trail_jbe,...
least_square_sum_eq1_jbe,least_square_sum_const_jbe,const_CP_jbe]...
= fitMonkeyData(size,populations,'OFF');

[size populations] = Load_Task_Data( 'lem');
[const_ll_lem,eqn1_ll_lem,LR_neuron_eqn1_lem,LR_neuron_const_lem,...
llh_cst_all_lem,llh_eqn1_all_lem,number_of_trail_lem,...
least_square_sum_eq1_lem,least_square_sum_const_lem,const_CP_lem]...
= fitMonkeyData(size,populations,'OFF');

%% scatter plot of 2 monkeys
jbe_constant = const_ll_jbe;
jbe_eqn1 = eqn1_ll_jbe;
lem_constant = const_ll_lem;
lem_eqn1 = eqn1_ll_lem;
figure
scatter(lem_constant,lem_eqn1,number_of_trail_lem/10,'b')
hold on
scatter(jbe_constant,jbe_eqn1,number_of_trail_jbe/10,'r')
xlabel('Negative log likelihood for constant(Feed Back) model');
ylabel('Negative log likelihood for equation 1(Feed Forward) model')
title('Negative log likelihood distribution for 2 monkeys')
legend('Monkey lem','Monkey jbe')




%%combining data from 2 monkeys
const_ll = [const_ll_jbe,const_ll_lem];
eqn1_ll = [eqn1_ll_jbe,eqn1_ll_lem];
LR_neuron_eqn1 = [LR_neuron_eqn1_jbe,LR_neuron_eqn1_lem];
LR_neuron_const = [LR_neuron_const_jbe,LR_neuron_const_lem];
llh_cst_all = [llh_cst_all_jbe,llh_cst_all_lem];
llh_eqn1_all = [llh_eqn1_all_jbe,llh_eqn1_all_lem];
least_square_sum_eq1_all = [least_square_sum_eq1_jbe,least_square_sum_eq1_lem];
least_square_sum_const_all = [least_square_sum_const_jbe,least_square_sum_const_lem];
num_of_trial_all = [number_of_trail_jbe,number_of_trail_lem];
const_CP_all = [const_CP_jbe, const_CP_jbe];
%% histogram of log likelihood for the constant and equation1 model
figure
subplot(2,2,1);
hist(const_ll,20)
legend('Negative LogLikelihood Of Feed Back Model')
xlabel('Negative Loglikelihood Value')
ylabel('Frequency')
title('Histogram of Log Likelihood Distribution for Feed Back Model')
axis([0 1620 0 30]);

subplot(2,2,2);
hist(eqn1_ll,20)
legend('Negative LogLikelihood Of Feed Forward Model')
xlabel('Negative Loglikelihood Value')
ylabel('Frequency')
title('Histogram of Log Likelihood Distribution for Feed Forward Model')
axis([0 1620 0 30]);

%% calculate likelihoood ratio across all neurons
Likelihood_Ratio = -2*llh_cst_all +2*llh_eqn1_all;

%% plot likelihood ratio between neurons without sorting 
subplot(2,2,3);
plot(LR_neuron_eqn1)
xlabel('Neurons in Experimental Order');
ylabel('log likelihood ratio')
title('log likelihood ratio across neurons assuming null hypothesis of Feed Back Model')
axis([0 188 -17 31]);
% subplot(1,2,2);
% plot(LR_neuron_const)
% xlabel('neurons');
% ylabel('likelihood ratio')
% title('likelihood ratio across neurons assuming null hypothesis of Feed Forward Model')
% axis([0 200 -20 20]);
%% sort neurons by negaitve log likelihood and plot
[sortedValues_eqn1,sortIndex_eqn1] = sort(LR_neuron_eqn1(:),'descend');
maxIndex_eqn1 = sortIndex_eqn1(1:5);
[sortedValues_const,sortIndex_const] = sort(const_ll(:),'descend');
maxIndex_eqn1_const = sortIndex_const(1:5);
subplot(2,2,4);

plot(sortedValues_eqn1)
xlabel('neurons with log likelihood ratio sorted in decending order')
ylabel('log likelihood ratio')
axis([0 188 -17 31]);
title('Sorted log likelihood ratio across neurons assuming null hypothesis of Feed Back Model')

% subplot(1,2,2)
% plot(sortedValues_const)
% xlabel('neurons with likelihood ratio in decending order')
% ylabel('likelihood ratio')
% axis([0 200 -20 20]);
% title('likelihood ratio across neurons assuming null hypothesis of Feed Forward Model')


% compute "summary" things
% plot "summary" things