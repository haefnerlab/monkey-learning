%% PART I:  load data
% FF = load('Top-Downs0.0.mat');% top down
% FFFB = load('Top-Downs0.5.mat');% top down + bottom up
load('zwischenstand_S_Run_Experiment.mat');
O = cat(4,e{1,3}.O,e{1,4}.O,e{1,5}.O,e{1,6}.O,...
  e{1,7}.O,e{1,8}.O,e{1,9}.O,e{1,10}.O,e{1,11}.O);
X = cat(4,e{1,3}.X,e{1,4}.X,e{1,5}.X,e{1,6}.X,...
  e{1,7}.X,e{1,8}.X,e{1,9}.X,e{1,10}.X,e{1,11}.X);

Size= length(e);
%% PART II: fit samples

IndexMatrix = [
         508         981
         383        1002
         545         958
         492          18
         356         173
         442         968
         496         231
         675          51
         601         932
         387         919];
[llh_eq1_neurons,const_eq1_neurons,theta_FF_neurons,...
          theta_FB_neurons,CP_model_eqn1,CP_model_const,CP, choiceRatio,...
          error_up,error_down,error_bar] =  fitSampling(Size,O,X,IndexMatrix);% input structs with all the data variables
        
        
        
        % fitSampling(size, e)

           avg_CP = mean(CP,2);
%% PART III: Analysis
% b_c = 1;
%            subplot(2,2,b_c);
           %plot(cr,,'b*')
           for i = 1:10
           figure
           k = IndexMatrix(i,1);
           cr = squeeze(choiceRatio(k,:));
           cp = squeeze(CP(k,:));
           h1 =  errorbar(squeeze(choiceRatio(k,:)),squeeze(CP(k,:)),squeeze(CP(k,:)) - squeeze(error_down(k,:,1)),squeeze(error_up(k,:,1))-squeeze(CP(k,:)),'b*');
           hold on
           h2 = errorbar(squeeze(choiceRatio(k,:)),squeeze(CP_model_const(k,:)),squeeze(error_bar(k,:)),'r');
           hold on
           h3 = errorbar(squeeze(choiceRatio(k,:)),squeeze(CP_model_eqn1(k,:)),squeeze(error_bar(k,:)),'m');
           hold on
%            plot(cr,CP_model_const,'r')
%            errorbar(cr,CP_model_const,error_bar,'r')
           legend('Sampling Choice Probability','Constant Model Predicted Choice Probability','Equation1 Model Predicted Choice Proability');
           
           xlabel('CR')
           ylabel('CP')
           str = {'neurons with the ',num2str(i),'th biggest choice proability',...
             'Negative Log Likelihood of FF(Equation 1 fit)',llh_eq1_neurons(k),...
             'Negative Log Likelihood of FB(Constant fit)', const_eq1_neurons(k),...
             'mean choice probability',avg_CP(k)};
           title(str)
           axis([0 1 0.1 1]);
           
           
           figure
            k = IndexMatrix(i,2);
           cr = squeeze(choiceRatio(k,:));
           cp = squeeze(CP(k,:));
           h1 =  errorbar(squeeze(choiceRatio(k,:)),squeeze(CP(k,:)),squeeze(CP(k,:))-squeeze(error_down(k,:,2)),squeeze(error_up(k,:,2))-squeeze(CP(k,:)),'b*');
           hold on
           h2 = errorbar(squeeze(choiceRatio(k,:)),squeeze(CP_model_const(k,:)),squeeze(error_bar(k,:)),'r');
           hold on
           h3 = errorbar(squeeze(choiceRatio(k,:)),squeeze(CP_model_eqn1(k,:)),squeeze(error_bar(k,:)),'m');
           hold on
%            plot(cr,CP_model_const,'r')
%            errorbar(cr,CP_model_const,error_bar,'r')
           legend('Sampling Choice Probability','Constant Model Predicted Choice Probability','Equation1 Model Predicted Choice Proability');
           
           xlabel('CR')
           ylabel('CP')
           str = {'neurons with the ',num2str(i),'th smallest choice proability',...
             'Negative Log Likelihood of FF(Equation 1 fit)',llh_eq1_neurons(k),...
             'Negative Log Likelihood of FB(Constant fit)', const_eq1_neurons(k),...
             'mean choice probability',avg_CP(k)};
           title(str)
           axis([0 1 0.1 1]);
           end 
                      
           [biggestCP,biggestCP_index] = sort(avg_CP(:),'descend');
           biggest_10_CP = biggestCP_index(1:10)';
           
           [smallestCP,smallestCP_index] = sort(avg_CP(:),'descend');
           smallest_10_CP = smallestCP_index(1:10)';
           
           % find cp with important index
           % generate error bar
           CP(smallest_10_CP)
%            b_c = b_c + 1;


