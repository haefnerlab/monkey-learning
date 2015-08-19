% Created By Shuchen Wu
% 08/2015
% 
% calculate choice probability and choice ratio
% input:
% i - neuron row number
% condVec - vector specifiying specific conditions, fit in generic data
% TPR_allSignal - calculated the index where monkey makes the correct
% choice
% FPR_allSignal - false positive rate
% spikeCounts - matrix that specify neuron's spikecounts across all trials

% output:
% n - number of neurons
function [CP, CR,n,err_up,err_down] = calculateCPCR(l,condVec,TPR_allSignal, FPR_allSignal, i,spikeCounts,Err_on)
spikeCounts = spikeCounts';
if (condVec==-1)% ignore the case condVec = -1
  condVec = -1;
end

if (l == -1)
  TPR = TPR_allSignal;%True positive rate drawn at the end of the trails
  FPR = FPR_allSignal;%false positive rate with sorted signal
else
  SignalIndex = findOrientationIndex(l,condVec);
  %TPR returns array of column number correspond to monkey's choice = 1
  TPR = intersect(TPR_allSignal, SignalIndex);%True positive rate drawn at the end of the trails
  FPR = intersect(FPR_allSignal, SignalIndex);%false positive rate with sorted signal
end


pref=spikeCounts(i,TPR);%spike counts of that particular neuron across all trails
anti=spikeCounts(i,FPR);

if ~isempty(find(spikeCounts(i,:) == -1, 1))
  CP = nan;
  CR = nan;
  n = nan;
  err_up = nan;
  err_down = nan;
  %CP = zeros(1,neurons);%%cp of neurons across time
  %choiceRatio = zeros(1,neurons);
else
  L_pref = length(pref);
  L_anti = length(anti);
  CP = ROC(pref,anti);
  CR = L_pref/(L_anti+L_pref);
  if (l ~= -1)
    n = length(SignalIndex);
  else
    n = length(TPR);
  end
  if strcmp(Err_on,'ON')
       [err_up,err_down] = sampling_Data_Error(pref, anti,1000);
  else 
         err_up = 0;
         err_down = 0;
  end 
end