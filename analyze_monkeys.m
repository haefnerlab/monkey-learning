% top-level script for monkey data analysis

%% Load and preprocess
populations = Load_Data('lem');
populations = Split_Conditions( populations );

%% Scatter plots of moments and 
figure();
colors = hsv(length(populations));

nmoments = 3;

for moment = 1:nmoments
    subplotsquare(nmoments, moment);
    title(sprintf('Moment %d', moment));
    
    for pi=1:length(populations)
        
    end
end