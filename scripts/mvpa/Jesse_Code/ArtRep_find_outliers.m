function [zscore_outliers, delta_outliers, all_outliers, zscore_trial_nums, delta_trial_nums] = ArtRep_find_outliers(num_runs, zscore_cell, delta_cell, zscore_thresh, delta_thresh)


%These store the actual timepoints exceeding the specified threshold
zscore_outliers=[];
delta_outliers=[];

%Store the trial number containing the outlier timepoints
zscore_trial_nums=[];
delta_trial_nums=[];



num_trials_per_run=40;
timepoints_per_run=5;

%Loops through the "ith" run
for i=1:num_runs
    %Loops through the "jth" trial in the "ith" run
    for j=1:num_trials_per_run
        
        %Loops through the five timepoints in the jth trial in the ith run
        for k=1:timepoints_per_run
            
            
            %%If the "current" timepoint exceeds the threshold, add the
            %%timepoint to zscore_outliers and add the corresponding trial
            %%number to zscore_trial_nums
            if zscore_cell{i}(1+5*(j-1)+k)>=zscore_thresh || zscore_cell{i}(1+5*(j-1)+k)<=(-1)*zscore_thresh
                %Adds the flagged timepoint (the number corresponding to
                %the actual timepoint) to the array zscore_outliers.
                zscore_outliers(length(zscore_outliers)+1)=203*(i-1)+(j-1)*5+k; %#ok<AGROW>
                %Adds the trial number to the array zscore_outlier_runs;
                zscore_trial_nums(length(zscore_trial_nums)+1)=(i-1)*num_trials_per_run+j; %#ok<AGROW>
            end
            
            if delta_cell{i}(1+5*(j-1)+k)>=delta_thresh || delta_cell{i}(1+5*(j-1)+k)<=(-1)*delta_thresh
                delta_outliers(length(delta_outliers)+1)=203*(i-1)+(j-1)*5+k; %#ok<AGROW>
                delta_trial_nums(length(delta_trial_nums)+1)=(i-1)*num_trials_per_run+j; %#ok<AGROW>
            end
        end
    end
end

%Remove duplicates of the same trial number
zscore_trial_nums=unique(zscore_trial_nums);
delta_trial_nums=unique(delta_trial_nums);



all_outliers=union(zscore_outliers, delta_outliers);


 