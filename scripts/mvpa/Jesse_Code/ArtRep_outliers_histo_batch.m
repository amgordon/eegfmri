function ArtRep_outliers_histo_batch(subj_numbers)
%%%subj_numbers should be an integer or an array of integers corresponding
%%%to the subjects whose data you want to process. e.g. subj_numbers=1
%%%(then subj=s101), or if subj_numbers=[1 3 5 7] then the script will use
%%%a loop to calculate the outliers for each subject (s101, s103, s105,
%%%s107), one at a time.


for j=subj_numbers
    subj=['s1' prepend(num2str(j), 2)];


    movement=[];
    global_signal=[];
    
    clear delta_cell;
    clear zscoreA_cel;
    
    
    eval(['load /Users/Jesse/fMRI/matlab/VMB_data/art_global_modified_' subj '.mat'])
    
    
    %%Calculates the outliers for movement threshold .2 .25 .3 .35 .4 .45
    %%.5 and for gloval signal threshold 2 2.5 3 3.5 4 4.5 5
    for i=1:7
            [zscore_outliers, delta_outliers, all_outliers, zscore_trial_nums, delta_trial_nums] = ArtRep_find_outliers(zscoreA_cell, delta_cell, 2+(i-1)*.5,.2+(i-1)*.05);
            movement(i)=length(delta_trial_nums);
            global_signal(i)=length(zscore_trial_nums);
           
    end
    
    movementArrName=['movement_' subj ];
    global_signalArrName=['global_signal' subj];
    
    eval([movementArrName '= movement']);
    eval([global_signalArrName '=global_signal']);
    
    
    save([subj '_outlier_histo.mat'], movementArrName, global_signalArrName);
    
end
    
    
    
    
    
    
    