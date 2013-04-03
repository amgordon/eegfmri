function ArtRep_outlier_indices_batch(subpar)
%%%subj_numbers should be an integer or an array of integers corresponding
%%%to the subjects whose data you want to process. e.g. subj_numbers=1
%%%(then subj=s101), or if subj_numbers=[1 3 5 7] then the script will use
%%%a loop to calculate the outliers for each subject (s101, s103, s105,
%%%s107), one at a time.


% ---load par params if need be---
if isstruct(subpar) % if it is par_params struct
    par = subpar;
else % assume subject string
    par = par_params(subpar);
end

if ~exist(par.artrepdir,'dir')
    mkdir(par.artrepdir)
end
cd(par.artrepdir);

ArtRep_find_artifact_timepoints(par);  % run the ArtRepair routines to find the outliers

movement=[];
global_signal=[];

clear delta_cell;
clear zscoreA_cel;


eval(['load art_global_modified_' par.substr '.mat'])

num_runs = par.numscans;

%%Calculates the outliers for movement threshold .2 .25 .3 .35 .4 .4 .5 and for global signal threshold 2 2.5 3 3.5 4 4.5 5
for i=1:7
    [zscore_outliers, delta_outliers, all_outliers, zscore_trial_nums, delta_trial_nums] = ArtRep_find_outliers(num_runs, zscoreA_cell, delta_cell, 2+(i-1)*.5,.2+(i-1)*.05);
    movement_outlier_trials{i}= delta_trial_nums;
    global_signal_outlier_trials{i}=zscore_trial_nums;
    all_outlier_timepoints{i} = all_outliers;
    
end
%
%     movementArrName=['movement_' subj ];
%     global_signalArrName=['global_signal' subj];
%
%     eval([movementArrName '= movement']);
%     eval([global_signalArrName '=global_signal']);


save_cmd = ['save ' par.substr '_outlier_indices.mat movement_outlier_trials global_signal_outlier_trials all_outlier_timepoints'];
eval(save_cmd);





