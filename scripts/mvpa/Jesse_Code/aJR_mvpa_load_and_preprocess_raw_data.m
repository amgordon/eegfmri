function [subj] = aJR_mvpa_load_and_preprocess_raw_data(subj_id, exp_name, roi_name, roi_file, raw_filenames, num_runs, num_TP_per_run)

    % initialize subj structure
    subj = init_subj(exp_name,subj_id);

    % load mask file
    subj = load_spm_mask(subj,roi_name,roi_file);
    
    % load functional data
    subj = load_spm_pattern(subj,'spiral',roi_name,raw_filenames,'single',true);

    % move pattern to hard drive
    % subj = move_pattern_to_hd(subj, 'spiral');

    % make runs vector
    subj = init_object(subj,'selector','runs');

    trial_idx = 1;                  
    for r = 1:num_runs
        runs(trial_idx:sum(num_TP_per_run(1:r))) = r;           %%%%%%%% this will have to change for RRI
        trial_idx = trial_idx + num_TP_per_run(r);
    end 

    subj = set_mat(subj,'selector','runs',runs);

    % detrend the timeseries data
    subj = detrend_runs(subj,'spiral','runs');  % not in mvpa tutorial, but seems important to do

    % move pattern to hard drive
    %subj = move_pattern_to_hd(subj, 'spiral_d');
    
    % clean up workspace
    subj = remove_mat(subj,'pattern','spiral');
    
    % high-pass filter the timeseries data
    subj = hpfilter_runs(subj,'spiral_d','runs',100,2); % remove frequencies below .01 Hz

     % clean up workspace
    subj = remove_mat(subj,'pattern','spiral_d');
    
    % move pattern to hard drive
    %subj = move_pattern_to_hd(subj, 'spiral_d_hp');
    
    % zscore the data from each run
    subj = zscore_runs(subj,'spiral_d_hp','runs'); % Z-score the data
    
    % clean up workspace
    subj = remove_mat(subj,'pattern','spiral_d_hp');
    
    % save final pattern in single precision form (8 sig figs) to save RAM and HD space    
    % subj.patterns{end}.mat = single(subj.patterns{end}.mat);