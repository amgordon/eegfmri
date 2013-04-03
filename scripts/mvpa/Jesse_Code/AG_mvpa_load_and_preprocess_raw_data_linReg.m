function [subj] = AG_mvpa_load_and_preprocess_raw_data(S)

    % initialize subj structure
    subj = init_subj(S.exp_name,S.subj_id);

    % load mask file
    subj = load_spm_mask(subj,S.roi_name,S.roi_file);
    
    % load functional data
    subj = load_analyze_pattern(subj,'spiral',S.roi_name, S.img_files,'single',true);

    % move pattern to hard drive
    % subj = move_pattern_to_hd(subj, 'spiral');

    % make runs vector
    subj = init_object(subj,'selector','runs');

%     trial_idx = 1;
%     for r = 1:num_runs
%         runs(trial_idx:num_TP_per_run*r) = r;
%         trial_idx = trial_idx + num_TP_per_run;
%     end 

    trial_idx = 1;
    for r = 1:length(S.runs_vector)
        runs(trial_idx:trial_idx+S.runs_vector(r)-1)=r;
        trial_idx = trial_idx+S.runs_vector(r);
    end
    
    
    
    subj = set_mat(subj,'selector','runs',runs);
    
    
    
    subj = AG_remove_zero_voxels(subj, 'spiral', 'runs',S.roi_name);
    
    subj = remove_mat(subj,'pattern','spiral');
    
    % detrend the timeseries data
    subj = detrend_runs(subj,'spiral_n0','runs');  % not in mvpa tutorial, but seems important to do
    
    % move pattern to hard drive
    %subj = move_pattern_to_hd(subj, 'spiral_d');
    
    % clean up workspace
    subj = remove_mat(subj,'pattern','spiral_n0');
    
    % high-pass filter the timeseries data
    subj = hpfilter_runs(subj,'spiral_n0_d','runs',100,2); % remove frequencies below .01 Hz
    
    % clean up workspace
    subj = remove_mat(subj,'pattern','spiral_n0_d');
    
    % move pattern to hard drive
    %subj = move_pattern_to_hd(subj, 'spiral_d_hp');
    
    % zscore the data from each run
    subj = zscore_runs(subj,'spiral_n0_d_hp','runs'); % Z-score the data
    
    % clean up workspace
    subj = remove_mat(subj,'pattern','spiral_n0_d_hp');  

    
    
    
    if ~exist(S.mvpa_dir)
        mkdir(S.mvpa_dir);
    end
    
    save (S.workspace);
    %save (fullfile(S.mvpa_dir, [S.subj_id '_' S.roi_name '.mat']));
    
    %save_cmd = ['save ' S.mvpa_dir '/' S.subj_id '_' S.roi_name '.mat'];
    %eval(save_cmd);
    
    % save final pattern in single precision form (8 sig figs) to save RAM and HD space    
    % subj.patterns{end}.mat = single(subj.patterns{end}.mat);