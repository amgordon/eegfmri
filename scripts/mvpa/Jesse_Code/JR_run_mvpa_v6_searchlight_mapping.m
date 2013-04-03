function []= JR_run_mvpa_v6_searchlight_mapping(subj_array, condition1, condition2, balance_per_subj_bin);

%%%%%%% specify user-defined variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for b=subj_array
    tic
    if b<10
        subj_id=strcat('s10', num2str(b))
    else
        subj_id=strcat('s1', num2str(b))
    end

    if balance_per_subj_bin==1
        balanced_or_unbal='balanced_bins';
    else
        balanced_or_unbal='unbalanced_bins';
    end

    PAST_dir = '/Users/Jesse/fMRI/data/PAST/fMRI';
    mvpa_dir = [PAST_dir '/' subj_id '/mvpa'];
    searchlight_dir = '/Users/Jesse/fMRI/data/PAST/fMRI/mvpa_results/searchlight_maps/searchlight_maps_unsmoothed_FIXED_AAL';
    
    if ~exist(searchlight_dir,'dir')
        mkdir(searchlight_dir);
    end

    load([PAST_dir '/vol_info.mat']); %get functional data resolution info for spm .img writing
    load(['/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/' subj_id '/onsets.mat']);

    roi_name = 'merged_AAL_ROIs_FIXED_HOLES';
    data_imgs_to_use = 'raw_filenames_wa.mat';
    %data_imgs_to_use = 'raw_filenames_wa.mat';

    mvpa_workspace = [PAST_dir '/' subj_id '/mvpa/' subj_id '_merged_AAL_ROIs_FIXED_HOLES_wa.m.mat'];

    num_results_iter = 5;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set flags (% 1 = yes please, 0 = no thanks)
    flags.equate_number_of_old_new_trials_per_subjective_bin = balance_per_subj_bin; % equate_per_subjective_bin;
    flags.equate_number_of_trials_in_cond_1_and_2 = 1;
    flags.anova_p_thresh = 1;  %p threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)
    flags.perform_second_round_of_zscoring = 0;
    flags.remove_artdetect_outliers = 1; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
    flags.artdetect_motion_thresh = 4;
    flags.artdetect_global_signal_thresh = 4;
    flags.remove_outlier_trials = 4;  % how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)
    flags.plot_ROC_curve = 0;
    flags.display_performance_breakdown = 1;
    flags.generate_importance_maps = 1;
    flags.write_data_log_to_text_file=1;

    % specify which conditions to use for classification (must correspond to the names of conditions specified below)
    condnames =  {condition1, condition2};

    TRs_to_average_over = [1 2 3 4 5]; %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier
    TR_weights = [0 0 .5 .5 0];
    %TR_weights = [.0072 .2168 .3781 .2742 .1237];  % from SPM canonical values at 1,3,5,7,and 9 sec post-stimulus

    % classifier parameters
    class_args.train_funct_name = 'train_bp';
    class_args.test_funct_name = 'test_bp';
    class_args.nHidden = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if ~exist(mvpa_workspace,'file')
        [subj num_runs num_TP_per_run]= JR_generate_mvpa_workspace_mat_file(subj_id, roi_name, data_imgs_to_use, mvpa_dir); % generate and save workspace
    else
        eval(['load ' mvpa_workspace])  %load workspace
    end

    subj_orig = subj; % save copy of original subj struct

    start_run = 1;

    for n = 1: num_results_iter

        subj = subj_orig; % overwrite subj struct w/ original

        % Extract info about conditions from onsets file
        num_conds = size(onsets,2);
        all_regs = zeros(num_conds,num_runs*num_TP_per_run); % initialize regs matrix as conditions x timepoints

        for cond = 1: num_conds-1 %(exclude last condition ("no_response" trials)
            for trial = 1: length(onsets{cond})
                time_idx = onsets{cond}(trial)/2+1; % divide by 2 and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
                all_regs(cond,time_idx) = 1;
            end
        end

        % SPECIAL FIX FOR because reg #6 (NEW_recollect) is a fake
        % placeholder trial for some subjects
        if ~isempty(find(strcmp(subj_id,{'s104','s105','s108','s113','s117'})))
            all_regs(6,:)=0;
        end

        % SPECIAL FIX FOR because reg #7 (NEW_HC_old) is a fake
        % placeholder trial for some subjects
        if ~isempty(find(strcmp(subj_id,{'s117'})))
            all_regs(7,:)=0;
        end

        % condense regs by removing zeros
        % initialize variables
        condensed_regs_all = [];
        condensed_runs = [];

        trial_counter = 1;
        for i = 1: size(all_regs,2)
            if ~isempty(find(all_regs(:,i))) % if not a rest timepoint
                %condensed_regs_of_interest(:,trial_counter) = regs_of_interest(:,i);
                condensed_regs_all(:,trial_counter) = all_regs(:,i);
                condensed_runs(1,trial_counter) = subj.selectors{1}.mat(i);
                trial_counter = trial_counter + 1;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % select TRs of interest (to correspond with peak post-stim BOLD response)

        all_trials = sum(all_regs,1); % vector of all trial   (patterns{4} contains fully preprocessed data)
        data_by_TR(1,:,:) = TR_weights(1)*subj.patterns{end}.mat(:,find(all_trials)+0); % 1st TR (0-2 sec)
        data_by_TR(2,:,:) = TR_weights(2)*subj.patterns{end}.mat(:,find(all_trials)+1); % 2nd TR (2-4 sec)
        data_by_TR(3,:,:) = TR_weights(3)*subj.patterns{end}.mat(:,find(all_trials)+2); % 3rd TR (4-6 sec)
        data_by_TR(4,:,:) = TR_weights(4)*subj.patterns{end}.mat(:,find(all_trials)+3); % 4th TR (6-8 sec)
        data_by_TR(5,:,:) = TR_weights(5)*subj.patterns{end}.mat(:,find(all_trials)+4); % 5th TR (8-10 sec)
        temporally_condensed_data = squeeze(sum(data_by_TR(TRs_to_average_over,:,:),1));

        clear data_by_TR; %clean up matlab workspace to save memory
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Exclude trials determined to be outliers by ArtDetect script
        % Guide to outlier file cell arrays...
        % Movement thresholds: .2 .25 .3 .35 .4 .4 .5
        % Global signal thresholds: 2 2.5 3 3.5 4 4.5 5

        if flags.remove_artdetect_outliers == 1
            load([PAST_dir '/outlier_indices/' subj_id '_outlier_indices']); %load outlier indices

            m_outliers = movement_outlier_trials{flags.artdetect_motion_thresh};  % remove trials with more than .35mm/TR of movement
            gs_outliers = global_signal_outlier_trials{flags.artdetect_global_signal_thresh}; % remove trials with global signal change of +/- 3.5 SD from mean
            combined_outliers = union(m_outliers,gs_outliers);

            %temporally_condensed_data(:,combined_outliers)=[]; % remove outlier trials from fmri data
            %condensed_regs_of_interest(:,combined_outliers) = [];
            condensed_regs_all(:,combined_outliers) = 0;
            %condensed_runs(combined_outliers) = [];

            display([num2str(length(m_outliers)) ' movement outlier trials flagged']);
            display([num2str(length(gs_outliers)) ' global signal outlier trials flagged']);
            display([num2str(length(combined_outliers)) ' total outlier trials excluded']);

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Artificially balance the number of trials given each subjective memory response
        if flags.equate_number_of_old_new_trials_per_subjective_bin == 1
            for j = 1:5
                OLD_trials = find(condensed_regs_all(j,:));
                NEW_trials = find(condensed_regs_all(j+5,:));

                num_OLD = length(OLD_trials);
                num_NEW = length(NEW_trials);

                if num_OLD > num_NEW
                    rand_array = rand(1,num_OLD);
                    [sorted inds]= sort(rand_array);
                    trials_to_cut = OLD_trials(inds(1:num_OLD-num_NEW));
                    condensed_regs_all(j,trials_to_cut)=0;
                elseif num_OLD < num_NEW
                    rand_array = rand(1,num_NEW);
                    [sorted inds]= sort(rand_array);
                    trials_to_cut = NEW_trials(inds(1:num_NEW-num_OLD));
                    condensed_regs_all(j+5,trials_to_cut)=0;
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % define conditions of interest
        % specify names 'OLD_recollect'    'OLD_hc_old'    'OLD_lc_old'    'OLD_lc_new'    'OLD_hc_new'    'NEW_recollect' 'NEW_hc_old'    'NEW_lc_old'    'NEW_lc_new'    'NEW_hc_new'    'no_resp'

        Objective_old = sum(condensed_regs_all(1:5,:));
        Objective_new = sum(condensed_regs_all(6:10,:));

        Subjective_old = sum(condensed_regs_all([1 2 3 6 7 8],:));
        Subjective_new = sum(condensed_regs_all([4 5 9 10],:));

        Subjective_old_HC_only = sum(condensed_regs_all([1 2 6 7],:));
        Subjective_new_HC_only = sum(condensed_regs_all([5 10],:));

        Hits = sum(condensed_regs_all([1 2 3],:));
        Misses = sum(condensed_regs_all([4 5],:));
        CRs = sum(condensed_regs_all([9 10],:));
        FAs = sum(condensed_regs_all(6:8,:));

        R_hits = condensed_regs_all(1,:);
        HC_hits = condensed_regs_all(2,:);
        LC_hits = condensed_regs_all(3,:);
        LC_misses = condensed_regs_all(4,:);
        HC_misses = condensed_regs_all(5,:);
        R_FAs = condensed_regs_all(6,:);
        HC_FAs = condensed_regs_all(7,:);
        LC_FAs = condensed_regs_all(8,:);
        LC_CRs = condensed_regs_all(9,:);
        HC_CRs = condensed_regs_all(10,:);

        R_and_HC_hits = sum(condensed_regs_all([1 2],:));

        %no_resp = condensed_regs_all(11,:); %excluded from analysis

        %assign conditions to train/test classifier on
        condensed_regs_of_interest = [];
        eval(['condensed_regs_of_interest(1,:) = ' condnames{1} ';'])
        eval(['condensed_regs_of_interest(2,:) = ' condnames{2} ';'])

        if flags.equate_number_of_trials_in_cond_1_and_2 == 1

            cond1_trials = find(condensed_regs_of_interest(1,:));
            cond2_trials = find(condensed_regs_of_interest(2,:));
            num_cond1 = length(cond1_trials);
            num_cond2 = length(cond2_trials);

            if num_cond1 > num_cond2
                rand_array = rand(1,num_cond1);
                [sorted inds]= sort(rand_array);
                trials_to_cut = cond1_trials(inds(1:num_cond1-num_cond2));
                condensed_regs_of_interest(1,trials_to_cut) = 0;
                display([num2str(length(trials_to_cut)) ' trials cut from ' condnames{1}]);
            elseif num_cond1 < num_cond2
                rand_array = rand(1,num_cond2);
                [sorted inds]= sort(rand_array);
                trials_to_cut = cond2_trials(inds(1:num_cond2-num_cond1));
                condensed_regs_of_interest(2,trials_to_cut) = 0;
                display([num2str(length(trials_to_cut)) ' trials cut from ' condnames{2}]);
            else
                display('Trial numbers are already balanced');
            end
        end

        display([num2str(count(condensed_regs_of_interest(1,:)==1)) ' trials in condition ' condnames{1}])
        display([num2str(count(condensed_regs_of_interest(2,:)==1)) ' trials in condition ' condnames{2}])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % initialize regressors object
        subj = init_object(subj,'regressors','conds');
        subj = set_mat(subj,'regressors','conds',condensed_regs_of_interest);
        subj = set_objfield(subj,'regressors','conds','condnames',condnames);

        % add new condensed activation pattern
        subj = duplicate_object(subj,'pattern','spiral_d_hp_z','spiral_d_hp_z_condensed');
        subj = set_mat(subj,'pattern','spiral_d_hp_z_condensed',temporally_condensed_data,'ignore_diff_size',true);

        zhist = sprintf('Pattern ''%s'' created by JR custom code','spiral_d_hp_z_condensed');
        subj = add_history(subj,'pattern','spiral_d_hp_z_condensed',zhist,true);

        % clean up workspace
        subj = remove_mat(subj,'pattern','spiral_d_hp_z');
        clear mean_data;

        % update run vector to condensed format
        subj.selectors{1}.mat = condensed_runs;
        subj.selectors{1}.matsize = size(condensed_runs);

        % "activate" only those trials of interest (from regs_of_interest) before
        % creating cross-validation indices
        active_trials = find(sum(condensed_regs_of_interest));

        if flags.remove_outlier_trials ~= 0
            % remove outlier trials (timepoints)
            mean_across_voxels = mean(subj.patterns{end}.mat(:,active_trials),1);
            z_mean_across_voxels = zscore(mean_across_voxels);
            upper_outliers = find(z_mean_across_voxels> flags.remove_outlier_trials);
            lower_outliers = find(z_mean_across_voxels< -1 * flags.remove_outlier_trials);
            all_outliers = union(upper_outliers,lower_outliers)
            active_trials(all_outliers) = [];
        end

        actives_selector = zeros(1,size(condensed_regs_all,2)); % intialize vector of all zeros
        actives_selector(active_trials) = 1; % remove all non-"regs_of_interst" trials (set to one)
        subj = init_object(subj,'selector','conditions_of_interest'); %initialize selector object
        subj = set_mat(subj,'selector','conditions_of_interest',actives_selector);
        subj = create_xvalid_indices(subj,'runs','actives_selname','conditions_of_interest');

        for q = 3:num_runs+2;
            subj.selectors{q}.mat(subj.selectors{q}.mat==2)=3; % change 2's to 3's
        end

        if flags.perform_second_round_of_zscoring == 1
            % zscore temporally-condensed data; active trials only (second round of z-scoring)
            %             for r=1:length(active_trials)
            %                 subj.patterns{end}.mat(:,active_trials(r)) = subj.patterns{end}.mat(:,active_trials(r))-mean(subj.patterns{end}.mat(:,active_trials),2);
            %             end
            subj.patterns{end}.mat(:,active_trials) = zscore(subj.patterns{end}.mat(:,active_trials)')';
            display('Performing second round of z-scoring')
        end

        searchlight_radius = 3;
        subj.adj_sphere = create_adj_list(subj,roi_name,'radius',searchlight_radius);

        %specify training/testing functions for within-sphere classification
        class_args.train_funct_name = 'train_gnb';
        class_args.test_funct_name = 'test_gnb';

        scratch.class_args = class_args;
        scratch.perfmet_funct = 'perfmet_maxclass';
        scratch.perfmet_args = struct([]);

        statmap_srch_arg.adj_list = subj.adj_sphere;
        statmap_srch_arg.obj_funct = 'statmap_classify';
        statmap_srch_arg.scratch = scratch;

        subj = feature_select( ...
            subj, ...
            'spiral_d_hp_z_condensed', ... % data
            'conds', ... % binary regs (for GNB)
            'runs_xval', ... % selector
            'statmap_funct','statmap_searchlight', ... % function
            'statmap_arg',statmap_srch_arg, ...
            'new_map_patname','spiral_d_hp_z_condensed_srch', ...
            'thresh',[]);

        % create masks from the statmaps, by picking the best N values (e.g., 500) in each
        %statmap to use for classification on the remaining run of test trials
        %NOTE: keeps top N voxels; not top N spheres
        %         subj = create_sorted_mask( ...
        %             subj,'spiral_d_hp_z_condensed_srch', ...
        %             'spiral_d_hp_z_condensed_srch_1000',1000, ...
        %             'descending',true);


        %% WRITE OUT MEAN SEARCHLIGHT MAP TO .IMG FILE
        % average searchlight performance maps across runs
        for r=1:num_runs
            run_index = start_run + r - 1;
            sl_voxel_values(:,run_index)=subj.patterns{end-num_runs+r}.mat;
        end

        start_run = start_run+num_runs; % update start run for next results iter
    end

    sl_mean_voxel_values = mean(sl_voxel_values,2);

    vol_info.fname = [searchlight_dir '/' subj_id '_' condnames{1} '_vs_' condnames{2} '_3vox_radius_searchlight.img'];

    sl_map = zeros(vol_info.dim);
    included_voxels = find(subj.masks{1}.mat);
    sl_map(included_voxels) = sl_mean_voxel_values.*100; %multiply by 100 to improve scaling for visualization

    spm_write_vol(vol_info, sl_map);

    clear sl_voxel_values;

    time2finish = toc/60;
    display(['Finished ' subj_id ' in ' num2str(time2finish) ' minutes']);
end
