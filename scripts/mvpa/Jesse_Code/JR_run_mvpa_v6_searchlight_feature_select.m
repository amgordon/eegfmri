function [subj results] = JR_run_mvpa_v6_batch(subj_arr, cond1, cond2, balance_per_subj_bin)

tic
for b=subj_arr
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

    importance_maps_dir=[PAST_dir '/mvpa_results/improved_importance_maps/' cond1 '_vs_' cond2 '_' balanced_or_unbal];
    xls_results_data_logs_dir=[PAST_dir '/mvpa_results/xls_results_data_logs/' cond1 '_vs_' cond2 '_' balanced_or_unbal];

    if ~exist(importance_maps_dir,'dir')
        mkdir(importance_maps_dir);
    end
    if ~exist(xls_results_data_logs_dir, 'dir')
        mkdir(xls_results_data_logs_dir)
    end

    mvpa_workspace = [PAST_dir '/' subj_id '/mvpa/' subj_id '_merged_AAL_ROIs_8mm_smoothing.mat'];
    %%%%%%% specify user-defined variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    flags.save_workspace = 1; % 1 = yes please, 0 = no thanks

    num_results_iter = 10; % number of times to run the cross validation process
    for x = 1: num_results_iter

        % unless a previously created workspace is specified, load in the preprocessed data
        if ~exist('mvpa_workspace')


            exp_name = 'PAST';
            roi_file = '/Users/Jesse/fMRI/data/PAST/fMRI/4x4x4/merged_AAL_ROIs.nii';
            roi_name = 'merged_AAL_ROIs'
            num_TP_per_run = 203;

            % load user-created filename and onsets lists into workspace
            data_imgs_to_use = 'raw_filenames_s8mm_wa.mat';
            load([mvpa_dir '/' data_imgs_to_use]); %loads predefined cell array called raw_filenames into memory
            num_runs = length(raw_filenames)/num_TP_per_run; %calculate the number of runs (allows script to work flexibly for subjects with missing runs)
            load(['/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/' subj_id '/onsets.mat']);
            load('/Users/Jesse/fMRI/data/PAST/fMRI/vol_info.mat'); %get functional data resolution info for spm .img writing
            [subj] = JR_mvpa_load_and_preprocess_raw_data(subj_id, exp_name, roi_name, roi_file, raw_filenames, num_runs, num_TP_per_run);

            if flags.save_workspace == 1
                save_cmd = ['save ' subj_id '_' roi_name '_8mm_smoothing.mat'];
                eval(save_cmd);
            end
        else
            eval(['load ' mvpa_workspace])
            if exist('log')
                clear log;  % just in case an old 'data_log' variable exists in memory
                display('data_log variable cleared from memory');
            end
            load('/Users/Jesse/fMRI/data/PAST/fMRI/vol_info.mat'); %get functional data resolution info for spm .img writing

        end

        % Set flags (% 1 = yes please, 0 = no thanks)
        flags.equate_number_of_old_new_trials_per_subjective_bin =balance_per_subj_bin; % equate_per_subjective_bin;
        flags.equate_number_of_trials_in_cond_1_and_2 = 1;
        flags.perform_second_round_of_zscoring = 0;
        flags.remove_artdetect_outliers = 1; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
        flags.remove_outlier_trials = 0;  % how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)
        flags.plot_ROC_curve = 0;
        flags.display_performance_breakdown = 1;
        flags.generate_importance_maps = 0;
        flags.write_data_log_to_text_file=1;


        % specify which conditions to use for classification (must correspond to the names of conditions specified below)
        condnames =  {cond1 , cond2};

        data_imgs_to_use = 'raw_filenames_s4mm_wa.mat';

        TRs_to_average_over = [1 2 3 4 5]; %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier
        TR_weights = [0 0 .5 .5 0];
        %TR_weights = [.0072 .2168 .3781 .2742 .1237];  % from SPM canonical values at 1,3,5,7,and 9 sec post-stimulus

        num_results_iter = 5; % number of times to run the cross validation process

        anova_p_thresh = 1;  %p threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)


        % classifier parameters
        class_args.train_funct_name = 'train_bp';
        class_args.test_funct_name = 'test_bp';
        class_args.nHidden = 0;



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

        condensed_regs_of_interest = []; % initialize variables
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

        %     if length(TRs_to_average_over)==1
        %         temporally_condensed_data = squeeze(data_by_TR(TRs_to_average_over,:,:)); %remove singleton dimension
        %     else
        %         temporally_condensed_data = squeeze(mean(data_by_TR(TRs_to_average_over,:,:),1));
        %     end

        clear data_by_TR; %clean up matlab workspace to save memory
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % Exclude trials determined to be outliers by ArtDetect script
        % Guide to outlier file cell arrays...
        % Movement thresholds: .2 .25 .3 .35 .4 .4 .5
        % Global signal thresholds: 2 2.5 3 3.5 4 4.5 5

        if flags.remove_artdetect_outliers == 1
            load([PAST_dir '/outlier_indices/' subj_id '_outlier_indices']); %load outlier indices

            m_outliers = movement_outlier_trials{3};  % remove trials with more than .3mm/TR of movement
            gs_outliers = global_signal_outlier_trials{3}; % remove trials with global signal change of +/- 3 SD from mean
            combined_outliers = union(m_outliers,gs_outliers);

            temporally_condensed_data(:,combined_outliers)=[]; % remove outlier trials from fmri data
            %condensed_regs_of_interest(:,combined_outliers) = [];
            condensed_regs_all(:,combined_outliers) = [];
            condensed_runs(combined_outliers) = [];
            
            display([num2str(length(m_outliers)) ' movement outlier trials excluded']);
            display([num2str(length(gs_outliers)) ' global signal outlier trials excluded']);

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
        subj = JR_create_xvalid_indices_searchlight(subj,'runs','actives_selname','conditions_of_interest');

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
        subj = create_sorted_mask( ...
            subj,'spiral_d_hp_z_condensed_srch', ...
            'spiral_d_hp_z_condensed_srch_1000',1000, ...
            'descending',true);


        % run classifier (No hidden layer NN Toolbox backprop algorithm)
        class_args.train_funct_name = 'train_bp';
        class_args.test_funct_name = 'test_bp';
        class_args.nHidden = 0;
              
        classifier_mask = 'spiral_d_hp_z_condensed_srch_1000'; % use original mask
                   
        
        for q = 3:num_runs+3;
            subj.selectors{q}.mat(subj.selectors{q}.mat==3)=1;
        end

        %%%%%%%%%%%%%%%%%%%%%% RUN THE CLASSIFIER (CROSS-VALIDATION)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [subj results{x}] = cross_validation(subj,'spiral_d_hp_z_condensed','conds','runs_xval',classifier_mask,class_args);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % do some important RAM clean-up and data archiving
        for y = 1:num_runs
            if flags.generate_importance_maps == 1
                results_IW{x}.iterations(y).scratchpad.net.IW{1} = results{x}.iterations(y).scratchpad.net.IW{1}; % save weights to pass to JR_interpret_weights
            end
            results{x}.iterations(y).scratchpad.net.inputs{1}.exampleInput=[]; % delete huge data object from results scratchpad to free up RAM
        end

        if flags.display_performance_breakdown == 1
            % analyze the results in more detail
            correct_vector = [];
            desireds_vector = [];
            guesses_vector = [];
            acts_diff_vector = [];
            for a = 1:num_runs
                correct_vector = horzcat(correct_vector,results{x}.iterations(a).perfmet.corrects);
                desireds_vector = horzcat(desireds_vector,results{x}.iterations(a).perfmet.desireds);
                guesses_vector = horzcat(guesses_vector,results{x}.iterations(a).perfmet.guesses);
                acts_diff_vector = horzcat(acts_diff_vector, results{x}.iterations(a).acts(1,:)-results{x}.iterations(a).acts(2,:));
            end       
            
            
            overall_accuracy = mean(correct_vector)
            overall_hit_rate = mean(correct_vector(desireds_vector==1))
            overall_fa_rate = 1-mean(correct_vector(desireds_vector==2))
            overall_d_prime = norminv(overall_hit_rate)-norminv(overall_fa_rate);


            % classification_accuracy_by_resp (for full classification only)
            for b = 1:10
                classification_accuracy_by_resp(b) = mean(correct_vector(find(condensed_regs_all(b,active_trials))));
                mean_acts_diffs_by_resp(b) = mean(acts_diff_vector(find(condensed_regs_all(b,active_trials))));
                number_of_trials_per_bin(b) = length(find(condensed_regs_all(b,active_trials)));
            end

            % display(classification_accuracy_by_resp);
            % display(mean_acts_diffs_by_resp);
            % display(number_of_trials_per_bin);

            % sort by absolute value of classifier "confidence"
            [abs_sorted_diffs abs_ind] = sort(abs(acts_diff_vector),2,'descend');
            abs_correct_sorted = correct_vector(abs_ind);
            abs_desireds_sorted = desireds_vector(abs_ind);

            % print results of top 50,100, 150, etc. trials, sorted by
            % classifier "confidence"
            num_trials = length(abs_correct_sorted);
            acc_sorted_by_classifier_confidence = [0 0 0 0 0 0 0 0];

            if num_trials>50

                for j= 1:floor(num_trials/50)
                    acc_sorted_by_classifier_confidence(j)=mean(abs_correct_sorted(1:50*j));
                end
                acc_sorted_by_classifier_confidence(j+1)=mean(abs_correct_sorted(1:end));

                %display(acc_sorted_by_classifier_confidence)
            end
        end
     end

        if flags.plot_ROC_curve == 1

            % sort by signed classifier "confidence" (for ROI curves)
            [sorted_diffs ind] = sort(acts_diff_vector,2,'descend');
            correct_sorted = correct_vector(ind);
            desireds_sorted = desireds_vector(ind);

            % create continuous ROC function
            for i = 1:length(sorted_diffs);
                hit_rate(i) = length(correct_sorted(intersect(find(desireds_sorted == 1),[1:i]))) / length(find(desireds_sorted == 1));
                fa_rate(i) = length(correct_sorted(intersect(find(desireds_sorted == 2),[1:i]))) / length(find(desireds_sorted == 2));
            end

            figure
            plot(fa_rate,hit_rate,'.-')
            hold on
            plot([0 1],[0 1],'r')
            xlabel('P(Old|New)')
            ylabel('P(Old|Old)')
        end

        if flags.write_data_log_to_text_file==1

            data_log.overall_acc(x)=overall_accuracy;
            data_log.hits(x)=overall_hit_rate;
            data_log.FAs(x)=overall_fa_rate;
            data_log.d_prime(x)=overall_d_prime;
            data_log.classification_accuracy_by_resp(x,:)=classification_accuracy_by_resp;
            data_log.number_trials_per_bin(x,:)=number_of_trials_per_bin;
            data_log.acc_sorted_by_classifier_confidence(x,:)=acc_sorted_by_classifier_confidence;
        end
    end


    if flags.write_data_log_to_text_file==1

        filename= [xls_results_data_logs_dir '/' subj_id '_' condnames{1} '_vs_' condnames{2} '.xls'];
        fid=fopen(filename, 'wt');
        fprintf(fid, '%s\r\n', ['subj_id = ' subj_id]);
        fprintf(fid, '%s\r\n', ['ROI_name = ' roi_name]);
        fprintf(fid, '%s\r\n', ['data_imgs_to_use =' data_imgs_to_use]);
        fprintf(fid, '%s\r\n', ['TR_weights = ' num2str(TR_weights)]);
        fprintf(fid, '%s\r\n', ['classification:' condnames{1} ' vs. ' condnames{2}]);
        fprintf(fid, '%s\r\n', ['flags.equate_number_of_trials_in_cond_1_and_2 = ' num2str(flags.equate_number_of_trials_in_cond_1_and_2)]);
        fprintf(fid, '%s\r\n', ['flags.equate_number_of_old_new_trials_per_subjective_bin = ' num2str(flags.equate_number_of_old_new_trials_per_subjective_bin)]);
        fprintf(fid, '%s\r\n', ['flags.perform_second_round_of_zscoring = ' num2str(flags.perform_second_round_of_zscoring)]);
        fprintf(fid, '%s\r\n', ['flags.remove_outlier_trials (std dev) = ' num2str(flags.remove_outlier_trials)]);


        fprintf(fid, '\n\n');
        fprintf(fid, 'results_iter\toverall_acc\tHits\tFAs\td_prime\tOLD_recollect\tOLD_hc_old\tOLD_lc_old\tOLD_lc_new\tOLD_hc_new\tNEW_recollect\tNEW_hc_old\tNEW_lc_old\tNEW_lc_new\tNEW_hc_new\tTop50\tTop100\tTop150\tTop200\tTop250\n');
        fprintf(fid, '#trials\t\t\t\t\t');
        fprintf(fid, '%3.0f\t', number_of_trials_per_bin(1:10));
        fprintf(fid, '\n');

        for q=1:x
            fprintf(fid, '%3.3f\t', q);
            fprintf(fid, '%3.3f\t', data_log.overall_acc(q));
            fprintf(fid, '%3.3f\t', data_log.hits(q));
            fprintf(fid, '%3.3f\t', data_log.FAs(q));
            fprintf(fid, '%3.3f\t', data_log.d_prime(q));
            fprintf(fid, '%3.3f\t', data_log.classification_accuracy_by_resp(q,:)); %[data_log.classification_accuracy_by_resp(1,q) '\t' data_log.classification_accuracy_by_resp(q,1) '\t' data_log.classification_accuracy_by_resp(q,1) '\t' data_log.classification_accuracy_by_resp(q,4) '\t' data_log.classification_accuracy_by_resp(q,5) '\t' data_log.classification_accuracy_by_resp(q,6) '\t' data_log.classification_accuracy_by_resp(q,7) '\t' data_log.classification_accuracy_by_resp(q,8) '\t' data_log.classification_accuracy_by_resp(q,9) '\t' data_log.classification_accuracy_by_resp(q,10) '\t']);
            fprintf(fid, '%3.3f\t', data_log.acc_sorted_by_classifier_confidence(q,:));
            fprintf(fid, '\n');
            %for p=1:size(data_log.acc_sorted_by_classifier_confidence, 1)
            %    fprintf(fid, '%3.3f\t', data_log.acc_sorted_by_classifier_confidence(q,p));
            %end

        end
        fprintf(fid, '%s\t', 'mean');
        fprintf(fid, '%3.3f\t', mean(data_log.overall_acc));
        fprintf(fid, '%3.3f\t', mean(data_log.hits));
        fprintf(fid, '%3.3f\t', mean(data_log.FAs));
        fprintf(fid, '%3.3f\t', mean(data_log.d_prime));
        fprintf(fid, '%3.3f\t', mean(data_log.classification_accuracy_by_resp,1));
        fprintf(fid, '%3.3f\t', mean(data_log.acc_sorted_by_classifier_confidence,1));
        fprintf(fid, '\n');
        fclose(fid);
    end



    if flags.generate_importance_maps == 1;

        subj = JR_interpret_weights_NEW(subj, results,results_IW);

        impmap1 = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map 1
        impmap2 = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map 2

        if anova_p_thresh == 1 % NO ANOVA VERSION

            voxel_inds = find(subj.masks{end}.mat); %get mask voxel indices
            for j = 1:num_runs
                temp1 = zeros(vol_info.dim); %initialize appropriately sized matrix
                temp2 = zeros(vol_info.dim); %initialize appropriately sized matrix
                temp1(voxel_inds)=subj.patterns{end-num_runs+j}.mat(:,1); %store impmap values at appropriate voxel indices
                temp2(voxel_inds)=subj.patterns{end-num_runs+j}.mat(:,2); %store impmap values at appropriate voxel indices
                impmap1 = impmap1+temp1; %add values cumulatively across iterations
                impmap2 = impmap2+temp2;
            end

            impmap1 = impmap1/num_runs*1000; %compute average and multiply by 1000 for scaling
            impmap2 = impmap2/num_runs*1000; %compute average and multiply by 1000 for scaling

            vol_info.fname = [importance_maps_dir '/' subj_id '_' condnames{1} '.img'];
            spm_write_vol(vol_info,impmap1);
            vol_info.fname = [importance_maps_dir '/' subj_id '_' condnames{2} '.img'];
            spm_write_vol(vol_info,impmap2);

        else % ANOVA VERSION

            for j = 1:num_runs
                temp1 = zeros(vol_info.dim); %initialize appropriately sized matrix
                temp2 = zeros(vol_info.dim); %initialize appropriately sized matrix
                voxel_inds{j} = find(subj.masks{end-num_runs+j}.mat); %get mask voxel indices
                temp1(voxel_inds{j})=subj.patterns{end-num_runs+j}.mat(:,1); %store impmap values at appropriate voxel indices
                temp2(voxel_inds{j})=subj.patterns{end-num_runs+j}.mat(:,2); %store impmap values at appropriate voxel indices
                impmap1 = impmap1+temp1; %add values cumulatively across iterations
                impmap2 = impmap2+temp2;
            end

            %sum across masks to get composite mask (where value of each voxel =
            %number of runs for which that voxel was included)
            composite_mask = zeros(vol_info.dim);
            for i = 2:size(subj.masks,2)  %exclude first mask (it's the starting ROI)
                composite_mask = composite_mask+subj.masks{i}.mat;
            end
            voxels_to_exclude = find(composite_mask<6);  % exclude voxels that exist for less than 6 of the ANOVA masks
            impmap1(voxels_to_exclude)=0;
            impmap2(voxels_to_exclude)=0;

            impmap1_avg = impmap1./composite_mask * 1000;  % divide by number of observations contributing to each sum (to get avg) and multiply by 1000 for scaling
            impmap2_avg = impmap2./composite_mask * 1000;  % divide by number of observations contributing to each sum (to get avg) and multiply by 1000 for scaling

            vol_info.fname = [importance_maps_dir '/' subj_id '_' condnames{1} '_p' num2str(anova_p_thresh) '.img'];
            spm_write_vol(vol_info,impmap1_avg);
            vol_info.fname = [importance_maps_dir '/' subj_id '_' condnames{2} '_p' num2str(anova_p_thresh) '.img'];
            spm_write_vol(vol_info,impmap2_avg);
        end

    end

    time2finish = toc/60;
    display(['Finished ' subj_id ' in ' num2str(time2finish) ' minutes']);

end
