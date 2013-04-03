function [qqq, results]= DKSort_run_mvpa_general(subj_array, saveName)





if nargin==2
    S.saveName= saveName;
else
    S.saveName = input('What do you want to name the results structure?', 's');
end




S.subj_array = subj_array;


%S.condsOfInterestList = {{'CorTrial_match_Session1', 'CorTrial_nonmatch_Session1'}; {'CorTrial_same_Session1', 'CorTrial_diff_Session1'};  {'CorTrial_shape_Session1', 'CorTrial_color_Session1'}; {'CorTrial_shape_Session1', 'CorTrial_pattern_Session1'}; {'CorTrial_color_Session1', 'CorTrial_pattern_Session1'}};
S.condsOfInterestList = {{'CorTrial_match', 'CorTrial_nonmatch'}; {'CorTrial_same', 'CorTrial_diff'};{'CorTrial_shape', 'CorTrial_color'}; {'CorTrial_shape', 'CorTrial_pattern'};{'CorTrial_color', 'CorTrial_pattern'}};
%S.condsOfInterestList = {{'finger1' 'finger2'}} ;
%S.condsOfInterestList = {{'cues' 'trials'}} ;
%S.condsOfInterestList = {{'cue_shape_same', 'cue_shape_diff'} {'cue_color_same', 'cue_color_diff'}  {'cue_pattern_same', 'cue_pattern_diff'}};
S.PenaltyParams = [100];
S.featureSels = [1000];
%%%%%%% specify user-defined variables
%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for FS = 1:length(S.featureSels)
    
    for penal = 1:length(S.PenaltyParams)
        
        for CI = 1:length(S.condsOfInterestList)
            
            for b=(1:length(subj_array))  % allows subject ID #'s to be specified in short format
                tic
                
                par = DKSortParams(subj_array{b});
                
                S.wascanfiles = mat2cell(par.allwascanfiles, [ones(1,size(par.allwascanfiles,1))], [size(par.allwascanfiles,2)]);
                S.swascanfiles = mat2cell(par.allswascanfiles, [ones(1,size(par.allswascanfiles,1))], [size(par.allswascanfiles,2)]);
                
                
     
                
                
                
                S.runs_vector =  par.numvols;%[594 594 594 594];
                S.meta_runs = S.runs_vector;%par.numvols;%[594 594 594 594]; %change this to reflect different amount of runs in different DKSort subjects
                %%%%%changed May 13 2010
                
                
                S.subj_id = subj_array{b};
                
                S.condsOfInterest = S.condsOfInterestList{CI};
                %s
                S.expt_dir = '/Users/alangordon/Studies/AG3/DK_sort/fmri_data';
                S.mvpa_dir = [S.expt_dir '/' S.subj_id '/mvpa_main_effects'];
                %S.univar_dir = [S.expt_dir '/' S.subj_id '/mvpa_analysis_withCues'];
                S.importance_maps_dir=[S.expt_dir '/mvpa_results/ImpMaps_' date  ];
                S.impType = {'pos' 'neg' 'both' 'raw'};
                S.group_mvpa_dir = [S.expt_dir '/mvpa_files'];
                
                
                
                S.condnames = {S.condsOfInterest{1}  S.condsOfInterest{2}};
                % create these directories if they don't already exist
                if ~exist(S.importance_maps_dir,'dir')
                    mkdir(S.importance_maps_dir);
                end
                
                
                S.num_runs = length(S.runs_vector);
                S.num_vols = sum(S.runs_vector);
                
                %1 for unsmoothed, 0 for smoothed;
                S.use_unsmoothed = 1;
                S.smoothTxt = {'smoothed' 'unsmoothed'};
                

                
                
                % OPTIONAL:  specify previously saved mvpa workspace to bypass
                % time-consuming data extraction and preprocessing
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Set flags (% unless otherwise indicated: 1 = yes please, 0 = no thanks)
                flags.use_premade_workspace = 1;
                flags.num_results_iter = 1; % number of times to run the entire classification process (select subset of the data and train/test classifier)
                flags.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data
                S.equate_number_of_trials_in_cond_1_and_2 = 1; % equate number of trials in conditions 1 and 2 (RECOMMENDED)
                flags.anova_p_thresh = 1;  % p-value threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)
                flags.perform_second_round_of_zscoring = 0;  % z-score data again immediately prior to classification (NOT RECOMMENDED)
                flags.remove_artdetect_outliers = 0; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
                flags.artdetect_motion_thresh = 0; % specify ArtDetect bin for motion outliers
                flags.artdetect_global_signal_thresh = 0; % specify ArtDetect bin for global signal outliers
                flags.remove_outlier_trials = 0;  % on-the-fly outlier detection/removal; specify how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)
                flags.plot_ROC_curve = 0;
                flags.display_performance_breakdown = 0;
                flags.generate_importance_maps = 0;
                flags.write_data_log_to_text_file=0;
                flags.mask = 2;
                flags.useIdealParams = 0;
                flags.betaSeriesAnalysis = 0;
                %anova_nVox_thresh = 100;
                
                % load some .mat files into memory
                S.vol_info = load([S.expt_dir '/vol_info.mat']); %get functional data resolution info for spm .img writing
                
                [S.onsets, S.img_files] = DKSort_mvpa_onsets_and_images(S);
                qqq.onsets{b} = S.onsets;
                
                S.num_conds = size(S.onsets,2);
                
                S.exp_name = 'DKSort';
                
                
                if flags.mask==1
                    S.roi_name = 'mask.img';
                    S.roi_file = [S.univar_dir '/' S.roi_name];
                elseif flags.mask == 2
                    S.roi_name = 'rPFC1.img';
                    S.roi_file = ['/Users/alangordon/Studies/AG3/DK_sort/masks/PFCMask_051110/' S.roi_name];
                end
                
                S.TR_weights_idx = [0 0  .5 .5]; %try with [0 0 0 .5 .5]
                
                S.TRs_to_average_over_train = 1:length(S.TR_weights_idx); %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier
                S.TR_weights_train = S.TR_weights_idx; % should sum to 1
                
                S.TRs_to_average_over_test = 1:length(S.TR_weights_idx); %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier
                S.TR_weights_test = S.TR_weights_idx; % should sum to 1
                %S.TR_weights = [.0072 .2168 .3781 .2742 .1237];  % from SPM canonical
                %values at 1,3,5,7,and 9 sec post-stimulus
                
                % classifier parameters
                S.class_args.train_funct_name = 'train_pLR';
                S.class_args.test_funct_name = 'test_pLR';
                S.class_args.nHidden = 0;
                S.class_args.penalty = S.PenaltyParams(penal);
                S.class_args.nVox = S.featureSels(FS);  %S.Ps(pP); % 0 = no feature selection
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                S.workspace = fullfile(S.mvpa_dir, [S.subj_id '_' S.roi_name '_' S.smoothTxt{S.use_unsmoothed + 1} '.mat']);
                
                existWorkspace = exist(S.workspace);
                
                Sequal = 0;
                if existWorkspace
                    WorkCheck = load(S.workspace);
                    Sequal = isequal(struct2cell(WorkCheck.S), struct2cell(S));
                end
                
                
                if (flags.use_premade_workspace&&existWorkspace) %(flags.use_premade_workspace && Sequal) %
                    load(S.workspace, 'subj');
                else
                    if flags.betaSeriesAnalysis == 1;
                        [subj]= AG_mvpa_load_betaSeries;
                    else
                        
                        [subj] = AG_mvpa_load_and_preprocess_raw_data(S);
                    end
                end
                
                S.thisSelector = 'meta_runs_condensed_xval';
                
                
                subj_orig = subj; % save copy of original subj struct
                
                x = 0; % initialize results counter x to 0
                
                for n = 1: flags.num_results_iter
                    
                    subj = subj_orig; % overwrite subj struct w/ original
                    
                    % Extract info about conditions from onsets file
                    
                    all_regs = zeros(S.num_conds,S.num_vols); % initialize regs matrix as conditions x timepoints
                    
                    for cond = 1: S.num_conds
                        for trial = 1: length(S.onsets{cond})
                            time_idx = round(S.onsets{cond}(trial))/2 + 1; % divide by 2 and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
                            all_regs(cond,time_idx) = 1;
                        end
                    end
                    
                    % condense regs by removing zeros
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
                    
                    idx_condense =find(sum(all_regs));
                    
                    
                    %Create meta-runs selector (e.g. one meta run is 'Enc', one is
                    %'Ret')
                    trial_idx = 1;
                    m_runs = 0;
                    for r = 1:length(S.meta_runs)
                        m_runs(trial_idx:trial_idx+S.meta_runs(r)-1)=r;
                        trial_idx = trial_idx+S.meta_runs(r);
                    end
                    
                    meta_runs_condensed = m_runs(idx_condense);
                    
                    trainTrials =  meta_runs_condensed==1;
                    testTrials = meta_runs_condensed==2;
                    
                    all_trials = sum(all_regs,1);
                    
                    
                    meta_runs_train = find(all_trials);
                    meta_runs_test = [];
                    
                    
                    subj = init_object(subj,'selector','meta_runs_condensed');
                    subj = set_mat(subj,'selector','meta_runs_condensed', meta_runs_condensed);
                    subj = create_xvalid_indices(subj,'meta_runs_condensed');
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % select TRs of interest (to correspond with peak post-stim BOLD response)
                    
                    
                    data_by_TR_train = [];
                    for dt = 1:length(S.TR_weights_idx)
                        data_by_TR_train(dt,:,:) = S.TR_weights_train(dt)*subj.patterns{end}.mat(:,meta_runs_train+(dt-1));
                    end
                    
                    temporally_condensed_data_train = squeeze(sum(data_by_TR_train(S.TRs_to_average_over_train,:,:),1));
                    
                    clear data_by_TR_train
                    
                    
                    
                    data_by_TR_test = [];
                    for dt = 1:length(S.TR_weights_idx)
                        data_by_TR_test(dt,:,:) = S.TR_weights_test(dt)*subj.patterns{end}.mat(:,meta_runs_test+(dt-1));
                    end
                    
                    temporally_condensed_data_test = squeeze(sum(data_by_TR_test(S.TRs_to_average_over_test,:,:),1));
                    
                    
                    clear data_by_TR_test
                    
                    finalTrialCorrector = zeros(1,length(meta_runs_test))+length(meta_runs_test);
                    finalTrialCorrector(meta_runs_test<=2100-5) = meta_runs_test(meta_runs_test<=2100-5)+5;
                    
                    
                    temporally_condensed_data = horzcat(temporally_condensed_data_train, temporally_condensed_data_test);
                    
                    clear temporally_condensed_data_train;
                    clear temporally_condensed_data_test;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % Exclude trials determined to be outliers by custom ArtDetect script
                    % Guide to outlier file cell arrays...
                    % Movement thresholds: .2 .25 .3 .35 .4 .4 .5
                    % Global signal thresholds: 2 2.5 3 3.5 4 4.5 5
                    
                    if flags.remove_artdetect_outliers == 1
                        load([expt_dir '/outlier_indices/' subj_id '_outlier_indices']); %load outlier indices
                        
                        m_outliers = movement_outlier_trials{flags.artdetect_motion_thresh};  % remove trials with more than .35mm/TR of movement
                        gs_outliers = global_signal_outlier_trials{flags.artdetect_global_signal_thresh}; % remove trials with global signal change of +/- 3.5 SD from mean
                        combined_outliers = union(m_outliers,gs_outliers);
                        
                        condensed_regs_all(:,combined_outliers) = 0;
                        
                        display([num2str(length(m_outliers)) ' movement outlier trials flagged']);
                        display([num2str(length(gs_outliers)) ' global signal outlier trials flagged']);
                        display([num2str(length(combined_outliers)) ' total outlier trials excluded']);
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    

                    
                    %assign conditions to train/test classifier on
                    condensed_regs_of_interest = condensed_regs_all;
                    
                    
%                     if flags.equate_number_of_trials_in_cond_1_and_2 == 1
%                         
%                         DK1_mvpaEquateConds;
%                     end
                    
                    
                    % initialize regressors object
                    subj = init_object(subj,'regressors','conds');
                    subj = set_mat(subj,'regressors','conds',condensed_regs_of_interest);
                    subj = set_objfield(subj,'regressors','conds','condnames',S.condnames);
                    
                    % add new condensed activation pattern
                    subj = duplicate_object(subj,'pattern','spiral_d_hp_z','spiral_d_hp_z_condensed');
                    subj = set_mat(subj,'pattern','spiral_d_hp_z_condensed',temporally_condensed_data,'ignore_diff_size',true);
                    
                    zhist = sprintf('Pattern ''%s'' created by JR custom code','spiral_d_hp_z_condensed');
                    subj = add_history(subj,'pattern','spiral_d_hp_z_condensed',zhist,true);
                    
                    % clean up workspace to save RAM
                    subj = remove_mat(subj,'pattern','spiral_d_hp_z');
                    
                    % update run vector to condensed format
                    subj.selectors{1}.mat = condensed_runs;
                    subj.selectors{1}.matsize = size(condensed_runs);
                    
                    % "activate" only those trials of interest (from regs_of_interest) before creating cross-validation indices
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
                    
                    trainTestOneIter_1 = subj.selectors{4}.mat;
                    
                    subj = initset_object(subj, 'selector', 'trainTestOneIter_1', trainTestOneIter_1, 'group_name', 'trainTestOneIter');
                    
                    if flags.perform_second_round_of_zscoring == 1
                        
                        display('Performing second round of z-scoring')
                        
                        if strcmp(S.thisTrain, S.thisTest)
                            subj.patterns{end}.mat(:,active_trials) = zscore(subj.patterns{end}.mat(:,active_trials));
                        else
                            subj.patterns{end}.mat(:,intersect(active_trials, find(trainTrials))) = zscore(subj.patterns{end}.mat(:,intersect(active_trials, find(trainTrials)))')';
                            subj.patterns{end}.mat(:,intersect(active_trials, find(testTrials))) = zscore(subj.patterns{end}.mat(:,intersect(active_trials, find(testTrials)))')';
                        end
                    end
                    
                    
                    % run feature selection ANOVA: specify pvalue (if desired)
                    
                    statmap_arg.use_mvpa_ver = true; % use mvpa anova code instead of matlab's (for speed reasons)
                    
                    classifier_mask = subj.masks{1}.name; % use original mask
                    
                    if flags.anova_p_thresh ~= 1
                        subj = feature_select(subj,'spiral_d_hp_z_condensed','conds',S.thisSelector,'thresh',flags.anova_p_thresh,'statmap_arg',statmap_arg);
                        classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
                        
                    end
                    
                    % run feature selection ANOVA: specify #of voxels (if desired)
                    if S.class_args.nVox>0
                        subj = AG_feature_select_top_N_vox(subj,'spiral_d_hp_z_condensed','conds',S.thisSelector,'nVox_thresh',S.class_args.nVox,'statmap_arg',statmap_arg);
                        classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
                    end
                    %%%%%%%%
                    %FIGURE THIS PENALTY OPTIMIZATION STUFF OUT LATER
                    %[subj best_penalties penalty_iteration_results] = optimal_penalty_search(subj,'spiral_d_hp_z_condensed','conds','meta_runs_condensed_xval','meta_runs_condensed',classifier_mask,'none');
                    
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%% RUN THE CLASSIFIER (CROSS-VALIDATION)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    for p = 1:flags.num_iter_with_same_data
                        x=x+1; % increment results iteration counter
                        
                        %optimal_penalty_search(subj,patin,regsname,selgroup,runs_selname,mask,actives_selname,varargin)
                        %optimal_penalty_search(subj,'spiral_d_hp_z_condensed','conds','meta_runs_condensed_xval','classifier_mask',
                        %[subj results{x}] = optimal_penalty_search(subj,'spiral_d_hp_z_condensed','conds','meta_runs_condensed_xval',classifier_mask,S.class_args);
                        %optimal_penalty_search(subj,patin,regsname,selgroup,runs_selname,mask,actives_selname,varargin)
                        %optimal_penalty_search(subj,'spiral_d_hp_z_condensed','conds','meta_runs_condensed_xval','runs',classifier_mask,'actives_selector')
                        
                        
                        [subj results{x}] = cross_validation(subj,'spiral_d_hp_z_condensed','conds', S.thisSelector ,classifier_mask,S.class_args); %modified specifically for Enc/Ret, 081709
                        %cross_validation(subj,patin,regsname,selgroup,maskgroup,class_args,varargin)
                        
                        
                        
                        %qqq.task{tT}.subj{b}.subjStruct = subj;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        % do some important RAM clean-up and data archiving
                        %for y = 1:length(S.meta_runs)
                        if flags.generate_importance_maps == 1
                            %results_IW{x}.iterations(y).scratchpad.net.IW{1} = results{x}.iterations(y).scratchpad.net.IW{1}; % save weights to pass to JR_interpret_weights
                            for rif = 1:length(results{1}.iterations);
                                results_IW{rif}.iterations(1).scratchpad.net.IW{1} = results{1}.iterations(rif).scratchpad.weights';
                            end
                            
                        end
                        %results{x}.iterations(y).scratchpad.net.inputs{1}.exampleInput=[]; % delete huge data object from results scratchpad to free up RAM
                        %end
                        
                        
                        
                        qqq.test{CI}.penalty(penal).FeatureSel(FS).subj{b}.iter{n} = results{1};
                        qqq.subjArray = subj_array;
                        qqq.test{CI}.penalty(penal).FeatureSel(FS).subj{b}.S = S;
                        
                        
                    end
                end
                
                
                save (fullfile(S.group_mvpa_dir, S.saveName), 'qqq');
                
                
                
                if flags.generate_importance_maps == 1;
                    AG_generate_importance_maps(subj, results, results_IW, S)
                    %Create mean importance maps
                end
                
                
                clear data_log acc_percentiles
                
                time2finish = toc/60;
                display(['Finished ' S.subj_id ' in ' num2str(time2finish) ' minutes']);
            end
            
            
            if flags.generate_importance_maps == 1;
                AG_mean_impmaps(S)
            end
        end
    end
end