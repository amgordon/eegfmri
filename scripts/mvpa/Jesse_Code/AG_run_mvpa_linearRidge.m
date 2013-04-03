function [qqq, results]= AG_run_mvpa_linearRidge(subj_array)

% example usage:  JR_run_mvpa_general([1 2 4 5 8 12 13 16], 'Hits', 'CRs')

S.saveName = input('What do you want to name the results structure?', 's');
qqq.notes = input('notes about this classification', 's');
%S.trainTasks = {'Enc'  'Ret' 'respSel' 'Enc' 'respSel'};
%S.testTasks = {'Enc' 'Ret' 'respSel' 'Ret' 'Ret'};

S.subj_array = subj_array;
S.trainTasks = {'Ret'};
S.testTasks = {'Ret'};

%S.Ps = [0 10 100 1000 10000];

S.taskIdx = {'Enc', 'Ret', 'respSel', 'loc'};

S.TR_weights_idx = {1 1 1 1};  %{[0 0 0 .5 .5   ], [0 0  .5 .5 ], [0 .5 .5 0 ] ,[0 0 0 0]}; 

S.onsetsIdxSet = {[1 2]};

S.condnames = {'RT'};

%S.iterPenalty = [10 500 1000 5000 10000 50000];

S.class_args.nVoxSet = [0]; 

for tT = 1:length(S.trainTasks);
for sOI = 1:length(S.onsetsIdxSet)

    for nV = 1:length(S.class_args.nVoxSet);
        
    S.class_args.nVox = S.class_args.nVoxSet(nV);
    
    S.onsetsIdx = S.onsetsIdxSet{sOI};
    
    S.thisTrain = S.trainTasks{tT};
    S.thisTest = S.testTasks{tT};
    
    
    S.RVIdx{1} = [210 210 210 210 210 210];
    S.RVIdx{2} = [140 140 140 140 140 140];
    S.RVIdx{3} = [120 120];
    S.RVIdx{4} = [147];
    
    S.mRVIdx{1} = [210 210 210 210 210 210];
    S.mRVIdx{2} = [140 140 140 140 140 140];
    S.mRVIdx{3} = [40 40 40 40 40 40];
    S.mRVIdx{4} = [34 34 34 45];    
    
    
    
    S.idxThisTrain = find(strcmp(S.taskIdx, S.trainTasks{tT}));
    S.idxThisTest = find(strcmp(S.taskIdx, S.testTasks{tT}));
    

    %%%%%%% specify user-defined variables
    %%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for b=1:length(subj_array)  % allows subject ID #'s to be specified in short format
        tic
        
        S.subj_id = subj_array{b};
       
        if strcmp('ag1_071409', S.subj_id)
            %special fix for this subject, who doesn't have fMRI data for
            %session4
           S.RVIdx{1} = [210 210 210 210 210];
           S.mRVIdx{1} = [210 210 210 210 210];
        else
           S.RVIdx{1} = [210 210 210 210 210 210];
           S.mRVIdx{1} = [210 210 210 210 210 210];
        end
        
        %s
        S.expt_dir = '/Users/alangordon/Studies/AG1/Accumulator_fMRI/AG1/fmri_data2';
        S.mvpa_dir = [S.expt_dir '/' S.subj_id '/mvpa'];
        S.univar_dir.Enc = [S.expt_dir '/' S.subj_id '/encAnalysis'];
        S.univar_dir.Ret = [S.expt_dir '/' S.subj_id '/retAnalysis'];
        S.univar_dir.loc = [S.expt_dir '/' S.subj_id '/locAnalysis'];
        S.univar_dir.respSel = [S.expt_dir '/' S.subj_id '/respSelAnalysis'];
        S.univar_dir.BS = [S.expt_dir '/' S.subj_id '/CorGuessesBetaSeries'];
        S.importance_maps_dir=[S.expt_dir '/mvpa_results/ImpMaps_Ret_by_RT_FS' date '/' [S.thisTrain S.thisTest] ];
        S.impType = {'RT_betas'};
        S.xls_results_data_logs_dir=[S.expt_dir '/mvpa_results/xls_results_data_logs/' S.condnames{1}];
        S.group_mvpa_dir = [S.expt_dir '/mvpa_files'];
        S.rtReg_dir = [S.expt_dir '/' S.subj_id '/RTReg'];
        
        % create these directories if they don't already exist
        if ~exist(S.importance_maps_dir,'dir')
            mkdir(S.importance_maps_dir);
        end
        if ~exist(S.xls_results_data_logs_dir, 'dir')
            mkdir(S.xls_results_data_logs_dir)
        end
        
        if S.idxThisTrain==S.idxThisTest
            S.runs_vector = S.RVIdx{S.idxThisTrain};
            S.meta_runs = S.mRVIdx{S.idxThisTrain};
        else
            S.runs_vector = horzcat(S.RVIdx{S.idxThisTrain}, S.RVIdx{S.idxThisTest});
            S.meta_runs = [sum(S.RVIdx{S.idxThisTrain}), sum(S.RVIdx{S.idxThisTest})];
        end
        

        S.num_runs = length(S.runs_vector);
        S.num_vols = sum(S.runs_vector);
        
        %1 for unsmoothed, 0 for smoothed;
        S.use_unsmoothed = 1;
        S.smoothTxt = {'smoothed' 'unsmoothed'};
        
        % load some .mat files into memory
        S.vol_info = load([S.expt_dir '/vol_info.mat']); %get functional data resolution info for spm .img writing
        
        [S.onsets_unsorted, S.img_files] = AG_mvpa_onsets_and_images_RTReg(S);
       
        idx.outOfTimeBounds =  find(ismember(S.onsets_unsorted{1}, [0 280  560 840 1120 1400]));
        
        S.onsets_unsorted{1}(idx.outOfTimeBounds) = []; %remove onsets corresponding to trials that we don't have relevant TRs for.
        %S.onsets_unsorted{1} = S.onsets_unsorted{1}(S.onsets_unsorted{1}>9);
        
        [S.onsets{1}, ix_onsets_sort] = sort(S.onsets_unsorted{1});
        qqq.onsets{b} = S.onsets;
        
        S.num_conds = size(S.onsets,2);
        
        
        % OPTIONAL:  specify previously saved mvpa workspace to bypass time-consuming data extraction and preprocessing
        
        %mvpa_workspace = [S.expt_dir '/' S.subj_id '/mvpa/' S.subj_id '_666' S.roi_name '.mat'];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set flags (% unless otherwise indicated: 1 = yes please, 0 = no thanks)
        flags.use_premade_workspace = 1;
        flags.num_results_iter = 1; % number of times to run the entire classification process (select subset of the data and train/test classifier)
        flags.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data
        flags.equate_number_of_old_new_trials_per_subjective_bin = 0; % equate number of trials per subjective bin
        flags.equate_number_of_trials_in_cond_1_and_2 = 0; % equate number of trials in conditions 1 and 2 (RECOMMENDED)
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
        flags.mask = 19; %1 = 'merged_AAL_ROIs_FIXED_HOLES.nii'; 2 = subject-specific analysis mask. 3 = noLPFC jesse mask. 4 = AlanMaskNoLPFCNoMotor.img  5 = AlanMaskNoIPSNoLPFCNoMotor.img 
        flags.useIdealParams = 0;
        flags.betaSeriesAnalysis = 0;
        %anova_nVox_thresh = 100;
        
        
        S.exp_name = 'Acc1';    
        
        
        if flags.mask==1
            S.roi_name = 'merged_AAL_ROIs_FIXED_HOLES.nii';
            S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/fmri_data/mvpa_files/' S.roi_name ];
        elseif flags.mask==2
            S.roi_name = 'mask.img';
            S.roi_file = [S.univar_dir.(S.thisTrain) '/' S.roi_name];
        elseif flags.mask==3
            S.roi_name = 'LPFC_improved.img';
            S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/mvpa/maskGeneration/' S.roi_name ];
        elseif flags.mask==4
            S.roi_name = 'AlanMaskNoLPFCNoMotor.img';
            S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/mvpa/newMask/' S.roi_name ];
        elseif flags.mask==5
            S.roi_name = 'noLPFCNoMotorNoIPS.img';
            S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/mvpa/noLPFCnoMotorNoIPS/' S.roi_name ];
        elseif flags.mask==6
            S.roi_name = 'rIPS_sphere.img';
            S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/fmri_data/group_analyses/ADCorGuessesNoIPS/CorByGreaterEvidence_guessCor_vsFix/SFN09ROIs/' S.roi_name ];
        elseif flags.mask==7
            S.roi_name = 'noMotornoLPFCnoLIPyesTempPoleFixed.img';
            S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/mvpa/noLPFCnoMotorNoIPSWithTemporalPoles/' S.roi_name];
        elseif flags.mask==8
            S.roi_name = 'EncEncImpMapBoth0005.img';
            S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/mvpa/impMapMask/' S.roi_name];
        elseif flags.mask==9
            S.roi_name = 'tsrmask23.img';
            S.roi_file = ['/Users/alangordon/Studies/AG1/Accumulator_fMRI/AG1/mvpa/mask_012510/' S.roi_name];
        elseif flags.mask==10
            S.roi_name = 'mask23.img';
            S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/mvpa/mask_012510/' S.roi_name];
        elseif flags.mask==11
            S.roi_name = 'Conjp0005tsrmask23.img';
            S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/mvpa/ContrastDefinedMask/' S.roi_name];
        elseif flags.mask==12
            S.roi_name = 'mask.img';
            S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/fmri_data2/ag1_080509/respSelAnalysis' S.roi_name];
        elseif flags.mask==13
            S.roi_name = 'RVsL_BothTails_p0001.img';
            S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/mvpa/RSEvMask/' S.roi_name];
        elseif flags.mask==14
            S.roi_name = 'mvpa_mask.img';
            S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/fmri_data2/ag1_021509/RTReg/' S.roi_name];
        elseif flags.mask==15
            S.roi_name = 'TAbsValEncRetImportanceMaps.img';
            S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/Modeling/RT_and_Evidence_Analysis/' S.roi_name];
        elseif flags.mask==16
            S.roi_name = 'top1000vox_RetByRT.img';
            S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/Modeling/RT_and_Evidence_Analysis/' S.roi_name];    
        elseif flags.mask==17
            S.roi_name = 'ImpMapFaceSceneGreaterThan0005.img';
            S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/fmri_data2/mvpa_results/ImpMaps_08-Jul-2010/EncEnc/pos/' S.roi_name];
        elseif flags.mask==18
            S.roi_name = ['thresh1000' S.subj_id '.img'];
            S.roi_file = ['/Users/alangordon/Studies/AG1/Accumulator_fMRI/AG1/fmri_data2/mvpa_results/ImpMaps_08-Jul-2010/EncEnc/pos/' S.roi_name];
        elseif flags.mask==19
            S.roi_name = [S.subj_id '_posImpThresh1000.img'];
            S.roi_file = ['/Users/alangordon/Studies/AG1/Accumulator_fMRI/AG1/fmri_data2/mvpa_results/ImpMaps_19-Aug-2010/EncRet/pos/' S.roi_name];
        elseif flags.mask==20
            S.roi_name = [S.subj_id '_Occ_ImpThresh1000.img'];
            S.roi_file = ['/Users/alangordon/Studies/AG1/Accumulator_fMRI/AG1/mvpa/occipital_mask/' S.roi_name];        
        end
        
        
        
        S.TRs_to_average_over_train = 1:length(S.TR_weights_idx{S.idxThisTrain}); %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier
        S.TR_weights_train = S.TR_weights_idx{S.idxThisTrain}; % should sum to 1
        
        S.TRs_to_average_over_test = 1:length(S.TR_weights_idx{S.idxThisTest}); %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier
        S.TR_weights_test = S.TR_weights_idx{S.idxThisTest}; % should sum to 1
        %S.TR_weights = [.0072 .2168 .3781 .2742 .1237];  % from SPM canonical values at 1,3,5,7,and 9 sec post-stimulus
        
%         S.TRs_to_average_over_train = [-2 -1];
%         S.TRs_to_average_over_test = [-2 -1];
%         
%         S.TR_weights_train = [.5 .5];
%         S.TR_weights_test = [.5 .5];
%         
        
        % classifier parameters
        S.class_args.train_funct_name = 'AG_train_ridge';
        S.class_args.test_funct_name = 'AG_test_ridge';
        S.class_args.nHidden = 0;
        S.class_args.perfmet_functs =  'perfmet_xcorr';
        
        
        if flags.useIdealParams
            S.class_args.penalty = idealParams.penalty.(S.thisTrain).(S.thisTest);
            %S.class_args.nVox = idealParams.nVox.(S.thisTrain).(S.thisTest);
        else
            S.class_args.penalty = 1000;
            %S.class_args.nVox = 1000;  %S.Ps(pP); % 0 = no feature selection
        end
        
                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S.workspace = fullfile(S.mvpa_dir, [S.subj_id '_linReg_train' S.thisTrain '_test' S.thisTest '_' S.roi_name '_' S.smoothTxt{S.use_unsmoothed + 1} '.mat']);

        existWorkspace = exist(S.workspace);
        
        Sequal = 0;
%         if existWorkspace
%             WorkCheck = load(S.workspace);
%             Sequal = isequal(struct2cell(WorkCheck.S), struct2cell(S));
%         end
        
        
        if (flags.use_premade_workspace&&existWorkspace) %(flags.use_premade_workspace && Sequal) %
            load(S.workspace, 'subj');
        else
            if flags.betaSeriesAnalysis == 1;
                [subj]= AG_mvpa_load_betaSeries;
            else

                [subj] = AG_mvpa_load_and_preprocess_raw_data(S);
            end
        end
        
        S.preprocPatName = subj.patterns{end}.name;
        S.thisMaskName =  subj.masks{end}.name;
        S.preprocPatNameCondensed = [subj.patterns{end}.name '_condensed'];
        
        if strcmp(S.thisTrain, S.thisTest)
            S.thisSelector = 'meta_runs_condensed_xval';
        else
            S.thisSelector = 'trainTestOneIter';
        end
        
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
            
            idx_condense =find(all_regs);
            
            
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
            
            if strcmp(S.thisTrain, S.thisTest)
                meta_runs_train = find(all_trials);
                meta_runs_test = [];
            else
                            %applies only when train and test data are different
                meta_runs_train = idx_condense(find(meta_runs_condensed==1));
                meta_runs_test = idx_condense(find(meta_runs_condensed==2));
            end
            

            

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % select TRs of interest (to correspond with peak post-stim BOLD response)
                        
            
            data_by_TR_train = [];
            for dt = 1:length(S.TR_weights_idx{S.idxThisTrain})

                data_by_TR_train(dt,:,:) = S.TR_weights_train(dt)*subj.patterns{end}.mat(:, abs(meta_runs_train+(dt-3))); %change back to dt_pos - 1 Hack alert! remove abs
            end
            
            temporally_condensed_data_train = squeeze(sum(data_by_TR_train(S.TRs_to_average_over_train,:,:),1));
            
            clear data_by_TR_train
            


            data_by_TR_test = [];
            for dt = 1:length(S.TR_weights_idx{S.idxThisTest})
                data_by_TR_test(dt,:,:) = S.TR_weights_test(dt)*subj.patterns{end}.mat(:, abs(meta_runs_test+(dt-3))); %change back to dt_pos - 1 %Hack alert! remove abs
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
            
            % Artificially balance the number of trials given each subjective memory response
            if flags.equate_number_of_old_new_trials_per_subjective_bin == 1
                for j = 1:2
                    Face_trials = find(condensed_regs_all(j,:));
                    Scene_trials = find(condensed_regs_all(j+5,:));
                    
                    num_Face = length(OLD_trials);
                    num_Scene = length(NEW_trials);
                    
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
            
%             face = condensed_regs_all(1,:);
%             scene = condensed_regs_all(2,:);
%             
%             left = condensed_regs_all(1,:);
%             right = condensed_regs_all(2,:);
            
  
            
            %assign conditions to train/test classifier on
            condensed_regs_of_interest = [];
            trainTestOneIter_1 = meta_runs_condensed;
            %eval(['condensed_regs_of_interest(1,:) = ' S.condnames{1} ';'])
            %eval(['condensed_regs_of_interest(2,:) = ' S.condnames{2} ';'])
            
            mvpa_ons = load([S.rtReg_dir, '/mvpa_ons']);
            condensed_regs_of_interest_h = [mvpa_ons.RTs{S.onsetsIdx}]; %make this dynamic
            
            condensed_regs_of_interest_h(idx.outOfTimeBounds) = []; %remove RTs corresponding to trials that we don't have relevant TRs for.
            
            condensed_regs_of_interest = condensed_regs_of_interest_h(ix_onsets_sort);
            
            
            
            S.RTs = mvpa_ons.RTs;
            
            S.theseRTs = condensed_regs_of_interest;
            
            if flags.equate_number_of_trials_in_cond_1_and_2 == 1
                
                AG1_mvpaEquateConds;
            end
            
%             display([num2str(count(condensed_regs_of_interest(1,:)==1)) '
%             trials in condition ' S.condnames{1}])
%             display([num2str(count(condensed_regs_of_interest(2,:)==1)) ' trials in condition ' S.condnames{2}])
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subj = init_object(subj,'selector','meta_runs_condensed');
            subj = set_mat(subj,'selector','meta_runs_condensed', meta_runs_condensed);
            subj = create_xvalid_indices(subj,'meta_runs_condensed');
            
            % initialize regressors object
            subj = init_object(subj,'regressors','conds');
            subj = set_mat(subj,'regressors','conds',condensed_regs_of_interest);
            subj = set_objfield(subj,'regressors','conds','condnames',S.condnames);
            
            % add new condensed activation pattern
            subj = duplicate_object(subj,'pattern', S.preprocPatName, S.preprocPatNameCondensed);
            subj = set_mat(subj,'pattern', S.preprocPatNameCondensed,temporally_condensed_data,'ignore_diff_size',true);
            
            zhist = sprintf('Pattern ''%s'' created by JR custom code', S.preprocPatNameCondensed);
            subj = add_history(subj,'pattern', S.preprocPatNameCondensed,zhist,true);
            
            % clean up workspace to save RAM
            subj = remove_mat(subj,'pattern', S.preprocPatName);
            
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
            
            %trainTestOneIter_1 = subj.selectors{4}.mat;
            
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
            
            %statmap_arg.use_mvpa_ver = true; % use mvpa anova code instead of matlab's (for speed reasons)
            statmap_funct = 'AG_statmap_random';
            
            classifier_mask = subj.masks{1}.name; % use original mask
            
            if flags.anova_p_thresh ~= 1
                subj = feature_select(subj, S.preprocPatNameCondensed,'conds',S.thisSelector,'thresh',flags.anova_p_thresh,'statmap_arg',statmap_arg);
                classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
                
            end
            
            % run feature selection ANOVA: specify #of voxels (if desired)
            if S.class_args.nVox>0
                subj = AG_feature_select_top_N_vox(subj, S.preprocPatNameCondensed,'conds',S.thisSelector,'nVox_thresh',S.class_args.nVox, 'statmap_funct', statmap_funct);
                classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
            end
            %%%%%%%%

           
            
            %FIGURE THIS PENALTY OPTIMIZATION STUFF OUT LATER
            %[subj best_penalties penalty_iteration_results] = optimal_penalty_search(subj,'spiral_d_hp_z_condensed','conds','meta_runs_condensed_xval','meta_runs_condensed',classifier_mask,'none');
          
            S.meanPats = mean(subj.patterns{5}.mat,1); %Just a test
            
            
            %%%%%%%%%%%%%%%%%%%%%% RUN THE CLASSIFIER (CROSS-VALIDATION)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for p = 1:flags.num_iter_with_same_data
                x=x+1; % increment results iteration counter
                
                %optimal_penalty_search(subj,patin,regsname,selgroup,runs_selname,mask,actives_selname,varargin)
                %optimal_penalty_search(subj,'spiral_d_hp_z_condensed','conds','meta_runs_condensed_xval','classifier_mask',  
                %[subj results{x}] = optimal_penalty_search(subj,'spiral_d_hp_z_condensed','conds','meta_runs_condensed_xval',classifier_mask,S.class_args);
                %optimal_penalty_search(subj,patin,regsname,selgroup,runs_selname,mask,actives_selname,varargin)
                %optimal_penalty_search(subj,'spiral_d_hp_z_condensed','conds','meta_runs_condensed_xval','runs',classifier_mask,'actives_selector')
                
                S.class_args.use_matlab = true;
                S.class_args.no_constant = 0;
                
                [subj results{x}] = cross_validation(subj, S.preprocPatNameCondensed,'conds', S.thisSelector ,classifier_mask,S.class_args, 'perfmet_functs', S.class_args.perfmet_functs); %modified specifically for Enc/Ret, 081709
                                    %cross_validation(subj,patin,regsname,selgroup,maskgroup,class_args,varargin)

                
                
                %qqq.task{tT}.subj{b}.subjStruct = subj;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % do some important RAM clean-up and data archiving
                %for y = 1:length(S.meta_runs)
                    if flags.generate_importance_maps == 1
                        %results_IW{x}.iterations(y).scratchpad.net.IW{1} = results{x}.iterations(y).scratchpad.net.IW{1}; % save weights to pass to JR_interpret_weights
                        for rif = 1:length(results{1}.iterations);
                        results_IW{rif}.iterations(1).scratchpad.net.IW{1} = results{1}.iterations(rif).scratchpad.ridge.betas;
                        end
                    
                    end
                    %results{x}.iterations(y).scratchpad.net.inputs{1}.exampleInput=[]; % delete huge data object from results scratchpad to free up RAM
                %end
                
                
%                save(fullfile(S.mvpa_dir, 'abs_acts_diff_ret'), 'abs_acts_diff_ret');
             
                %results{n}.iterations.scratchpad.W = [];
                %results{n}.iterations.scratchpad.weights = [];

                qqq.test{nV}.subj{b}.iter{n} = results{n};
                qqq.subjArray = subj_array;
                qqq.test{nV}.subj{b}.S = S;
                

%                qqq.actsDiff{tT}.subj{b} = acts_diff_vector;
                
               
                
            end
        end
        
        save (fullfile(S.group_mvpa_dir, S.saveName), 'qqq');
        
        
        
        
        if flags.generate_importance_maps == 1;
            AG_generate_importance_maps_linreg(subj, results, results_IW, S)
            
            
            %Create mean importance maps
            
            
            
        end
        %             if flags.anova_p_thresh == 1 % NO ANOVA VERSION
        
        clear data_log acc_percentiles
        
        time2finish = toc/60;
        display(['Finished ' S.subj_id ' in ' num2str(time2finish) ' minutes']);
    end
    end
end
end

if flags.generate_importance_maps == 1;
    AG_mean_impmaps(S)  
end
