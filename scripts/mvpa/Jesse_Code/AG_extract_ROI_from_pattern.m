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

S.TR_weights_idx = {[0 0 0 .5 .5   ], [0 0  .5 .5 ], [0 .5 .5 0 ] ,[0 0 0 0]}; 

S.onsetsIdxSet = {[1 2]};

S.condnames = {'RT'};

%S.iterPenalty = [10 500 1000 5000 10000 50000];

for tT = 1:length(S.trainTasks);
    for sOI = 1:length(S.onsetsIdxSet)
        
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
            S.expt_dir = '/Users/alangordon/Accumulator_fMRI/AG1/fmri_data2';
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
            S.actsToCorrelateDir = [S.expt_dir '/mvpa_files/LinReg18Subs_FaceSceneCombined'];
            
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
            flags.mask = 6; %1 = 'merged_AAL_ROIs_FIXED_HOLES.nii'; 2 = subject-specific analysis mask. 3 = noLPFC jesse mask. 4 = AlanMaskNoLPFCNoMotor.img  5 = AlanMaskNoIPSNoLPFCNoMotor.img
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
                S.roi_name = 'rBen_IPS.nii';
                S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/fmri_data/group_analyses/ADCorGuessesNoIPS/CorByGreaterEvidence_guessCor_vsFix/SFN09ROIs/' S.roi_name ];
            elseif flags.mask==7
                S.roi_name = 'noMotornoLPFCnoLIPyesTempPoleFixed.img';
                S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/mvpa/noLPFCnoMotorNoIPSWithTemporalPoles/' S.roi_name];
            elseif flags.mask==8
                S.roi_name = 'EncEncImpMapBoth0005.img';
                S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/mvpa/impMapMask/' S.roi_name];
            elseif flags.mask==9
                S.roi_name = 'tsrmask23.img';
                S.roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/mvpa/mask_012510/' S.roi_name];
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
            end
            
            
            
            S.TRs_to_average_over_train = 1:length(S.TR_weights_idx{S.idxThisTrain}); %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier
            S.TR_weights_train = S.TR_weights_idx{S.idxThisTrain}; % should sum to 1
            
            S.TRs_to_average_over_test = 1:length(S.TR_weights_idx{S.idxThisTest}); %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier
            S.TR_weights_test = S.TR_weights_idx{S.idxThisTest}; % should sum to 1
            %S.TR_weights = [.0072 .2168 .3781 .2742 .1237];  % from SPM canonical values at 1,3,5,7,and 9 sec post-stimulus
            
            % classifier parameters
            S.class_args.train_funct_name = 'AG_train_ridge';
            S.class_args.test_funct_name = 'AG_test_ridge';
            S.class_args.nHidden = 0;
            S.class_args.perfmet_functs =  'perfmet_xcorr';
            
            
            if flags.useIdealParams
                S.class_args.penalty = idealParams.penalty.(S.thisTrain).(S.thisTest);
                S.class_args.nVox = idealParams.nVox.(S.thisTrain).(S.thisTest);
            else
                S.class_args.penalty = 10;
                S.class_args.nVox = 0;  %S.Ps(pP); % 0 = no feature selection
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
            
            mean_vals = mean(subj.patterns{end}.mat(:,S.onsets{1}/2+1));
            
            rawActs = load(S.actsToCorrelateDir);
            
            ActsVec = [rawActs.qqq.test{1}.subj{b}.iter{:}.iterations.acts];
            
            [qqq.r(b), qqq.p(b)] = corr(mean_vals', ActsVec');
        end
    end
end