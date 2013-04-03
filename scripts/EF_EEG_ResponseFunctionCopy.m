function [S results] = EF_EEG_ResponseFunction(subj_id)
%to do:

%% Immediately

% look at IPS regions instead of AnG.

% replicate classification of motor cortex, with onset locked instead of
% RT.

% 1) classify memory states with EEG.  2) read out trialwise certainty.
% 3) use this certainty as a parametric modulator.  

%       -Can we do it without controlling for subjective response? (I HOPE
%       SO)
%       -Can we do it controlling for subjective response? (more
%       interesting, less likely...)
%       -Can we do each of these by binning the erp classifications?

% why can't we classify BOLD in some subjects e.g #2, #3?

% Can we really not classify parietal BOLD activity with parietal
% electrodes?  What parameters can we tweak to make this so?



% 1) classify memory states with EEG.  2) classify memory states with BOLD
% 3) is trialwise classifier certainty with these two modalities related?

%still have interp1 errors

%also, is using premade workspaces really worth it??

% ef_092711 sucks at predicting memory with EEG. Why? (Because it only has
% two crappy runs.)

%% Eventually

% put dynamic ability to modulate which diagnostics to perform in Params.

% standardize results structure at both the group and individual level.  

% use mixed models in R to look at accuracy

% use non-normalized, non-smoothed patterns for deconvolution

% direct correlation of eeg activity and BOLD? (e.g. what anthony wants).

% test out ridge.  can we include regressors of no interest (e.g. model
% each specific button?


%%
results = [];

%generate params
[S par idx] = EF_EEGResponseParamsCopy(subj_id);

%% Begin Analysis Here
if S.loadBOLDData
    existWorkspace = exist(S.workspace);
    
    if strcmp(S.lock, 'RT')
        S.onsets = {idx.allTrials + idx.RT};
    else
        S.onsets = {idx.allTrials};
    end
    
    S.num_conds = size(S.onsets,2);
    
    if (S.use_premade_workspace&&existWorkspace) % if we are supposed to use premade workspace, and one with the correct name exists
        load(S.workspace, 'subj','condensed_regs_of_interest', 'condensed_runs');
    else
        
        if S.useBetaMapPatterns % read beta map patterns
            
            % initialize subj structure
            subj = init_subj(S.exp_name,S.subj_id);
            
            % load mask file
            subj = load_spm_mask(subj,S.roi_name,S.roi_file);
            
            for f=1:length(S.img_files)
                
                % load betaMaps
                subj = load_analyze_pattern(subj,S.preprocPatCondensedName{f},S.roi_name, S.img_files{f},'single',true);
                loadedPat = get_mat(subj, 'pattern', S.preprocPatCondensedName{f});
                
                thisPat = single(nan(size(loadedPat,1), size(idx.TRsExist,1)));
                
                thisPat(:,idx.TRsExist(:,f)) = loadedPat;
                subj = set_mat(subj, 'pattern', S.preprocPatCondensedName{f}, thisPat);
            end
            
            % make runs vector
            subj = init_object(subj,'selector','runs');
            
            subj = set_mat(subj,'selector','runs',idx.sess');
            
            clear thisPat loadedPat
            
            save (S.workspace)
            
        else %load raw data patterns
            
            [subj] = EF_mvpa_load_and_preprocess_raw_data(S);
            
            % turn onsets into regressors matrix
            all_regs = zeros(S.num_conds,S.num_vols); % initialize regs matrix as conditions x timepoints
            for cond = 1: S.num_conds
                for trial = 1: length(S.onsets{cond})
                    time_idx = round(S.onsets{cond}(trial)/S.TR) + 1; % divide by 2 and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
                    all_regs(cond, round(time_idx)) = 1;
                end
            end
            
            % condense regs by removing zeros
            condensed_runs = [];
            condensed_regs_of_interest = [];
            trial_counter = 1;
            for i = 1: size(all_regs,2)
                if ~isempty(find(all_regs(:,i))) % if not a rest timepoint
                    condensed_regs_of_interest(:,trial_counter) = all_regs(:,i);
                    condensed_runs(1,trial_counter) = subj.selectors{1}.mat(i);
                    trial_counter = trial_counter + 1;
                end
            end
            
            % create meta runs structure
            idx_condense =find(sum(all_regs));
            all_trials = sum(all_regs,1);
            meta_runs = find(all_trials);
            
            %% sessions
            regs = condensed_regs_of_interest(1,:);
            
            for i = 1:length(regs)
                allSess_h(condensed_runs(i),i) = 1;
            end
            
            allSess_h(end,:) = 1;
            allSess = allSess_h;
            
            %%
            if S.condensePatterns
                preprocPat = get_mat(subj, 'pattern', S.preprocPatName);
                for tr = 1:length(S.TR_weights)
                    
                    TR_weights = S.TR_weights{tr};
                    if S.interp
                        temporally_condensed_data = interp1(preprocPat',find(TR_weights)-1+S.onsets{1}/S.TR)';
                    else
                        %condense data
                        data_by_TR = [];
                        for dt = 1:length(S.TR_weights)
                            data_by_TR(dt,:,:) = S.TR_weights(dt)*preprocPat(:,meta_runs+(dt-1));
                        end
                        temporally_condensed_data = squeeze(sum(data_by_TR(TRs_to_average_over,:,:),1));
                    end
                    
                    thisPat = [S.preprocPatCondensedName '_TR' num2str(tr)];
                    
                    % add new condensed activation pattern
                    subj = duplicate_object(subj,'pattern',S.preprocPatName,thisPat);
                    subj = set_mat(subj,'pattern',thisPat,temporally_condensed_data,'ignore_diff_size',true);
                    subj = set_objfield(subj, 'pattern', thisPat, 'group_name', S.TRPatName);
                    
                    clear temporally_condensed_data;
                    clear data_by_TR;
                end
                
                subj = remove_mat(subj,'pattern', S.preprocPatName);
                subj = remove_mat(subj,'pattern','spiral');
                clear preprocPat;
                save (S.workspace)
            end
        end
    end
    
    
    
    %% post workspace loading
    
    
    % mask a workspace mask with another mask.
    if ~isempty(S.secondaryMask)
        subj = load_spm_mask(subj, 'secondaryMask', S.secondaryMask);
        subj = intersect_masks(subj,S.roi_name,'secondaryMask');
        for f = 1:length(S.img_files)
            subj = create_pattern_from_mask(subj, S.preprocPatCondensedName{f}, subj.masks{end}.name , S.preprocPatCondensedMaskedName{f});
            subj = set_objfield(subj, 'pattern', S.preprocPatCondensedMaskedName{f}, 'group_name', S.patFinalGroupName);
        end
    end
    
    allPats = get_group_as_matrix(subj,'pattern',S.patFinalGroupName);
    
    
    %interpolate between the TR-space beta values, to get to the
    %onset-space beta values.
    szPats = size(allPats);
    
    warning off all
    if S.interp
        interpolatedPats = nan(szPats(1)-1, szPats(2), szPats(3));
        for t=1:size(allPats,3)
            interp_offset = 1+(idx.allTrials(t)-floor(idx.allTrials(t)))/S.TR;
            interp_vals = interp_offset:1:((S.num_FIR_bins-2) + interp_offset);
            interpolatedPats(:,:,t) = interp1(squeeze(allPats(:,:,t)),interp_vals);
        end
        allPats = interpolatedPats;
    end
    warning on all
    
    idx.trialHasEstimateTRWise = ~isnan(squeeze(mean(allPats,2)))';
    idx.trialHasEstimate = all(idx.trialHasEstimateTRWise');
end

%% load eeg data
if S.loadEEGData
    eeg = load(fullfile(S.eeg_dir, S.eegDataType));
    
    if strcmp(S.eegDataType, 'spectraldata3_RT.mat')
        eegDat_h = squeeze(eeg.data.normpower);
    else
        eegDat_h = eeg.data.trialdata;
    end
    
    %% define bad channels and trials
    channel_rois;
    S.chnls = chnls;
    clear chnls;
    
    eegfmri_ON_trialch_info
    
    S.goodChannels = true(eeg.data.Nchan,1);
    S.goodChannels(S.badch{par.subNo}) = false;
    
    %% channel and timepoint selection
    ChOI_h = cell(size(S.ChOI));
    for c = 1:length(S.ChOI)
        thisChan = ['S.chnls.' S.ChOI{c}];
        ChOI_h{c} = eval(thisChan);
    end
    ChOI = [ChOI_h{:}];
    chansToInclude = sort(setdiff(ChOI, S.badch{par.subNo}));
    
    for t = 1:length(S.TOI)
        thisTOI = S.TOI{t};
        eegDat = eegDat_h(chansToInclude,:,thisTOI);
        
        %% reshape and downsample
        sz = size(eegDat);
        S.szEEGInit = sz;
        
        %
        eegDatReSam = zeros(sz(1), sz(2), ceil(sz(3)/S.decimationFactor));
        for j = 1:sz(1)
            for k = 1:sz(2)
                thisEEGTrial = squeeze(eegDat(j,k,:));
                thisResEEGTrial = single(decimate(double(thisEEGTrial), S.decimationFactor));
                eegDatReSam(j,k,:) = thisResEEGTrial;
            end
        end
        
        eegDatSf = shiftdim(eegDatReSam,2);
        eegDatResh = reshape(eegDatSf, size(eegDatSf,1)*size(eegDatSf,2), size(eegDatSf,3));
        
        eegDat_pats = eegDatResh';
        S.szEEGPats = size(eegDat_pats);
        %
        %
        %     % reshape the data such that features now combine across time and channels
        %     eegDatSf = shiftdim(eegDat,2);
        %     eegDatResh = reshape(eegDatSf, sz(1)*sz(3), sz(2));
        %
        %     % downsample the eeg data
        %     eegDatReSam = zeros(ceil(( sz(1)*sz(3))/S.decimationFactor), sz(2));
        %     for k = 1:size(eegDatResh,2)
        %         eegDatReSam(:,k) = single(decimate(double(eegDatResh(:,k)),S.decimationFactor))';
        %     end
        %
        %     eegDat_pats = eegDatReSam';
        %     S.szEEGPats = size(eegDat_pats);
        
        %% select good trials
        S.goodtr = true(size(idx.noResp));
        S.goodtr(S.badtr{par.subNo}) = false;
        
        thisModIdxEEG = (~idx.noResp & S.goodtr);
        %     %% validation test of BOLD patterns
        %     for m=1:size(allPats,1)
        %
        %         BOLDButtonPress = squeeze(mean(allPats(m,:,(~idx.noResp & idx.TRsExist(:,m)'))));
        %         BOLDFix = squeeze(mean(allPats(m,:,(idx.noResp & idx.TRsExist(:,m)'))));
        %
        %         [~, ~, ~, results.tTestBOLDButtonPressVsFix.TR(m).stats] = ttest2(BOLDButtonPress', BOLDFix);
        %     end
        %
        %     %     %can you predict memory state?
        %     for m=1:size(allPats,1)
        %         idxCorMemory = idx.cor & S.goodtr & idx.TRsExist(:,m)';
        %         sparseX = sparse(squeeze(allPats(m,:,idxCorMemory,:)))';
        %         hitVsCRs = 2*idx.respOld(idxCorMemory)'-1;
        %         results.classifyMemoryStateWithBOLD(m) = EF_x_validation(sparseX,hitVsCRs,S,S.ValidationLambda,'discrete', 'liblinear');
        %     end
        
        % Can you predict trialwise EEG with BOLD?
        for m=1:size(allPats,1)
            nanpats_h = isnan(squeeze(allPats(m,:,:)));
            idx.nanpats = sum(nanpats_h)' > 0; %index of BOLD patterns that have nans in them
            thisTrialSubset = (idx.TRsExist(:,m) & ~idx.nanpats); % which trials should be included in the classification?
            
            BOLD = squeeze(allPats(m,:,thisTrialSubset))';
            eegMean = mean(eegDat_pats(thisTrialSubset,:),2);
            
            results.EEGWithBOLD.BOLDPat(m).EEGBin(t) = EF_x_validation(BOLD,eegMean,S,S.ValidationLambda,'continuous', 'svr');
            results.EEGWithBOLD.BOLDPat(m).EEGBin(t).mod = [];
        end
        %% validation tests of EEG
        %     EEGButtonPress = squeeze(mean(mean(eegDat(:, ~idx.noResp, :),3)));
        %     EEGFix = squeeze(mean(mean(eegDat(:, idx.noResp, :),3)));
        %
        %     [~, ~, ~, results.tTestEEGButtonPressVsFix.stats] = ttest2(EEGButtonPress', EEGFix);
        %
        %     YRTs = eeg.data.RTs(thisModIdxEEG);
        %
        %     XRTs = eegDat_pats(thisModIdxEEG, :);
        %
        %     %can you predict the RT?
        %     results.classifyRTWithEEG = EF_x_validation(XRTs,YRTs,S,S.ValidationLambda,'continuous', 'svr');
        %
        %     %can you predict a button press?
        %     sparseX = sparse(eegDat_pats(S.goodtr,:));
        %     respsNegPos = 2*idx.noResp(S.goodtr)'-1;
        %     results.classifyButtonPressWithEEG = EF_x_validation(sparseX,respsNegPos,S,S.ValidationLambda,'discrete', 'liblinear');
        %
        
        %         %can you predict memory state?
        S.idxCorMemory = idx.cor & S.goodtr & ismember(idx.sess', par.scans_to_include);
        S.idxCorMemoryHCWithR = idx.cor & S.goodtr & idx.highConf & ismember(idx.sess', par.scans_to_include);
        S.idxCorMemoryHCNoR = idx.cor & S.goodtr & (~idx.recollect) & ismember(idx.sess', par.scans_to_include);
        S.idxCorMemoryHCHitsAndR = idx.cor & S.goodtr & idx.highConf & idx.old & ismember(idx.sess', par.scans_to_include);
        
        %results.classifyMemoryState = EF_EEG_classifyOldNew(S.idxCorMemory,idx.respOld,S,eegDat_pats);
        %results.classifyHCMemoryState = EF_EEG_classifyOldNew(S.idxCorMemoryHCWithR,idx.respOld,S,eegDat_pats);
        %results.classifyHCHitsVsHCCRs = EF_EEG_classifyOldNew(S.idxCorMemoryHCNoR,idx.respOld,S,eegDat_pats);
        %results.classifyHCHitsVsRecollection = EF_EEG_classifyOldNew(S.idxCorMemoryHCHitsAndR,idx.recollect,S,eegDat_pats);
        
        %% classify EEG with BOLD
        
        
        %% classify BOLD with EEG
        
        %     if S.TRByTR
        %         YMat = squeeze(mean(allPats(:,:,:),2))';
        %     else
        %         YMat = squeeze(mean(allPats(:,:,:),2));
        %     end
        %
        %     %% classify
        %     for m=1:size(YMat,2)
        %         thisModIdx = (~idx.noResp & S.goodtr & idx.TRsExist(:,m)' & idx.trialHasEstimate);
        %         Y = double(YMat(thisModIdx,m));
        %         X = eegDat_pats(thisModIdx, :);
        %         if S.normalizeFeatures
        %             X = zscore(X);
        %         end
        %         results.classifyBOLDWithEEG.TR(m) = EF_x_validation(X,Y,S,S.lambda,'continuous', 'svr');
        %     end
        
    end
end
fprintf('finished subject %s', subj_id)
end

function [theseResults] = EF_EEG_classifyOldNew(thisIdx,classLabels,S,eegDat_pats)

        fullX = eegDat_pats(thisIdx,:);
        sparseX = sparse(fullX);
        
        hisVsCRs_1And2 = 2-classLabels(thisIdx)';
        hitVsCRs_pos1AndNeg1 = 2*classLabels(thisIdx)'-1;
        
        if strcmp(S.classifier, 'glmnet')
            theseResults = EF_x_validation(fullX,hisVsCRs_1And2,S,S.ValidationLambda,'discrete', S.classifier);
        else
            theseResults = EF_x_validation(sparseX,hitVsCRs_pos1AndNeg1,S,S.ValidationLambda,'discrete', S.classifier);
        end

end


