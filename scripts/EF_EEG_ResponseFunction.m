function [S results] = EF_EEG_ResponseFunction(subj_id)
%to do:

% separate some stuff into  mini functions

% run EEG ERP creation software to zscore each channel to itself, run-wise.

% have more flexible control over which trials get selected

% test out ridge.  can we include regressors of no interest (e.g. model
% each specific button?

% include classification of memory states

% classification across time points

% what's with bold data for subject 1?  try with a
% smaller roi.

% use mixed models in R to look at accuracy

% use non-normalized, non-smoothed patterns for deconvolution

% make github for this

% direct correlation of eeg activity and BOLD? (e.g. what anthony wants).

% what's going on with some subjects, e.g #2, #3?

results = [];

%% Params
par = EF_Params(subj_id);

%get behavioral info
[res idx] = EF_BehAnalyzer(par);

S.subj_id = subj_id;
S.expt_dir = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/';
S.univar_dir = [S.expt_dir S.subj_id '/' 'analysis_buttonPresses'];
S.eeg_dir = [S.expt_dir S.subj_id '/erpData'];
S.workspace_dir = [S.expt_dir S.subj_id '/' 'mvpa'];
S.patRegSet = {'EF'};
S.patReg = 'trialBetas';
S.patRegDir = [S.expt_dir S.subj_id '/' S.patReg];
S.exp_name = 'AG';

S.createAnalysisMap = false;

%% Condition Parameters
S.conds = 'finger1';
S.dur = sum(par.numvols) * par.TR;

%% Smoothing Parameters
S.use_unsmoothed = 0;
S.smoothTxt = {'smoothed' 'unsmoothed'};

if S.use_unsmoothed
    par.filesForPatterns = par.wascanfiles;
else
    par.filesForPatterns = par.swascanfiles;
end

%% Volume Params
S.roi_name = 'mask.hdr';
S.roi_file = ['/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/' subj_id '/analysis_buttonPresses_bySpectralPower/mask.hdr'];
%S.roi_file = [S.univar_dir '/' S.roi_name];
%S.roi_file = ['/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/' subj_id '/analysis_buttonPresses_byAmp_RTLocked/leftMotorCortex.nii'];
S.vol_info = spm_vol(['/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/' subj_id '/analysis_buttonPresses_bySpectralPower/beta_0001.hdr']);
%S.secondaryMask = [];
S.secondaryMask = ['/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/' subj_id '/analysis_buttonPresses_byAmp_RTLocked/leftMotorCortex.nii'];

%% locking and interpolation
S.interp = true;
S.lock = 'RT';
S.interpText = {'interp' 'noInterp'};

%% TRs
S.TR_weights = {[1] [0 1] [0 0 1] [0 0 0 1] [0 0 0 0 1] [0 0 0 0 0 1]};
S.TRs_to_average_over = 1:length(S.TR_weights);
S.TRByTR = 1;
S.TRPatName = 'patsByTR';

%% Workspace Parameters
S.use_premade_workspace = 1;

if ~iscell(S.TR_weights)
    S.TRDescription = num2str(find(S.TR_weights));
else
    S.TRDescription = 'TRs1Through6';
end

S.condensePatterns = false;
S.useBetaMapPatterns = true;

if S.useBetaMapPatterns
    S.betaMapTxt = 'betaMaps';
else
    S.betaMapTxt = '';
end

S.workspace = fullfile(S.workspace_dir, [S.subj_id '_' S.roi_name '_' S.smoothTxt{S.use_unsmoothed + 1} '_' S.TRDescription '_' S.lock '_' S.interpText{2-S.interp} S.betaMapTxt '.mat']);

%% EEG Params
S.eegToFMRIClassification = true;
S.eegDataType = 'trialdata_RTlock.mat';
S.decimationFactor = 5;
S.ChOI = {'frontal.FM' 'frontal.LFI' 'frontal.LFS' 'central.CM' 'central.LCI' 'central.LCS'};
S.TOI = 400:600; % in units of samples

%% Classifier Params

%S.lambda = [.01 .1 1 10 100 1000 10000];
S.lambda = 10;
S.trainOptsSVM = '-s 3 -t 0 -c ';
S.trainOptsLibLinear = '-q -s 0 -c ';
S.nXvals = 10;
S.classifier = 'svr';
S.nFeats = 1000;
S.FS = false;
S.normalizeFeatures = true;
S.ValidationLambda = 10;

%% classFile
S.classMatDir = ['/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/EEG_to_BOLD/classMats/'];
S.classMatType = 'EF';

%% beta map creation
S.modelSubSample = 1; %0 for no sub-sampling, 1 to subsample only the session of interest
S.betaMapPrefix = 'whole_brain_trial_';
S.bf = 'FIR'; %FIR or HRF
S.num_FIR_bins = 8;

%% images
if S.useBetaMapPatterns
    if S.TRByTR
        idx.TRsExist = false(length(idx.allTrials), S.num_FIR_bins);
        for TRN = 1:S.num_FIR_bins
            S.filenames{TRN} = [];
            for tN = 1:length(idx.allTrials)
                dImg = dir(fullfile(S.patRegDir,[ S.betaMapPrefix prepend(num2str(tN),4) '_TR' num2str(TRN) '.img']));
                if ~isempty(dImg)
                    idx.TRsExist(tN, TRN) = true;
                    S.filenames{TRN} = [S.filenames{TRN}; [S.patRegDir '/' dImg.name]];
                end
            end
            S.preprocPatCondensedName{TRN} = ['betaMaps_' num2str(TRN)];
        end
    else
        dImg = dir(fullfile(S.patRegDir,[ S.betaMapPrefix '*.img']));
        S.filenames{1} = [repmat([S.patRegDir '/'],length(dImg),1) vertcat(dImg.name)];
        S.preprocPatCondensedName = {'betaMaps'};
    end
    S.runs_vector = histc(idx.sess, unique(idx.sess));
    S.preprocPatName = 'spiral_hp_z';
    for f=1:length(S.filenames)
        S.preprocPatCondensedMaskedName{f} = [S.preprocPatCondensedName{f} '_masked'];
        S.img_files{f} = mat2cell(S.filenames{f}, [ones(1,size(S.filenames{f},1))], [size(S.filenames{f},2)]);
    end
    S.patFinalGroupName = [S.preprocPatCondensedName{1}(1:(end-1)) 'group'];
else
    S.preprocPatName = 'spiral_hp_z';
    S.preprocPatCondensedName = 'spiral_hp_z_condensed';
    S.filenames = vertcat(par.filesForPatterns);
    S.img_files = mat2cell(S.filenames, [ones(1,size(S.filenames,1))], [size(S.filenames,2)]);
    S.runs_vector =  par.usedVols;
end



%% Runs Parameters
S.meta_runs = S.runs_vector;
S.num_runs = length(S.runs_vector);
S.num_vols = sum(S.runs_vector);
S.TR = 2;

%% Begin Analysis Here
existWorkspace = exist(S.workspace);

if strcmp(S.lock, 'RT')
    S.onsets = {idx.allTrials + idx.RT};
else
    S.onsets = {idx.allTrials};
end

S.num_conds = size(S.onsets,2);

if S.useBetaMapPatterns
    if (S.use_premade_workspace&&existWorkspace) % if we are supposed to use premade workspace, and one with the correct name exists
        load(S.workspace, 'subj')
    else
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
        
    end
    
else
    if (S.use_premade_workspace&&existWorkspace) % if we are supposed to use premade workspace, and one with the correct name exists
        load(S.workspace, 'subj', 'condensed_regs_of_interest', 'condensed_runs');
    else
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
        
        %% create beta maps through deconvolultion
        if S.createAnalysisMap
            
            xSPMOrig = load(fullfile(S.univar_dir, 'SPM.mat'));
            xSPMOrig = xSPMOrig.SPM;
            
            %create new U struct for trial-specific regressor
            xSPMOrig.Sess.U(2:(end+1)) = xSPMOrig.Sess.U(1:(end));
            xSPMOrig.Sess.U(1).P.name = 'none';
            xSPMOrig.Sess.U(1).dur = 0;
            
            k = xSPMOrig.nscan;
            fMRI_T     = xSPMOrig.xBF.T;
            fMRI_T0    = xSPMOrig.xBF.T0;
            bf = xSPMOrig.xBF.bf;
            V = xSPMOrig.xBF.Volterra;
            
            pat = subj.patterns{end}.mat; %MAKE THIS MORE GENERAL
            sess = xSPMOrig.Sess.C.C;
            constant_col = ones(size(sess,1),1);
            
            idx.sessTR = subj.selectors{1}.mat;
                                
            for a = 1:length(idx.allTrials);
                
                xSPM = xSPMOrig;
                
                thisOnset = idx.allTrials(a);
                
                thisOnsetName = [S.betaMapPrefix prepend(num2str(a),4)];
                
                %remove thisOnset from its original location
                for u=2:length(xSPM.Sess.U);
                    idxToDelete = ismember(xSPM.Sess.U(u).ons, thisOnset);
                    xSPM.Sess.U(u).ons(idxToDelete) = [];
                    xSPM.Sess.U(u).dur(idxToDelete) = [];
                end
                
                    
                    xSPM.Sess.U(1).ons = thisOnset;
                    xSPM.Sess.U(1).dur = 0;
                    xSPM.Sess.U(1).name = {thisOnsetName};
                    
                    %remove conditions with no onsets.
                    idxToRemove = [];
                    for u=2:length(xSPM.Sess.U);
                        if isempty(xSPM.Sess.U(u).ons)
                            idxToRemove = [u idxToRemove];
                        end
                    end
                    xSPM.Sess.U(idxToRemove)=[];
                    
                    % Get inputs, neuronal causes or stimulus functions U
                    U = spm_get_ons(xSPM,1);
                    
                    % Convolve stimulus functions with basis functions
                    X = spm_Volterra(U,bf,V);
                    
                    % Resample regressors at acquisition times (32 bin offset)
                    X = X([0:(k - 1)]*fMRI_T + fMRI_T0 + 32,:);
                    
                    betaMaps = {'1'};
                    
                    %replace the first regressor with a set of FIR regs
                    if strcmp(S.bf, 'FIR')
                        X(:,1) = [];
                        
                        xFIR = zeros(size(X,1),S.num_FIR_bins);
                        
                        for h=1:S.num_FIR_bins
                            thisTR = h + floor(thisOnset/S.TR);
                            
                            % only add the regressor if its value is less
                            % than the number of volumes in a session.
                            if thisTR<=size(X,1)
                                xFIR(thisTR, h)=1;
                                betaMaps{h} = ['TR' num2str(h)];
                            end
                        end
                        
                        X = [xFIR X];

                    end
                
                % only put data and model from the session to which the
                % given trial belongs
                if S.modelSubSample==1

                    idx.thisSess = (idx.sessTR==idx.sess(a));
                    
                    XMat = [X(idx.thisSess,:) constant_col(idx.thisSess)];
                    
                    %remove all-zero columns
                    idx.notAllZero = sum(abs(XMat))~=0;
                    XMat = XMat(:,idx.notAllZero);
                    
                    patMat = pat(:,idx.thisSess);
                    
                    betaMaps(~idx.notAllZero(1:length(betaMaps))) = [];
                else
                    XMat = [X sess constant_col ];
                end
                
                %preallocate betas
                b = nan(size(XMat ,2), size(patMat,1));
                
                tic
                for i = 1:size(patMat,1)
                    thisPat = patMat(i,:);
                    b(:,i) = regress(thisPat', XMat);
                end
                toc
                
                vol_info = S.vol_info;
                voxel_inds = subj.masks{end}.mat;
               
                
                for m = 1:length(betaMaps)
                    
                    outMap=  zeros(vol_info.dim);
                    
                    outMap(voxel_inds==1) = b(m,:);
                    
                    vol_info.dir = S.patRegDir;
                    vol_info.fname = [ vol_info.dir '/' thisOnsetName '_' betaMaps{m} '.img'];
                    
                    S.pat{i}.map{m}.fname = vol_info.fname;
                    
                    if isempty(dir([vol_info.dir]))
                        mkdir(vol_info.dir);
                    end
                    
                    spm_write_vol(vol_info,outMap);
                    
                    fprintf('\n wrote out %s \n', vol_info.fname);
                    
                    d{i} = vol_info.fname;
                end
                
            end
        end
        
        
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
if S.eegToFMRIClassification
    
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
    if S.interp
        interpolatedPats = nan(szPats(1)-1, szPats(2), szPats(3));
        for t=1:size(allPats,3)
            interp_offset = 1+(idx.allTrials(t)-floor(idx.allTrials(t)))/S.TR;
            interp_vals = interp_offset:1:((S.num_FIR_bins-2) + interp_offset);
            interpolatedPats(:,:,t) = interp1(squeeze(allPats(:,:,t)),interp_vals);
        end
        allPats = interpolatedPats;
    end
    
    idx.trialHasEstimateTRWise = ~isnan(squeeze(mean(allPats,2)))';
    idx.trialHasEstimate = all(idx.trialHasEstimateTRWise');
    
    %% load eeg data
    eeg = load(fullfile(S.eeg_dir, S.eegDataType));
    
    if strcmp(S.eegDataType, 'spectraldata3_RT.mat')
        eegDat = squeeze(eeg.data.normpower);
    else
        eegDat = eeg.data.trialdata;
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
    
    eegDat = eegDat(chansToInclude,:,S.TOI);
    
    
    %% reshape and downsample
    sz = size(eegDat);
    
    % reshape the data such that features now combine across time and channels
    eegDatSf = shiftdim(eegDat,2);
    eegDatResh = reshape(eegDatSf, sz(1)*sz(3), sz(2));
    
    % downsample the eeg data
    eegDatReSam = zeros(ceil(( sz(1)*sz(3))/S.decimationFactor), sz(2));
    for k = 1:size(eegDatResh,2)
        eegDatReSam(:,k) = decimate(eegDatResh(:,k),S.decimationFactor)';
    end
    
    eegDat_pats = eegDatReSam';
    
    
    
    %% validation test of BOLD patterns
    for m=1:size(allPats,1)
        
        BOLDButtonPress = squeeze(mean(allPats(m,:,(~idx.noResp & idx.TRsExist(:,m)'))));
        BOLDFix = squeeze(mean(allPats(m,:,(idx.noResp & idx.TRsExist(:,m)'))));
        
        [~, ~, ~, results.tTestBOLDButtonPressVsFix.TR(m).stats] = ttest2(BOLDButtonPress', BOLDFix);
    end
    
     %% select good trials
    S.goodtr = true(size(idx.noResp));
    S.goodtr(S.badtr{par.subNo}) = false;  
    
    thisModIdxEEG = (~idx.noResp & S.goodtr);
    
    %% validation test of EEG
    EEGButtonPress = squeeze(mean(mean(eegDat(:, ~idx.noResp, :),3)));
    EEGFix = squeeze(mean(mean(eegDat(:, idx.noResp, :),3)));
    
    [~, ~, ~, results.tTestEEGButtonPressVsFix.stats] = ttest2(EEGButtonPress', EEGFix);
        
    YRTs = eeg.data.RTs(thisModIdxEEG);
    
    XRTs = eegDat_pats(thisModIdxEEG, :);
    
    %can you predict the RT?
    results.classifyRTWithEEG = x_validation(XRTs,YRTs,S,S.ValidationLambda,'continuous', 'svr');
    
    %can you predict a button press?
    sparseX = sparse(eegDat_pats(S.goodtr,:));
    respsNegPos = 2*idx.noResp(S.goodtr)'-1;
    
    results.classifyButtonPressWithEEG = x_validation(sparseX,respsNegPos,S,S.ValidationLambda,'discrete', 'liblinear');

    if S.TRByTR
        YMat = squeeze(mean(allPats(:,:,:),2))';
    else
        YMat = squeeze(mean(allPats(:,:,:),2));
    end
    
    %% classify

    
    for m=1:size(YMat,2)
        thisModIdx = (~idx.noResp & S.goodtr & idx.TRsExist(:,m)' & idx.trialHasEstimate);
        Y = double(YMat(thisModIdx,m));
        X = eegDat_pats(thisModIdx, :);
        if S.normalizeFeatures
            X = zscore(X);
        end
        results.classifyBOLDWithEEG.TR(m) = x_validation(X,Y,S,S.lambda,'continuous', 'svr');
    end
    
end

end


function [out] = x_validation(X,Y,S,lambda,predictionType,classifier)
    % cross validation index
    cv_tr = shuffle(ceil(S.nXvals*(1:size(X,1))/size(X,1)));
    
    yPred = nan(size(Y));

    %penalty range
    for l = 1:length(lambda)
        thisL = lambda(l);
        thisLStr = num2str(thisL);        
        
        for i=1:10
            idxTR = (cv_tr~=i);
            idxTE = (cv_tr==i);
            
            XTrain = X(idxTR, :);
            YTrain = Y(idxTR);
            
            XTest = X(idxTE, :);
            YTest = Y(idxTE);
            
            %feature selection based on mutual information
            mi = nan(size(XTrain,2),1);
            if S.FS
                for t = 1:size(XTrain,2)
                    mi(t) = mutualinfo(YTrain, XTrain(:,t));
                end
                
                mi_sorted = sort(mi, 'descend');
                thresh = mi_sorted(S.nFeats);
                XTrain = XTrain(:,mi>=thresh);
                XTest = XTest(:,mi>=thresh);
            end
            
            if strcmp(classifier, 'svr')
                % svr
                trainOpts = [S.trainOptsSVM thisLStr];
                model = svmtrain(YTrain, XTrain, trainOpts);
                [yPred(idxTE)] = svmpredict(YTest, XTest, model);
            elseif strcmp(classifier, 'ridge')
                %ridge regression
                trainOpts = [S.trainOptsSVM thisLStr];
                b = ridge(YTrain,XTrain,thisL, 0);
                yPred(idxTE) = [ones(size(XTest,1),1) XTest]*b;
            elseif strcmp(classifier, 'liblinear')
                trainOpts = [S.trainOptsLibLinear thisLStr];
                model = train(YTrain, XTrain, trainOpts);
                yPred(idxTE) = predict(YTest, XTest, model);
            end
        end
        
        out.mod.Y = Y;
        out.mod.YPred = yPred;
        
        if strcmp(predictionType, 'continuous')
            [out.rObsVsPred(l), out.pObsVsPred(l)] = corr(Y,yPred);
        elseif strcmp(predictionType, 'discrete')
            accuracy= Y==yPred;
            [out.meanAccuracyOverClasses(l)] = mean([mean(accuracy(Y==1)) mean(accuracy(Y==-1))]); 
        end
    end
end

