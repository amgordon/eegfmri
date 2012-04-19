function [S par] = EF_EEG_ResponseParams(subj_id)

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