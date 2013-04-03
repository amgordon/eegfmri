% for each subject calculate the bad trials, downsample and store all the
% data in a single file.

clearvars

% add eeglab path
addpath(genpath('/Users/alexgonzalez/Documents/MATLAB/Mylib/eeglab'))

% declare main working directory path
mainpath = '/biac4/wagner/biac3/wagner5/alan/eegfmri/';
datapath = [mainpath 'fmri_data/'];

% declare data struct
data=[];
%subjects folder names
data.subjs = {'ef_091211','ef_091511','ef_092111','ef_092211','ef_092711',...
    'ef_092911','ef_100511','ef_101411', 'ef_040412',  'ef_040512', ...
    'ef_040712', 'ef_040712_2' , 'ef_041112', 'ef_042912'};
% number of subejects
data.nsubj = numel(data.subjs);
% all possible conditions of interest
data.conds = {'old','new','hits','cr','Rem_hits','HC_hits','LC_hits','HC_cr',...
    'LC_cr','FA','misses'};
data.trial_type = cell2struct(cell(numel(data.conds),1),data.conds,1);

% internal variables
bc=1; % baseline correction
type = 'ztrialdata'; % {'ztrialdata', 'trialdata'}
lock = 'evonset'; % {'RT', 'evonset'}
ext = 'LP30Hz_zscored.mat';
tpr =80; % trials per run
tps = 400; % trials per subject
%nruns = 5; % number of runs per subject
nchan = 256;
cmp=5;
ns = data.nsubj; % number of subjects

% other pre-allocations
data.(type) = nan(ns,nchan,tps,200);
data.goodtrials =[];
data.RTs=[];
data.baseline_corrected = bc;
% trial rejection threshold
data.trial_rejection_threh = 2;
% use these channels to classify a trial as bad
data.trial_rejection_channels = 'LPS';

% channels sorted by region structure
data.rois = channel_rois;

% data compression flag
data.compress = 1; 
% downsamples the data to a factor, compression flag needs to be 1 for 
% this. Note that this is just downsampling and that it is assuming that
% appropiate filtering has been done on the continous time signal to avoid
% aliasing at the compression rate.
data.compress_factor=cmp; 

% for each subject.
for s = 1:ns
    
    temp = load([datapath data.subjs{s} '/erpData/trial_data_' ext]);
    
    for c = data.conds
        data.trial_type.(char(c))(:,s) = temp.data.([char(c) '_trials'])==1;
    end
    
    
    % zscore and find bad trials based on amplitude data
    channels = data.rois.(data.trial_rejection_channels);
    thr = data.trial_rejection_threh;
    X = temp.data.(type)(channels,:,:);
    [nch ntr nsamps]=size(X);
    
    % zscore the channels for trial rejection purposes
    if strcmp(type,'trialdata')
        XZ = zscore(X(:,:),[],2);
        XZ = reshape(XZ,[nch ntr nbins]);
        % collapse across channels
        X = permute(XZ,[2 1 3]); XZ=[];
    else
        X = permute(X,[2 1 3]);
    end
    
    
    % calculate the trial covariane matrix
    trial_var = var(X(:,:),0,2);
    high_var_bound = nan(size(trial_var));
    % calculate an upper bound on the variance for each trial per run
    for r = 1:nruns
        % index of the trials for a run in the matrix
        trial_idx = ((r-1)*tpr+1):(tpr*r);
        % calculate the bound by trial
        high_var_bound(trial_idx) = mean(trial_var(trial_idx))...
            + thr*std(trial_var(trial_idx));
    end
    
    % deterime good trials
    data.goodtrials(:,s) = (trial_var < high_var_bound);
    
    % store subjects dP and RT
    data.dur = temp.data.dur;
    data.fs = temp.data.fs;
    data.dP(:,s) = temp.data.dP;
    data.RTs(:,s) = temp.data.RTs;
    data.badch{s} = temp.data.badch;
    
    if data.compress
        % toss out samples at the end.
        nsamps_comp = floor((tps*nsamps/cmp)/tps);
        X = temp.data.(type)(:,:,1:nsamps_comp*cmp);
        nch = size(X,1);
        % downsample;
        X = downsample(X(:,:)',cmp)';
        X = reshape(X,[nch,ntr,nsamps_comp]);
        data.fs_comp = data.fs/cmp;
        data.trialtime = (data.dur(1):1/data.fs_comp:(data.dur(2)-1/data.fs_comp));
        data.trialsamps = data.trialtime*data.fs_comp;
    else
        X = temp.data.(type);
        data.trialtime = (data.dur(1):1/data.fs:data.dur(2));
        data.trialsamps = temp.data.samples;
    end
    
    if bc
        mean_all = mean(X(:,:,1:-data.trialsamps(1)),3);
        X = bsxfun(@minus,X,mean_all);
    end
    temp=[];
    
    data.(type)(s,:,:,:) = X;
    
end

save(sprintf('%s/erp_data/STGA_STPA_groupdata_%s_%s_bc%d',datapath,...
    lock,type, bc),'data','-v7.3')