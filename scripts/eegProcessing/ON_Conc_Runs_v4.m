% after running ON_EEGFMRI_PRA_PreProcess
% run this to concatenate the trials of all the runs



clearvars
%mainpath = '/Volumes/EXT_HD/ON_eegfmri/';

subjs = [4 6 8:13];
lock = 'RT';
%foldernames = {'ef_040412'  'ef_040512' 'ef_040712' 'ef_041112' 'ef_042912' 'ef_050112'};
foldernames = {'ef_091211' 'ef_091511' 'ef_092111' 'ef_092211' 'ef_092711' 'ef_092911' 'ef_100511' 'ef_101411'};

savepath = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/';

% isOpen = matlabpool('size') > 0;
% 
% if isOpen
%     matlabpool close
% end

%matlabpool open 2

%parfor s =  1:numel(subjs);

for s =  1:numel(subjs);    
    sub = subjs(s);
    
    mainpath1 = fullfile('/biac4/wagner/biac3/wagner5/alan/eegfmri/fmri_data/', foldernames{s}, '/erpData/');
   
   % mainpath1 = ['/Volumes/EXT_HD/ON_eegfmri/s' num2str(sub) '/clean_data2/'];
   % mainpath2 = ['/Volumes/EXT_HD2/ON_eegfmri/s' num2str(sub) '/clean_data2/'];
    
    data =[];
    data.subjs = sub;
    data.runs = 1:5;
    
    if strcmp(lock,'RT')
        data.dur = [-1 0.5];
    elseif strcmp(lock,'evonset')
        data.dur = [-0.5 1.5];
    end
    
    data.fs = 500;
    data.trtime = data.dur(1):1/data.fs:data.dur(2);
    data.dur_samp = data.dur*data.fs;
    data.Nchan = 256;
    data.trialsperrun = 80;
    %data.LPfilter = 40;
    
    data.oldtrials =[];
    data.newtrials =[];
    
    data.hitstrials =[];
    data.crtrials =[];
    
    data.REMtrials =[];
    data.HChitstrials =[];
    data.LChitstrials =[];
    
    data.HCcrtrials =[];
    data.LCcrtrials =[];
    
    data.FAtrials =[];
    data.MISStrials =[];
    
    data.RTs =[];
    data.trialdata =[];
    
    for r = data.runs
        datapath =['r' num2str(r)];
        temp1 = load([mainpath1 datapath '/data']);
        temp2 = load([mainpath1 datapath '/behdata.mat']);
        S = temp1.S;
        S.behavdata = temp2.data; 
        temp1 = [];
        temp2 = [];
        
        %S.recons_signal = eegfmri_filt_v3(S.recons_signal,data.fs,data.LPfilter);
        % duration of events
        dur = numel(data.dur_samp(1):data.dur_samp(2));
        
        % trial markers
        if strcmp(lock,'RT')
            marks = S.ev_markers' + ceil(S.behavdata.RT*S.fs);
        elseif strcmp(lock,'evonset')
            marks = S.ev_markers';
        end
        
        % trial markers
        S.all_trials_markers = [marks+data.dur_samp(1), ...
            marks, marks+data.dur_samp(2)];
        
        %zscore each run independently.  
        recons_signal = (S.recons_signal')';
        
        S.all_trials_clean = zeros(S.Nchan,size(S.all_trials_markers,1),dur);
        % trial matrix
        for tr = 1: size(S.all_trials_markers,1)
            S.all_trials_clean(:,tr,:) = recons_signal(:,S.all_trials_markers(tr,1):S.all_trials_markers(tr,3));
        end
        
        % baseline correct
%         mean_all = mean(S.all_trials_clean(:,:,1:-data.dur_samp(1)),3);
%         bsl_mean = repmat(mean_all,[1,1,dur]);
%         S.all_trials_clean = S.all_trials_clean - bsl_mean;
        
        % if for some reason data is not available, concantenate behavioral data
        %for trials that have eeg data
        end_marker = numel(S.ev_markers);
        
        % old and new trials
        data.oldtrials = vertcat(data.oldtrials,S.behavdata.old(1:end_marker));
        data.newtrials = vertcat(data.newtrials,S.behavdata.new(1:end_marker));
        
        % hits and cr trials
        data.hitstrials = vertcat(data.hitstrials,S.behavdata.hits(1:end_marker));
        data.crtrials = vertcat(data.crtrials,S.behavdata.cr(1:end_marker));
        
        % Rem hits, HC hits, LC hits trials
        data.REMtrials = vertcat(data.REMtrials,S.behavdata.Rem_hits(1:end_marker));
        data.HChitstrials = vertcat(data.HChitstrials,S.behavdata.HC_hits(1:end_marker));
        data.LChitstrials = vertcat(data.LChitstrials,S.behavdata.LC_hits(1:end_marker));
        
        % HC cr and LC cr
        data.HCcrtrials = vertcat(data.HCcrtrials,S.behavdata.HC_cr(1:end_marker));
        data.LCcrtrials = vertcat(data.LCcrtrials,S.behavdata.LC_cr(1:end_marker));
        
        %misses and FA
        data.FAtrials = vertcat(data.FAtrials,S.behavdata.FA(1:end_marker));
        data.MISStrials = vertcat(data.MISStrials,S.behavdata.misses(1:end_marker));
        
        % RTs
        data.RTs = vertcat(data.RTs,S.behavdata.RT(1:end_marker));
        
        % trials
        data.trialdata = cat(2,data.trialdata,S.all_trials_clean);
        
        %S = [];
    end

    
    if ~exist([mainpath1 '/results'],'dir')
        mkdir([mainpath1 '/results'])
        %mkdir([mainpath2 '/results'])
    end
    %parsave([mainpath1 '/results/trialdataLP' num2str(data.LPfilter)],data,'data','-v7.3')
    parsave([mainpath1 '/results/trialdataLP125Hz'],data,'data','-v7.3')
    %parsave([savepath foldernames{s} '/erpData/trialdata_' lock 'lock_LP125Hz'],data,'data','-v7.3')
end
%matlabpool close