
% this script is to average events by condition within a run and to put all
% runs together


clearvars

subjs = [4 6 8:13];
runs = 1:5;
nruns = 5;
lock = 'RT';
comp=4;

foldernames = {'ef_091211','ef_091511','ef_092111','ef_092211','ef_092711',...
    'ef_092911','ef_100511','ef_101411'};


conds = {'old','hits','HC_hits','Rem_hits','misses','new','cr','HC_cr', ...
    'LC_cr','FA'};
load bandpassfilters3.mat
trials_ind_start = 1:80:400;

eegfmri_ON_trialch_info

% isOpen = matlabpool('size') > 0;
% 
% if isOpen
%     matlabpool close
% end
% 
% matlabpool open 2


for s =  1:numel(subjs);
    
    savepath = fullfile('/biac4/wagner/biac3/wagner5/alan/eegfmri/fmri_data/', foldernames{s}, '/');

    sub = subjs(s);
    %mainpath = ['/Volumes/EXT_HD/ON_eegfmri/s' num2str(sub) '/clean_data2/'];
    
    % data and parameters
    data =[];
    data.subjs = sub;
    data.runs = runs;
    data.spv = 35; spv = 35;% slices per volume
    data.fs = 500;
    data.comp = comp;
    data.fs_comp = data.fs/comp;
    data.bandpassfilts = f;
    data.nbands = 1; nbands = 1;
    data.Nchan = 256; nchan = 256;
    data.badch = S.badch{sub};
    data.trialsperrun = 80; tpr = 80; 
    data.goodtrials =  true(nruns*data.trialsperrun,1);
    data.goodtrials(S.badtr{sub}) = false;
    
    if strcmp(lock,'RT')
        data.dur = [-1 1];
    elseif strcmp(lock,'evonset')
        data.dur = [-0.5 1.5];
    end
    
    data.samples =  data.dur(1)*data.fs_comp:data.dur(2)*data.fs_comp;
    data.nsamples = numel(data.samples); spt =numel(data.samples); % samples per trial
    
    % preallocate variables
    ntrials = nruns*tpr;
    data.normpower = zeros(data.Nchan,nbands,ntrials,spt,'single');
    data.phase = zeros(data.Nchan,nbands,ntrials,spt,'single');
    data.mP = zeros(nruns,data.Nchan,nbands,'single');
    data.marks = zeros(size(ntrials));
    
    % preallocate fields for conditition trials
    for c =1:numel(conds)
        data.([conds{c} '_trials']) = false(ntrials,1);
    end
   
   % for every run
    trial_count =1; 
    for r = runs
        
        % load the data for that run
        datapath =['r' num2str(r)];
        temp1 = load([savepath 'erpData/' datapath '/data']);
        temp2 = load([savepath 'erpData/' datapath '/behdata.mat']);
        contdata = temp1.S; temp1 =[];
        behavdata = temp2.data; temp2 =[];
        
        for c =1:numel(conds)
            data.([conds{c} '_trials'])(trial_count:trial_count+tpr-1) = ...
                behavdata.(conds{c});
        end
        % behavioral performance
        data.dP(r) = behavdata.dP;
        data.RTs(trial_count:trial_count+tpr-1) = behavdata.RT;  
        
        % samples to calculate power average power in the run
        % using 2 TRs into the run and 2 TR before the run is over
        signal_ind = contdata.MR_Pulses(2*spv):contdata.MR_Pulses(end-2*spv);
        
        %convert markers after resampling
        signal_ind = floor(signal_ind/comp);
        ev_markers = contdata.ev_markers - signal_ind(1);
        markers = floor(ev_markers'/comp);
        marks = zeros(size(markers));
        
        % trial markers
        if strcmp(lock,'RT')
            marks = contdata.ev_markers' + ceil(behavdata.RT*data.fs);
            marks = floor(marks/comp);
        elseif strcmp(lock,'evonset')
            marks = floor(contdata.ev_markers/comp)';
        end
        
        ev_onset{r} = marks;
        
        % sample matrix
        all_trials_markers = [marks+data.samples(1), ...
             marks, marks+data.samples(end)];
        samps = repmat(marks,1,spt) + repmat(data.samples,tpr,1);
        
        for ch = 1:data.Nchan
            display(sprintf('Subject %i, Run %i, Decomposing Channel %i',sub,r,ch))
            % create bandpass signal and compress
            x = eegfmri_multibandpass_v3(contdata.recons_signal(ch,:),f,comp);
            % calculate analytic signal
            x = hilbert(x')';
            
            % calculate power
            pw = abs(x).^2;            
            data.mP(r,ch,:) = mean(pw(:,signal_ind),2);
            
            % phase calculation
            data.phase(ch,:,trial_count:trial_count+tpr-1,:) = reshape(angle(x(:,samps)),nbands,tpr,spt);
            
            % calc meanPower; one volume into the run to avoid transients
            %data.mP(r,ch,:) = mean(pw(:,signal_ind),2);
          
            % normalize trials by mean power
            data.normpower(ch,:,trial_count:trial_count+tpr-1,:)  = reshape(...
                diag(1./squeeze(data.mP(r,ch,:)))*pw(:,samps) ...
                ,nbands,tpr,spt);
             
            x = []; pw =[];
            
        end
        contdata =[];
        trial_count = trial_count + tpr;
    end
    
    data.ev_onset = vertcat(ev_onset{:});
    
    % save
    %parsave([mainpath '/results/spectral_data3'],data,'data','-v7.3')
    parsave([savepath 'erpData/spectraldata3_' lock],data,'data','-v7.3')
end
matlabpool close