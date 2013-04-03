% datas together the subjects data


clearvars
data=[];
data.subjs = [4 6 8 9 10 11 12 13];
data.nsubj = numel(data.subjs);

data.conds = {'old','new','hits','cr','Rem_hits','HC_hits','LC_hits','HC_cr',...
    'LC_cr','FA','misses'};


type = 'ztrial'; % zscore data or regular data
lock = 'RT';

bc=0; % baseline correction
data.(type) = [];
data.goodtrials =[];
data.RTs=[];
channel_rois
data.chnls = chnls;
tpr =80; % trials per run
tps = 400; % trials per subject
ns = 8; % number of subjects


mainpath1 = '/Volumes/EXT_HD/ON_eegfmri/';

foldernames = {'ef_091211','ef_091511','ef_092111','ef_092211','ef_092711',...
    'ef_092911','ef_100511','ef_101411'};

data.trial_type = cell2struct(cell(numel(data.conds),1),data.conds,1);
eeg_fmri_on_subject_info

for s = 1:numel(data.subjs)
    
    %datapath =['s' num2str(s) '/processed_data/results/data'];
    temp = load([mainpath1 's' num2str(data.subjs(s)) ...
        '/clean_data2/results/trial_data_lp30hz_zscored_' lock]);
    
    for c = data.conds
        data.trial_type.(char(c)) = vertcat(data.trial_type.(char(c)),...
            temp.data.([char(c) '_trials']));
    end
    
    data.samples = temp.data.samples;
    data.goodtrials = vertcat(data.goodtrials,temp.data.goodtrials);
    data.dP{s} = temp.data.dP;
    data.RTs = vertcat(data.RTs,temp.data.RTs');
    
    data.badtrials{s} = S(data.subjs(s)).badtr;
    data.badch{s} = S(data.subjs(s)).badch;
    
    if bc
        mean_all = mean(temp.data.([type 'data'])(:,:,1:-data.samples(1)),3);
        bsl_mean = repmat(mean_all,[1,1,temp.data.nsamples]);
        trial_data = temp.data.([type 'data']) - bsl_mean;
        temp=[];
    else
        trial_data = temp.data.([type 'data']);
        temp=[];
    end
    
    data.(type) = cat(2,data.(type),trial_data);
    
end

if ~exist([mainpath1 '/group_data/'],'dir')
    mkdir([mainpath1 '/group_data/'])
end

save(sprintf('%s/group_data/STGA_STPA_groupdata_%s_%s_bc%d',mainpath1,...
    lock,type, bc),'data','-v7.3')