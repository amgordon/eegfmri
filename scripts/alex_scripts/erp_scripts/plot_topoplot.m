clear all
type = 'ztrial';
bc='0';
lock = 'RT';
mainpath = '/Volumes/EXT_HD/ON_eegfmri/';
load([mainpath '/group_data/STGA_STPA_groupdata_' lock '_' type '_bc' bc]);

eeglab_path = '/Volumes/EXT_HD/MATLAB/Mylib/eeglab';
addpath(genpath(eeglab_path))

locs_file = '/Volumes/EXT_HD/MATLAB/MemoryLabResearch/EEGFMRI_cleaning_scripts/eeg_fmri_fncs/GSN256.sfp';
chanlocs = readlocs(locs_file);

%%
allconds = {'old','new','hits','cr','Rem_hits','HC_hits','LC_hits','HC_cr',...
    'LC_cr','FA','misses'};

mainpath = '/Volumes/EXT_HD/ON_eegfmri/';
type = 'ztrial';
bc='0';
close all
conds = {'Rem_hits', 'HC_cr'};
fs = 500;


window = [0.6 1];
window_str = sprintf('0_%0.2gto%0.2gsec',window(1)*10,window(2));
if strcmp(lock,'RT')
    trial_time = -1:1/fs:1;
else
    trial_time = -0.5:1/fs:1.5;
end
samp_start = (window(1)+trial_time(1))*fs;
samp_end = samp_start+(window(2)+trial_time(1))*fs;
samps = round(samp_start:samp_end);

nsamps = numel(trial_time);
area = 'parietal';
rois = {'LPS'};%fieldnames(group.chnls.parietal)';
ch_analysis = 0;
smfactor = 20;
colors = {'r','b','g','k','m','c','y'};
allsubj  = [4,6,8:13];
subj = [4,6,8:9,11:12];
nsubj = 8;
ntrialspersubj = 400;

strial_all = false(nsubj*ntrialspersubj,1);


mean(mean(data.ztrial(:,[50 60],samps),3),2)

x_mean = zeros(numel(subj)+1,numel(conds),256,'single');
%x = zeros(numel(conds),256,'single');
s_cnt = 1;

for s = [subj 100]
    
    if s <= max(subj)
        strials = false(1,nsubj);
        strials(s==allsubj) = true;
        strials = reshape(repmat(strials,ntrialspersubj,1),1,[])';
        subj_str = sprintf('subj%i',s);
        strial_all = strial_all | strials;
        
    else
        strials = strial_all;
        subj_str = 'allsubj';
        
    end
    
    
    %title_str = sprintf('%s ROI Cleaned ERP: ',roi);
    cond_str =[];
    
    for c = 1:numel(conds)
        cond = logical(data.trial_type.(conds{c}))&strials;
        
        if s <= max(subj)
            x_mean(s_cnt,c,:) = mean(mean(data.(type)(:,cond,samps),3),2);
        else
            x_mean(s_cnt,c,:) = squeeze(mean(x_mean(1:(end-1),c,:),1));
        end
        cond_str = [cond_str ' ' conds{c}];
        
    end
    y = squeeze(x_mean(s_cnt,1,:)) - squeeze(x_mean(s_cnt,2,:));
    figure;topoplot([0 0 0 y' 0],chanlocs)
    title([ subj_str ' topoplot ' cond_str ' ' window_str])
    colorbar
    
    path = [mainpath 'group_data/topo_plots/' subj_str '/' lock '/'];
    mkdir(path)
    str = sprintf('%s_%s_%s_%s_bc%s',subj_str,cond_str,type,window_str,bc);
    print(gcf,'-loose','-dtiff',[path str])
    
    
    s_cnt= s_cnt+1;
    
end


%%

% [m,n] = size(W);
% figure;
% for i = 1: m
%     subplot(6,5,i)
%     figure;topoplot([0 0 0 W(i,:) 0],chanlocs)
% end