
% %
% clear all
% type = 'ztrialdata';
% bc='1';
% lock = 'evonset';
% mainpath = '/biac4/wagner/biac3/wagner5/alan/eegfmri/fmri_data/';
% load([mainpath '/erp_data/STGA_STPA_groupdata_' lock '_' type '_bc' bc]);

eeglab_path = '/Documents/MATLAB/Mylib/eeglab';
addpath(genpath(eeglab_path))

locs_file = 'GSN256.sfp';
chanlocs = readlocs(locs_file);
%%
data.trial_type.RemHC_hits =data.trial_type.Rem_hits|data.trial_type.HC_hits;
%%
type = 'ztrialdata';
bc='1';
lock = 'evonset';
mainpath = '/biac4/wagner/biac3/wagner5/alan/eegfmri/fmri_data/';
allconds = {'old','new','hits','cr','Rem_hits','HC_hits','LC_hits','HC_cr',...
    'LC_cr','FA','misses'};

close all

conds = {'RemHC_hits', 'HC_cr'};
trial_time = data.trialtime;

area = 'parietal';
rois = {'LFI'};%fieldnames(group.chnls.parietal)';

window = [0.3 0.5];

fs = data.fs_comp;
smfactor = 10;
colors = {'r','b','g','k','m','c','y'};

allsubj  = [4,6,8:13 15:21];
subj = [4,6,8:9,10:12 15:18 20:21];

nsubj = numel(subj);
[~,subidx]=intersect(allsubj,subj);

ntrialspersubj = 400;
strial_all = false(nsubj*ntrialspersubj,1);
window_str = strrep(sprintf('%0.2gto%0.2gsec',window(1),window(2)),'.','p'); 
[~,samps]=histc([window(1):1/fs:window(2)],trial_time);

x_mean = zeros(numel(subj)+1,numel(conds),256,'single');

%x = zeros(numel(conds),256,'single');
s_cnt = 1;

for s = [subidx 100]
    
   if s <= max(subj)
        splots=s;
        strials = false(1,nsubj);
        strials(s==allsubj) = true;
        strials = reshape(repmat(strials,ntrialspersubj,1),1,[])';
        subj_str = sprintf('subj%i',allsubj(s));
        strial_all = strial_all | strials;
    else
        splots=1:numel(s);
        strials = strial_all;
        subj_str = 'allsubj';
        
    end
    
    %title_str = sprintf('%s ROI Cleaned ERP: ',roi);
    cond_str =[];
    
    for c = 1:numel(conds)
         
        
        if s <= max(subj)
            cond = logical(data.trial_type.(conds{c})(:,s))&data.goodtrials(:,s);
            x_mean(s_cnt,c,:) = mean(mean(data.(type)(s,:,cond,samps),4),3);
        else
            x_mean(s_cnt,c,:) = squeeze(mean(x_mean(1:(end-1),c,:),1));
        end
        cond_str = [cond_str ' ' conds{c}];
        
    end
    y = squeeze(x_mean(s_cnt,1,:)) - squeeze(x_mean(s_cnt,2,:));
    figure;topoplot([0 0 0 y' 0],chanlocs,'shading','interp')
    title([ subj_str ' topoplot ' cond_str ' ' window_str])
    colorbar
    
    path = [mainpath '/erp_data/erp_figures/topo_plots/'];
    %mkdir(path)
    
    str = sprintf('%s_%s_%s_%s_%s_bc%s_%s',subj_str,window_str, ...
                conds{1},conds{2},type,bc,lock);
    print(gcf,'-loose','-dtiff',[path str])
    
    
    s_cnt= s_cnt+1;
    
end



x1=squeeze(x_mean(1:(end-1),1,:));
x2=squeeze(x_mean(1:(end-1),2,:));
[~,~,~,t]=ttest(x1,x2);
figure;topoplot([0 0 0 t.tstat 0],chanlocs,'shading','interp')
title([ subj_str ' topoplot tstat ' cond_str ' ' window_str])
colorbar
str = sprintf('%s_%s_%s_%s_%s_bc%s_%s_tstat',subj_str,window_str, ...
                conds{1},conds{2},type,bc,lock);
print(gcf,'-loose','-dtiff',[path str])
    
    

%%

% [m,n] = size(W);
% figure;
% for i = 1: m
%     subplot(6,5,i)
%     figure;topoplot([0 0 0 W(i,:) 0],chanlocs)
% end