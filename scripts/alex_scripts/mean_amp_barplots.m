% bar plot comparing amplitudes


% run this from my directory, I think it might need some functions in
% it.

% load data...
%load('/biac4/wagner/biac3/wagner5/alan/eegfmri/fmri_data/erp_data/STGA_STPA_groupdata_evonset_ztrialdata_bc1.mat')

ROI = 'PM';
chans = data.rois.(ROI);
smfactor = 5;

allsubj  = [4,6,8:13 15:21];
subj = [4,6,8:10, 11:12, 15:18, 20:21];
[~,subidx]=intersect(allsubj,subj);
subidxl = false(numel(allsubj),1);
subidxl(subidx)=1;

window = [0.4 0.8];
window_str = strrep(sprintf('%0.2gto%0.2gsec',window(1),window(2)),'.','p');
trial_time = data.trialtime;
fs = data.fs_comp;
[~,samps]=histc([window(1):1/fs:window(2)],trial_time);

data.trial_type.RemHC_hits =data.trial_type.Rem_hits|data.trial_type.HC_hits;
conds = {'RemHC_hits', 'LC_hits','LC_cr','HC_cr'};
condsshort = {'RHCH','LCH','LCCR','HCCR'};

temp = squeeze(mean(data.ztrialdata(:,chans,:,:),2));
nsubs = numel(subidx);
ntrials = 400;
nconds = numel(conds);

S=[];
X=[];
Y=[subj]';
condstr = [];
for c = 1:nconds
    S.([conds{c} '_trials']) = (logical(data.trial_type.(conds{c}))&data.goodtrials)';
    
    for s= 1:nsubs
        temp2 = smooth(squeeze(mean(temp(subidx(s),S.([conds{c} '_trials'])(subidx(s),:),:),2)),smfactor,'rloess');
        S.(conds{c})(s) = mean(temp2(samps));
    end
    
    X = [X; S.(conds{c})' repmat(c,nsubs,1) (1:nsubs)'];
    bar_vals(c)=mean(S.(conds{c}));
    condstr = [condstr '_' condsshort{c}];
    
    Y = [Y S.(conds{c})'];
end

T=ar_rmanova1(X);
bar_ers = sqrt(T.MSE)/sqrt(nsubs);


h=figure(1);clf;
colors={'r','g','b','c'};
%colors={'r','b','c'};

for c=1:nconds
    bar(c,bar_vals(c),1,colors{c}); hold on
    errorbar(c,bar_vals(c),bar_ers,'k','linewidth',2)
end

set(gca,'xtick',1:nconds);
set(gca,'xticklabel',condsshort);

ylabel('Mean Zscore','fontsize',14)
xlim([0.4 nconds+0.6])

mainpath = '/biac4/wagner/biac3/wagner5/alan/eegfmri/fmri_data/erp_data/';
str = ['erp_figures/mean_bar_plots/allsubj_' ROI 'Mean_Zscore' '_' window_str condstr];
print(gcf,'-loose','-dtiff',[mainpath str])

str = [ROI 'MeanZscore_' window_str condstr];
save([mainpath str],'Y','condsshort')
