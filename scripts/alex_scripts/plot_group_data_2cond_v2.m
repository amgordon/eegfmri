% plot conditions of interest, shaded error bar and two sample ttest to
% determine significance between 2 conditions

% %
% clear all
% type = 'ztrialdata';
% bc='1';
% lock = 'evonset';
% mainpath = '/biac4/wagner/biac3/wagner5/alan/eegfmri/fmri_data/';
% load([mainpath '/erp_data/STGA_STPA_groupdata_' lock '_' type '_bc' bc]);

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


conds = {'Rem_hits', 'HC_cr'};
trial_time = data.trialtime;

area = 'parietal';
rois = {'LPI'};%fieldnames(group.chnls.parietal)';

ch_analysis = 0;
fs = data.fs_comp;
smfactor = 10;
colors = {'r','b','g','k','m','c','y'};

allsubj  = [4,6,8:13 15:21];
subj = [4,6,8:9,10:12 15:18 20:21];
nsubj = numel(subj);
[~,subidx]=intersect(allsubj,subj);

ntrialspersubj = 400;

strial_all = false(nsubj*ntrialspersubj,1);

x_mean=[];
for s = [subidx]
    
    
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
    
    for r = rois
        roi = r{1};
        chnls = data.rois.(roi);
        
        if ch_analysis
            %%
            for ch =chnls
                figure;
                title_str = sprintf('Channel %i Cleaned ERP: ',ch);
                
                x =[];
                for c = 1:numel(conds)
                    cond = logical(data.trial_type.(conds{c})) &strials;
                    
                    x{c} = squeeze(data.(type)(ch,cond,:));
                    x_mean = smooth(mean(x{c},1),smfactor);
                    x_std = std(x{c},[],1);
                    x_sem = (x_std)/sqrt(sum(cond));
                    
                    shadedErrorBar(trial_time(10:end-10),x_mean(10:end-10),x_sem(10:end-10),...
                        {[colors{c} '-'],'LineWidth',2,'markerfacecolor',colors{c}},1);
                    hold on
                    title_str = [title_str sprintf('%s is %s ',colors{c},conds{c})];
                end
                
                title([subj_str title_str])
                xlabel('Time(s)')
                axis tight
                line(xlim,[0 0], 'color','k','linestyle','--')
                lims_y = ylim;
                t0 = find(trial_time==0);
                lims_y(2) = lims_y(2)+lims_y(2)/5;
                line([t0 t0],ylim, 'color','k','linestyle','--')
                
                [H] = ttest2(x{1},x{2});
                sig_points = find(H)/fs + trial_time(1);
                if ~isempty(sig_points)
                    nsegments = sum(diff(sig_points)>2/fs)+1;
                    seg_ind = [1 find(diff(sig_points)>2/fs)+1 numel(sig_points)];
                    
                    for i = 1:nsegments
                        lx_start = sig_points(seg_ind(i));
                        if seg_ind(i) ==numel(sig_points)
                            lx_end = sig_points(seg_ind(i+1));
                        else
                            lx_end = sig_points(seg_ind(i+1)-1);
                        end
                        line([lx_start lx_end], [lims_y(2) lims_y(2)], ...
                            'color', 'k', 'linewidth',2)
                    end
                    ylim([lims_y(1) lims_y(2)+lims_y(2)/5])
                    line([t0 t0],ylim, 'color','k','linestyle','--')
                end
                hold off
            end
        else
            %%
            figure;
            title_str = sprintf('%s ROI Cleaned ERP: ',rois{1});
            
            x =[];
            for c = 1:numel(conds)
                
                
                cond = logical(data.trial_type.(conds{c})(:,s))&data.goodtrials(:,s);
                
                x{c} = squeeze(mean(data.(type)(s,chnls,cond,:),2));
                %x_mean(:,c,s) = smooth(mean(x{c},1),smfactor,'rloess');
                x_mean(:,c,s) = smooth(mean(x{c},1),smfactor);
                x_std = std(x{c},[],1);
                x_sem = (x_std)/sqrt(sum(cond));
                
                shadedErrorBar(trial_time(10:end-10),x_mean(10:end-10,c,s),x_sem(10:end-10),...
                    {[colors{c} '-'],'LineWidth',2,'markerfacecolor',colors{c}},1);
                hold on
                title_str = [title_str ' ' colors{c} ' is ' conds{c}];
                
            end
            title([subj_str title_str])
            xlabel('Time(s)')
            if strcmp(type,'ztrialdata')
                ylabel('zscore')
            else
                ylabel('\mu Volts')
            end
            
            axis tight
            line(xlim,[0 0], 'color','k','linestyle','--')
            lims_y = ylim;
            
            t0 = trial_time(trial_time==0);
            
            lims_y(2) = lims_y(2)+lims_y(2)/5;
            [H] = ttest2(x{1},x{2});
            sig_points = find(H)/fs + trial_time(1);
            
            if ~isempty(sig_points)
                nsegments = sum(diff(sig_points)>2/fs)+1;
                seg_ind = [1 find(diff(sig_points)>2/fs)+1 numel(sig_points)];
                
                for i = 1:nsegments
                    lx_start = sig_points(seg_ind(i));
                    if seg_ind(i) ==numel(sig_points)
                        lx_end = sig_points(seg_ind(i+1));
                    else
                        lx_end = sig_points(seg_ind(i+1)-1);
                    end
                    line([lx_start lx_end], [lims_y(2) lims_y(2)], ...
                        'color', 'k', 'linewidth',2)
                end
            end
            ylim([lims_y(1) lims_y(2)+lims_y(2)/5])
            xlim([-0.2 1])
            line([t0 t0],ylim, 'color','k','linestyle','--')
            
            path = [mainpath '/erp_data/erp_figures/roi_plots/'];
            str = sprintf('%s_ROI_%s_%s_%s_%s_bc%s_%s',subj_str,roi, ...
                conds{1},conds{2},type,bc,lock);
            %print(gcf,'-loose','-dtiff',[path str])
            hold off
        end
    end
end

%% do statistics in group data
x=[];
subj_str = 'allsubj';
figure;
title_str = sprintf('%s ROI Cleaned ERP: ',rois{1});
for c = 1:numel(conds)
    x{c} = squeeze(x_mean(:,c,:))';
    %x_groupmean = smooth(mean(x{c},1),smfactor,'rloess');
    %x_groupmean = smooth(mean(x{c},1),smfactor);
    x_groupmean = mean(x{c},1);
    x_std = std(x{c},[],1);
    x_sem = (x_std)/sqrt(sum(cond));
    
    shadedErrorBar(trial_time(10:end-10),x_groupmean(10:end-10),x_sem(10:end-10),...
        {[colors{c} '-'],'LineWidth',2,'markerfacecolor',colors{c}},1);
    
    hold on
    title_str = [title_str ' ' colors{c} ' is ' conds{c}];
end


title([subj_str title_str])
xlabel('Time(s)')
if strcmp(type,'ztrialdata')
    ylabel('zscore')
else
    ylabel('\mu Volts')
end

axis tight
line(xlim,[0 0], 'color','k','linestyle','--')
lims_y = ylim;

t0 = trial_time(trial_time==0);

lims_y(2) = lims_y(2)+lims_y(2)/5;
ylim([lims_y(1) lims_y(2)+lims_y(2)/5])
xlim([-0.2 1])
line([t0 t0],ylim, 'color','k','linestyle','--')

binsize = 100;
sldwin = 25;
[~, nsamps] = size(x{1});
nsubj = numel(subidx);

% corresponding time points to the data
time_points = round(linspace(data.dur(1)*1000,data.dur(2)*1000-1,nsamps));

binpoints = time_points(1):sldwin:time_points(end)+1;
% start point of each bin
binStartsamps = binpoints(binpoints<=(time_points(end)+1- binsize));

% end point of each bin
binEndsamps = binpoints(binpoints>=(time_points(1) + binsize));

% for reference, this are the processed bins are stored
bins = [binStartsamps' binEndsamps'];

% number of bins in the classification window
nbins = size(bins,1);

[~,timebin_allocation]=histc(time_points,binpoints);

% create matrix with reduced data size
bin_trials1 = zeros(nsubj,nbins,'single');
bin_trials2 = zeros(nsubj,nbins,'single');

for bin = 1:nbins
    bin_logical = timebin_allocation==bin |timebin_allocation==bin+1;
    bin_trials1(:,bin) = nanmean(x{1}(subidx,bin_logical),2);
    bin_trials2(:,bin) = nanmean(x{2}(subidx,bin_logical),2);
end

% test for signficance
H = ttest(bin_trials1,bin_trials2);
nsegments = sum(H);
seg_idx = find(H);

if nsegments>0
    for i = 1:nsegments
        line([bins(seg_idx(i),:)/1000], [lims_y(2) lims_y(2)], ...
            'color', 'k', 'linewidth',2)
    end
end

path = [mainpath '/erp_data/erp_figures/roi_plots/'];
str = sprintf('%s_ROI_%s_%s_%s_%s_bc%s_%s',subj_str,roi, ...
    conds{1},conds{2},type,bc,lock);
print(gcf,'-loose','-dtiff',[path str])
hold off
%%


