% plot conditions of interest, shaded error bar and two sample ttest to
% determine significance between 2 conditions


clear all
type = 'ztrial';
bc='0';
lock = 'RT';
mainpath = '/Volumes/EXT_HD/ON_eegfmri/';
load([mainpath '/group_data/STGA_STPA_groupdata_' lock '_' type '_bc' bc]);

%%
allconds = {'old','new','hits','cr','Rem_hits','HC_hits','LC_hits','HC_cr',...
    'LC_cr','FA','misses'};

mainpath = '/Volumes/EXT_HD/ON_eegfmri/';
type = 'ztrial';
bc='0';
close all
conds = {'old', 'new'};
fs = 500;

if strcmp(lock,'RT')
    trial_time = -1:1/fs:1;
else
    trial_time = -0.5:1/fs:1.5;
end

area = 'parietal';
rois = {'LPS'};%fieldnames(group.chnls.parietal)';
ch_analysis = 0;
smfactor = 20;
colors = {'r','b','g','k','m','c','y'};
allsubj  = [4,6,8:13];
subj = [4,6,8:9,10:13];
nsubj = 8;
ntrialspersubj = 400;

strial_all = false(nsubj*ntrialspersubj,1);

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
    
    for r = rois
        roi = r{1};
        chnls = data.chnls.(area).(roi);
        
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
                cond = logical(data.trial_type.(conds{c}))&strials;
                
                x{c} = squeeze(mean(data.(type)(chnls,cond,:),1));
                x_mean = smooth(mean(x{c},1),smfactor);
                x_std = std(x{c},[],1);
                x_sem = (x_std)/sqrt(sum(cond));
                
                shadedErrorBar(trial_time(10:end-10),x_mean(10:end-10),x_sem(10:end-10),...
                    {[colors{c} '-'],'LineWidth',2,'markerfacecolor',colors{c}},1);
                hold on
                title_str = [title_str ' ' colors{c} ' is ' conds{c}];
                
            end
            title([subj_str title_str])
            xlabel('Time(s)')
            ylabel('\mu Volts')
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
                ylim([lims_y(1) lims_y(2)+lims_y(2)/5])
                line([t0 t0],ylim, 'color','k','linestyle','--')
            end
            path = [mainpath 'group_data/plots_ROI/'];
            str = sprintf('%s_ROI_%s_%s_%s_%s_bc%s_%s',subj_str,roi, ...
                conds{1},conds{2},type,bc,lock);
            print(gcf,'-loose','-dtiff',[path str])
            hold off
        end
    end
end


