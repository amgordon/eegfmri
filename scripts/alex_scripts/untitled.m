%his script puts all the results into a single file

subjs = { 'ef_040412'  'ef_040512' 'ef_040712' 'ef_040712_2'  'ef_041112' 'ef_042912' 'ef_050112'};
nsubjs = numel(subjs);
rois = channel_rois();
rois = rois.roinames;
nrois = numel(rois);
perfmetric = 'bac';
conds = {'hits','cr'};
solver = 'glmnetsolveralpha0p5';
dur = '0_1';
binsize = 'binsize50';

datapath = 'biac4/wagner/biac3/wagner5/alan/eegfmri/';

perfmat = zeros(nsubj,nrois);

for s = 1: nsubjs
    filepath = [datapath subjs{s} '/erpData/classificationResults/'];
    for r = 1: nrois
         filename = [subjs{s} '_' conds(1) '_' conds(2) '_' solver  ...
             '_time' dur '_' rois{r} '_' binsize '.mat']
         load([filepath filename])
         perfmat(s,r)=param_model.maxperf_BAC;
         
    end 
end