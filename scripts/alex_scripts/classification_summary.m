%his script puts all the results into a single file

subjs = { 'ef_040412'  'ef_040512' 'ef_040712' 'ef_040712_2'  'ef_041112' 'ef_042912' 'ef_050112'};
nsubjs = numel(subjs);
rois =  {'LPI','LPS','PM','RPS','RPI','LFI','LFS','FM','RFS','RFI'};% {'allchnls'};
nrois = numel(rois);
perfmetric = 'bac';
conds = {'hits','cr'};
solver = 'glmnetsolveralpha0p5'; % 'glmnetsolveralpha0p5';
dur = '0_1';
binsize = 'binsize50';

datapath = '/biac4/wagner/biac3/wagner5/alan/eegfmri/fmri_data/';

perfmat = zeros(nsubjs,nrois);

for s = 1: nsubjs
    filepath = [datapath subjs{s} '/erpData/classificationResults/'];
    for r = 1: nrois
         filename = [subjs{s} '_' conds{1} '_' conds{2} '_' solver  ...
             '_time' dur '_' rois{r} '_' binsize '.mat'];
         load([filepath filename])
         perfmat(s,r)=param_model.maxperf_BAC;
         
    end 
end

%%

% using all channels
% perfmat =
% 
%     0.5075
%     0.5579
%     0.4896
%     0.5179
%     0.6171
%     0.4797
%     0.4514
% 
% using rois
% 
% perfmat =
% 
%     0.4938    0.4073    0.5449    0.5373    0.4648    0.5118    0.5225    0.4518    0.5131    0.5429
%     0.5308    0.5098    0.4972    0.4190    0.4946    0.4783    0.5676    0.5258    0.4752    0.4644
%     0.5290    0.4604    0.4953    0.5568    0.4881    0.4705    0.4629    0.4793    0.5147    0.5157
%     0.4617    0.4960    0.5821    0.4922    0.5098    0.5090    0.5060    0.5641    0.5460    0.5456
%     0.5259    0.5056    0.4897    0.4983    0.4936    0.5723    0.5105    0.4330    0.5192    0.5067
%     0.5206    0.5775    0.5450    0.5256    0.5206    0.4459    0.4419    0.4919    0.4359    0.4644
%     0.4610    0.5315    0.4904    0.4840    0.4979    0.4562    0.5224    0.5401    0.5080    0.4931

