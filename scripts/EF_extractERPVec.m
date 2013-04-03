function [ erpidx ] = EF_extractERPVec(par)

eeg = load(par.erpFile);

channel_rois;

eegDat_h = eeg.data.trialdata;

%% define bad channels and trials
channel_rois;
S.chnls = chnls;
clear chnls;

eegfmri_ON_trialch_info

S.goodChannels = true(eeg.data.Nchan,1);
S.goodChannels(S.badch{par.subNo}) = false;

%% channel and timepoint selection

thisChan = ['S.chnls.' par.ChOI];
ChOI_h = eval(thisChan);

ChOI = ChOI_h;
chansToInclude = sort(setdiff(ChOI, S.badch{par.subNo}));

thisTOI = par.TOI;
eegDat = eegDat_h(chansToInclude,:,thisTOI);

erpidx.sigAll = mean(mean(eegDat),3);

% 
% erpidx.sigAll = vertcat(erpidx.sigAll_mean{:});
% erpidx.all = vertcat(erpidx.all_h{:})';
% 
% 
% eeg_fmri_on_subject_info;
erpidx.all = ones(size(idx.allTrials));
erpidx.all(S.badtr{par.subNo}) = 0;


if strcmp(par.substr, 'ef_072111')
    erpidx.all = [erpidx.all 0];
    erpidx.sigAll = [erpidx.sigAll; 0];
end


end

