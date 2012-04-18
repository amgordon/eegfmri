analysisdir = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/ef_072111/analysis_ON';
behavdir = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/ef_072111/behav';

clear idx run1 run2 run3 allons allcond allresp onsets pmod sessreg p1 p2 p3 p4;

oldResp_HC = {'4$' '9('};
newResp_HC = {'1!' };
oldResp_LC = {'3#'};
newResp_LC = {'2@' };

%49 trials for run1;
%eegDat = load('LPS_zavg_bins');

cd (behavdir);

yrun1 = load('Acc1_retrieve_1_21Jul11out(1).mat');

run1.cond = yrun1.theData.oldNew;
run1.resp = yrun1.theData.stimresp;
run1.resp(strcmp(yrun1.theData.stimresp, 'noanswer'))=yrun1.theData.judgeResp(strcmp(yrun1.theData.stimresp, 'noanswer'));
run1.ons = yrun1.theData.StimulusOnseTime - yrun1.startTime;
%run1.goodtr = eegDat.zsubj(2).run(1).new.goodtr + eegDat.zsubj(2).run(1).old.goodtr;

yrun2 = load('Acc1_retrieve_1_21Jul11out(2).mat');

run2.cond = yrun2.theData.oldNew;
run2.resp = yrun2.theData.stimresp;
run2.resp(strcmp(yrun2.theData.stimresp, 'noanswer'))=yrun2.theData.judgeResp(strcmp(yrun2.theData.stimresp, 'noanswer'));
run2.ons = yrun2.theData.StimulusOnseTime - yrun2.startTime;
%run2.goodtr = [eegDat.zsubj(2).run(2).new.goodtr + eegDat.zsubj(2).run(2).old.goodtr; 0];

yrun3 = load('Acc1_retrieve_1_21Jul11out(4).mat');

run3.cond = yrun3.theData.oldNew;
run3.resp = yrun3.theData.stimresp;
run3.resp(strcmp(yrun3.theData.stimresp, 'noanswer'))=yrun3.theData.judgeResp(strcmp(yrun3.theData.stimresp, 'noanswer'));
run3.ons = yrun3.theData.StimulusOnseTime - yrun3.startTime;
%run3.goodtr = eegDat.zsubj(2).run(4).new.goodtr + eegDat.zsubj(2).run(4).old.goodtr;

yrun4 = load('Acc1_retrieve_1_21Jul11out(5).mat');

run4.cond = yrun4.theData.oldNew;
run4.resp = yrun4.theData.stimresp;
run4.resp(strcmp(yrun4.theData.stimresp, 'noanswer'))=yrun4.theData.judgeResp(strcmp(yrun4.theData.stimresp, 'noanswer'));
run4.ons = yrun4.theData.StimulusOnseTime - yrun4.startTime;
%run4.goodtr = eegDat.zsubj(2).run(5).new.goodtr + eegDat.zsubj(2).run(5).old.goodtr;

allons = [(run1.ons) (run2.ons+368) (run3.ons+2*368) (run4.ons + 3*368)];
allcond = [run1.cond run2.cond run3.cond run4.cond];
allresp = [run1.resp run2.resp run3.resp run4.resp];
%idx.allgoodtr = [run1.goodtr; run2.goodtr; run3.goodtr; run4.goodtr ]';

for i = 1:length(allresp)
    idx.respOld_HC(i) = ismember(allresp(i), oldResp_HC);
    idx.respNew_HC(i) = ismember(allresp(i), newResp_HC);
    idx.respOld_LC(i) = ismember(allresp(i), oldResp_LC);
    idx.respNew_LC(i) = ismember(allresp(i), newResp_LC);
    
end

idx.old = allcond==2;
idx.new = allcond==1;

idx.respOld = idx.respOld_HC + idx.respOld_LC;
idx.respNew = idx.respNew_HC + idx.respNew_LC;
% 
%  onsets{1} = allons(find(idx.old .* idx.respOld_HC .* idx.allgoodtr));
%  onsets{2} = allons(find(idx.new.* idx.respNew_HC.* idx.allgoodtr));
%  onsets{3} = allons(find(idx.old .* idx.respOld_LC.* idx.allgoodtr));
%  onsets{4} = allons(find(idx.new.* idx.respNew_LC.* idx.allgoodtr));
%  onsets{5} = allons(find(idx.old .* ~idx.respOld_HC .* ~idx.respOld_LC + idx.new.* ~idx.respNew_HC .* ~idx.respNew_LC + ~idx.allgoodtr));

onsets{1} = allons(find(idx.old .* idx.respOld));
onsets{2} = allons(find(idx.new .* idx.respNew));
onsets{3} = allons(find(idx.old .* ~idx.respOld + idx.new .* ~idx.respNew));

subidx_h{1} = find(idx.old .* idx.respOld_HC);
subidx_h{2} = find(idx.new .* idx.respNew_HC);
subidx_h{3} = find(idx.old .* idx.respOld_LC);
subidx_h{4} = find(idx.new .* idx.respNew_LC);
 
%subidx{1} = find(idx.allgoodtr(subidx_h{1}));
%subidx{2} = find(idx.allgoodtr(subidx_h{2}));
%subidx{3} = find(idx.allgoodtr(subidx_h{3}));
%subidx{4} = find(idx.allgoodtr(subidx_h{4}));

% names = {'hits_HC' 'CRs_HC' 'hits_LC' 'CRs_LC' 'junk'};
names = {'hits' 'CRs' 'junk'};

durations = {0 0 0};

goodChannels=[6 8 9 10 11 12 16 17];

%for i = [1 2 4 5], p1{i} = squeeze(mean(eegDat.zsubj(2).run(i).HChit.bins(goodChannels,6,:))); end
%for i = [1 2 4 5], p2{i} = squeeze(mean(eegDat.zsubj(2).run(i).HCcr.bins(goodChannels,6,:))); end
%for i = [1 2 4 5], p3{i} = squeeze(mean(eegDat.zsubj(2).run(i).hit.bins(goodChannels,6,:))); end
%for i = [1 2 4 5], p4{i} = squeeze(mean(eegDat.zsubj(2).run(i).cr.bins(goodChannels,6,:))); end


% param_h{1} = [vertcat(p1{:})];
% param_h{2} = [vertcat(p2{:})];
% param_h{3} = [setxor(vertcat(p3{:}), vertcat(p1{:}))];
% param_h{4} = [setxor(vertcat(p4{:}), vertcat(p2{:}))];

% pmod(1).name = {'by_erp_amplitude'};
% pmod(1).param = {param_h{1}(subidx{1})};
% pmod(1).poly = {1};
% 
% pmod(2).name = {'by_erp_amplitude'};
% pmod(2).param = {param_h{2}(subidx{2})};
% pmod(2).poly = {1};
% 
% pmod(3).name = {'by_erp_amplitude'};
% pmod(3).param = {param_h{3}(subidx{3})};
% pmod(3).poly = {1};
% 
% pmod(4).name = {'by_erp_amplitude'};
% pmod(4).param = {param_h{4}(subidx{4})};
% pmod(4).poly = {1};

sessReg = zeros(sum(par.numvols),length(par.numvols)-1);
for i = 1:(length(par.numvols)-1)
    sessReg((i-1)*par.numvols(i)+1:i*par.numvols(i),i) = ones(par.numvols(i),1);
end

R = sessReg;

cd (analysisdir);

%save ons onsets names durations pmod;
save ons onsets names durations;
save regs R;