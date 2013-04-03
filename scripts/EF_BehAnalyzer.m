function [res idx] = EF_BehAnalyzer(par)


if strcmp(par.substr, 'ef_083111')
    d = dir(fullfile(par.behavdir, '*respsRecovered.mat'));
else
    d = dir(fullfile(par.behavdir, 'Acc*retrieve*.mat'));
end

oldNew = [];
item = [];
stimResp = [];
judgeResp = [];
stimRT = [];
judgeRT = [];
onset = [];
allTrials = [];
sess = [];

for i = 1:length(par.usedVols)
    R.sub(i).dat = load(fullfile(par.behavdir, d(i).name));
    oldNew  = [oldNew R.sub(i).dat.theData.oldNew];
    item  = [item; R.sub(i).dat.theData.item];
    stimResp  = [stimResp R.sub(i).dat.theData.stimresp];
    judgeResp  = [judgeResp R.sub(i).dat.theData.judgeResp];
    stimRT  = [stimRT R.sub(i).dat.theData.stimRT];
    judgeRT  = [judgeRT R.sub(i).dat.theData.judgeRT];
    onset = [onset R.sub(i).dat.theData.onset - par.dropvol*par.TR];
    allTrials = [allTrials R.sub(i).dat.theData.onset - par.dropvol*par.TR + par.TR*sum(par.usedVols(1:i-1))];
    sess = [sess; i*ones(size(R.sub(i).dat.theData.item))];
end

%the first run for this subject ended before all 64 trials were collected.
%omit the trials collected after the scanner had stopped running.
if strcmp(par.substr, 'ef_071411')
    oldNew(50:64) = [];
    item(50:64) = [];
    stimResp(50:64) = [];
    judgeResp(50:64) = [];
    stimRT(50:64) = [];
    judgeRT(50:64) = [];
    onset(50:64) = [];
    allTrials(50:64) = [];
    sess(50:64) = [];
end

if strcmp(par.substr, 'ef_083111')
    allTrials(end-65:end) = R.sub(2).dat.theData.onset - par.dropvol*par.TR + par.TR*sum(par.usedVols(1:3-1));
end

for i =1:length(stimResp)
    if strcmp(stimResp{i}, 'noanswer')
        recordedResp = judgeResp{i};
        recordedRT = 1+judgeRT{i};
    else
        recordedResp = stimResp{i};
        recordedRT = stimRT{i};
    end
    
    if iscell (recordedResp(1))
        firstResp{i} = recordedResp{1}(1);
        
    else
        firstResp{i} = recordedResp(1);
        
    end
    
    firstRT(i) = recordedRT(1);
end

idx.RT = firstRT;
idx.RT(idx.RT==1) = 0;

firstRespLengths = cellfun(@length, firstResp);
if (prod(firstRespLengths)~=1)
   error('a firstResp element is not of length 1') 
end

idx.onsets = onset;
idx.allTrials = allTrials;
idx.sess = sess;

idx.firstResp = firstResp;

if ismember(par.substr, {'ef_071411' 'ef_111512' 'ef_111412_1' 'ef_111412_2' 'ef_111112_1' 'ef_111112_2' 'ef_111112_3'} )
    idx.respNew = ismember(firstResp, {'6' '7'});
    idx.respOld = ismember(firstResp, {'4' '8' '9' 'a'});
    
    idx.highConf = ismember(firstResp, {'4' '6' '9' 'a'});
    idx.lowConf = ismember(firstResp, {'7' '8'});
    
    idx.old = (oldNew == 2);
    idx.new = (oldNew == 1);
    
    idx.recollect = ismember(firstResp, {'4' 'a'});
    
elseif strcmp(par.substr, 'ef_072111')
    idx.respNew = ismember(firstResp, {'1' '2'});
    idx.respOld = ismember(firstResp, {'3' '4' '9'});
    
    idx.highConf = ismember(firstResp, {'1' '4' '9'});
    idx.lowConf = ismember(firstResp, {'2' '3'});
    
    idx.old = (oldNew == 2);
    idx.new = (oldNew == 1);
    
    idx.recollect = ismember(firstResp, {'9'});
    
elseif strcmp(par.substr, 'ef_083111')
    idx.respNew = ismember(firstResp, {'3' '4'});
    idx.respOld = ismember(firstResp, {'1' '2' '7'});
    
    idx.highConf = ismember(firstResp, {'1' '4' '7'});
    idx.lowConf = ismember(firstResp, {'2' '3'});
    
    idx.old = (oldNew == 2);
    idx.new = (oldNew == 1);
    
    idx.recollect = ismember(firstResp, {'7'});
else
    idx.respNew = ismember(firstResp, {'5' '4'});
    idx.respOld = ismember(firstResp, {'1' '2' '3'});
    
    idx.highConf = ismember(firstResp, {'1' '2' '5'});
    idx.lowConf = ismember(firstResp, {'3' '4'});
    
    idx.old = (oldNew == 1);
    idx.new = (oldNew == 2);
    
    idx.recollect = ismember(firstResp, {'1'});
end

idx.item = item;
idx.noResp = strcmp(firstResp, 'n');
idx.fix = strcmp(item, '+')';
idx.cor = idx.old .* idx.respOld + idx.new .* idx.respNew;

idx.legitResp = idx.respOld + idx.respNew;
idx.memoryTrials = idx.old + idx.new;

idx.conf = 1*(idx.highConf.*idx.respNew) + 2*(idx.lowConf.*idx.respNew) ...
    + 3*(idx.lowConf.*idx.respOld) + 4*(idx.highConf.*idx.respOld.*~idx.recollect) + 5*(idx.recollect);

res.pctCor = sum(idx.cor) / sum(idx.legitResp);
res.pctLegit = sum(idx.legitResp) / sum(idx.memoryTrials);

res.pctCorRem = sum(idx.recollect .* idx.cor ) / sum(idx.recollect);
res.pctCorHiConfOld = sum(idx.highConf .* idx.cor .* idx.respOld) / sum(idx.respOld .* idx.highConf);
res.pctCorLoConfOld = sum(idx.lowConf .* idx.cor .* idx.respOld) / sum(idx.respOld .* idx.lowConf);
res.pctCorLoConfNew = sum(idx.highConf .* idx.cor .* idx.respNew) / sum(idx.respNew .* idx.highConf);
res.pctCorHiConfNew = sum(idx.lowConf .* idx.cor .* idx.respNew) / sum(idx.respNew .* idx.lowConf);

res.pctOldResp = sum(idx.respOld) / sum(idx.legitResp);