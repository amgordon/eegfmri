function onsets = AG_mvpa_onsets_maker(S)

yTrain = load(fullfile(S.univa_dir.(S.thisTrain), 'SPM'));
yTest = load(fullfile(S.univa_dir.(S.thisTest), 'SPM'));

% yEnc = load(fullfile(univar_dir.Enc, 'SPM'));
% 
% yRet = load(fullfile(univar_dir.Ret, 'SPM'));
% 
% yLoc = load(fullfile(univar_dir.Loc, 'SPM'));
% 
% yRS = load(fullfile(univar_dir.RS, 'SPM'));

encDur = length(yTrain.SPM.xY.VY)*2; %for training on enc, testing on ret

%this is only because the first and second regressors of each task happen
%to be that which we want to train and test.  modify this if you want to
%train and test on other things...
onsets = {[yTrain.SPM.Sess.U(1).ons ; yTest.SPM.Sess.U(1).ons + encDur], ...
    [yTrain.SPM.Sess.U(2).ons ; yTest.SPM.Sess.U(2).ons + encDur]};

%locDur = length(yLoc.SPM.xY.VY)*2; %for training on loc, testing on enc/ret

%-1 to allign the retrieval trials (whose onsets begin at t = 0 in the
%trial with the encode trials (whose onsets begin at t = 1)

%OnsP = [yEnc.SPM.Sess.U(1).ons ; yRet.SPM.Sess.U(1).ons + encDur+1]; %for training on enc, testing on ret

%OnsS = [yEnc.SPM.Sess.U(2).ons ; yRet.SPM.Sess.U(2).ons + encDur+1]; %train on enc, test on ret


%OnsP = [yEnc.SPM.Sess.U(1).ons]; %train on enc, test on enc

%OnsS = [yEnc.SPM.Sess.U(3).ons]; %train on enc, test on enc


%OnsP = [yLoc.SPM.Sess.U(1).ons+1]; % for training on loc, testing on loc

%OnsS = [yLoc.SPM.Sess.U(2).ons+1]; % for training on loc, testing on loc

%OnsP = [yLoc.SPM.Sess.U(1).ons+1 ; yEnc.SPM.Sess.U(1).ons + locDur];
%%for training on loc, testing on enc

%OnsS = [yLoc.SPM.Sess.U(2).ons+1 ; yEnc.SPM.Sess.U(3).ons + locDur];
%%train on loc, test on enc

%OnsP = [yLoc.SPM.Sess.U(1).ons+1 ; yRet.SPM.Sess.U(1).ons+locDur+1];
%%for training on loc, testing on ret

%OnsS = [yLoc.SPM.Sess.U(2).ons+1 ; yRet.SPM.Sess.U(3).ons+locDur+1];
%%train on loc, test on ret

%OnsR = [yRS.SPM.Sess.U(1).ons+1];
%%for training on RS, testing on RS

%OnsL = [yRS.SPM.Sess.U(2).ons+1];
%%train on RS, test on RS


%{Person correct, Scene correct}
%onsets = {OnsP, OnsS};

end





