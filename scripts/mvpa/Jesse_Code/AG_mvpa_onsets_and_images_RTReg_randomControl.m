function [onsets, filenames] = AG_mvpa_onsets_and_images(S)

yTrainSPM = load(fullfile(S.univar_dir.(S.thisTrain), 'SPM'));
yTestSPM = load(fullfile(S.univar_dir.(S.thisTest), 'SPM'));

yTrain = load(fullfile(S.univar_dir.(S.thisTrain), 'mvpa_ons'));
yTest = load(fullfile(S.univar_dir.(S.thisTest), 'mvpa_ons'));

yIDDir = fullfile(S.expt_dir, S.subj_id, 'behav');
yIDHelp = dir([yIDDir '/' 'AG1_retrieve*']);
yIDCheck = load(fullfile(yIDDir, yIDHelp.name));
yRetHand = yIDCheck.S.retHandNum ;

RegsOfInterest = {[1] [1 2] [1 2] [1 2]};

cF = [1 0 0 0];
% yEnc = load(fullfile(univar_dir.Enc, 'SPM'));
%
% yRet = load(fullfile(univar_dir.Ret, 'SPM'));
%
% yLoc = load(fullfile(univar_dir.Loc, 'SPM'));
%
% yRS = load(fullfile(univar_dir.RS, 'SPM'));

encDur = length(yTrainSPM.SPM.xY.VY)*2; %for training on enc, testing on ret

fnTrain = cellstr(yTrainSPM.SPM.xY.P);
fnTest = cellstr(yTestSPM.SPM.xY.P);


%convert from "swascans" to "wascans"

if S.use_unsmoothed
    for i = 1:length(fnTrain)
        fnTrain{i}(80) = [];
    end
    
    for i = 1:length(fnTest)
        fnTest{i}(80) = [];
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Standard MVPA Stuff
%---------------------
%

onsets_h= load(fullfile(S.rtReg_dir, 'mvpa_ons'));
onsets = {vertcat(onsets_h.onsets{S.onsetsIdx})};

  [useless_var, ix_randomizer] = sort(rand(1,length(onsets{1})));
  onsets{1} = onsets{1}(ix_randomizer);

trainRegs = RegsOfInterest{S.idxThisTrain};
testRegs = RegsOfInterest{S.idxThisTest};



if ((S.idxThisTrain == S.idxThisTest)&&(strcmp(S.thisTest,'Enc')))
    

    filenames = fnTrain;
    
%     onsets = {[yTrain.onsets{trainRegs(1)} - cF(S.idxThisTrain)], ...
%         [yTrain.onsets{trainRegs(2)} - cF(S.idxThisTrain)]};
    
    
elseif  ((S.idxThisTrain == S.idxThisTest))
%     onsets = {[yTrain.onsets{trainRegs(1)} - cF(S.idxThisTrain)], ...
%         [yTrain.onsets{trainRegs(2)} - cF(S.idxThisTrain)]};
    
    filenames = fnTrain;
    
elseif ((strcmp(S.thisTrain, 'Enc'))&&(strcmp(S.thisTest,'Ret')))
%     onsets = {[yTrain.onsets{trainRegs(1)} - cF(S.idxThisTrain); yTest.onsets{testRegs(1)} + encDur - cF(S.idxThisTest)], ...
%         [yTrain.onsets{trainRegs(2)} - cF(S.idxThisTrain); yTest.onsets{testRegs(2)} + encDur - cF(S.idxThisTest)]};
    filenames = vertcat(fnTrain, fnTest);
    
%  elseif ((strcmp(S.thisTrain, 'Enc'))&&(strcmp(S.thisTest,'Ret')))
%      onsets = {[yTrain.onsets{trainRegs(1)} - cF(S.idxThisTrain);  yTest.onsets{3} + encDur - cF(S.idxThisTest)], ...
%          [yTrain.onsets{trainRegs(2)} - cF(S.idxThisTrain); yTest.onsets{4} + encDur - cF(S.idxThisTest)]};
%      filenames = vertcat(fnTrain, fnTest);
%     
    
elseif ((yRetHand==1)&&(strcmp(S.thisTrain, 'respSel'))&&(strcmp(S.thisTest,'Ret')))
%     onsets = {[yTrain.onsets{trainRegs(1)} - cF(S.idxThisTrain); yTest.onsets{testRegs(1)} + encDur - cF(S.idxThisTest)], ...
%         [yTrain.onsets{trainRegs(2)} - cF(S.idxThisTrain); yTest.onsets{testRegs(2)} + encDur - cF(S.idxThisTest)]};
    filenames = vertcat(fnTrain, fnTest);
    
    
elseif ((yRetHand==2)&&(strcmp(S.thisTrain, 'respSel'))&&(strcmp(S.thisTest,'Ret')))
%     onsets = {[yTrain.onsets{trainRegs(1)} - cF(S.idxThisTrain); yTest.onsets{testRegs(2)} + encDur - cF(S.idxThisTest)], ...
%         [yTrain.onsets{trainRegs(2)} - cF(S.idxThisTrain); yTest.onsets{trainRegs(1)} + encDur - cF(S.idxThisTest)]};
    filenames = vertcat(fnTrain, fnTest);
    
else
%     onsets = {[yTrain.onsets{trainRegs(1)} - cF(S.idxThisTrain); yTest.onsets(testRegs(1)) + encDur - cF(S.idxThisTest)], ...
%         [yTrain.onsets{trainRegs(2)} - cF(S.idxThisTrain); yTest.onsets(testRegs(2)) + encDur - cF(S.idxThisTest)]};
    
    filenames = vertcat(fnTrain, fnTest);
end




%this is only because the first and second regressors of each task happen
%to be that which we want to train and test.  modify this if you want to
%train and test on other things...
% if S.idxThisTrain == S.idxThisTest
%
%
%
%     onsets = {[yTrain.SPM.Sess.U(1).ons - cF(S.idxThisTrain)], ...
%         [yTrain.SPM.Sess.U(2).ons - cF(S.idxThisTrain)]};
%
%     filenames = fnTrain;
% elseif ((yRetHand==2)&&(strcmp(S.thisTrain, 'respSel'))&&(strcmp(S.thisTest,'Ret')))
%     onsets = {[yTrain.SPM.Sess.U(1).ons - cF(S.idxThisTrain); yTest.SPM.Sess.U(2).ons + encDur - cF(S.idxThisTest)], ...
%         [yTrain.SPM.Sess.U(2).ons - cF(S.idxThisTrain); yTest.SPM.Sess.U(1).ons + encDur - cF(S.idxThisTest)]};
%     filenames = vertcat(fnTrain, fnTest);
% else
%     onsets = {[yTrain.SPM.Sess.U(1).ons - cF(S.idxThisTrain); yTest.SPM.Sess.U(1).ons + encDur - cF(S.idxThisTest)], ...
%         [yTrain.SPM.Sess.U(2).ons - cF(S.idxThisTrain); yTest.SPM.Sess.U(2).ons + encDur - cF(S.idxThisTest)]};
%
%     filenames = vertcat(fnTrain, fnTest);
% end



%Testing on all trials
%----------------------
%
% onsets = [];
%
% sceneAndJunkTrials = sort([yTest.SPM.Sess.U(1).ons; yTest.SPM.Sess.U(2).ons; yTest.SPM.Sess.U(3).ons]);
%
% onsets = {[yTrain.SPM.Sess.U(1).ons - cF(S.idxThisTrain); sceneAndJunkTrials + encDur - cF(S.idxThisTest)], ...
%         [yTrain.SPM.Sess.U(2).ons - cF(S.idxThisTrain)]};










