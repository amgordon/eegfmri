function [out] = EF_x_validation(X,Y,S,lambda,predictionType,classifier)
% cross validation index
cv_tr = shuffle(ceil(S.nXvals*(1:size(X,1))/size(X,1)));

yPred = nan(size(Y));
yProb = nan(size(Y));
%penalty range
for l = 1:length(lambda)
    thisL = lambda(l);
    thisLStr = num2str(thisL);
    
    for i=1:S.nXvals
        idxTR = (cv_tr~=i);
        idxTE = (cv_tr==i);
        
        XTrain = X(idxTR, :);
        YTrain = Y(idxTR);
        
        XTest = X(idxTE, :);
        YTest = Y(idxTE);
        
        %feature selection based on mutual information
        mi = nan(size(XTrain,2),1);
        if S.FS
            for t = 1:size(XTrain,2)
                mi(t) = mutualinfo(YTrain, XTrain(:,t));
            end
            
            mi_sorted = sort(mi, 'descend');
            thresh = mi_sorted(S.nFeats);
            XTrain = XTrain(:,mi>=thresh);
            XTest = XTest(:,mi>=thresh);
        end
        
        if strcmp(classifier, 'svr')
            % svr
            trainOpts = [S.trainOptsSVM thisLStr];
            model = svmtrain(YTrain, XTrain, trainOpts);
            [yPred(idxTE),~,yProb(idxTE)]  = svmpredict(YTest, XTest, model);
        elseif strcmp(classifier, 'ridge')
            %ridge regression
            trainOpts = [S.trainOptsSVM thisLStr];
            b = ridge(YTrain,XTrain,thisL, 0);
            yPred(idxTE) = [ones(size(XTest,1),1) XTest]*b;
        elseif strcmp(classifier, 'liblinear')
            trainOpts = [S.trainOptsLibLinear thisLStr];
            model = train(YTrain, XTrain, trainOpts);
            [yPred(idxTE),~,yProb(idxTE)] = predict(YTest, XTest, model);
        end
    end
    
%     if strcmp(classifier, 'svr')
%         model = svmtrain(YTrain, XTrain, trainOpts);
%     elseif strcmp(classifier, 'liblinear')
%         model = train(Y, X, trainOpts);
%     end
    
    out.mod.model = model;
    out.mod.Y = Y;
    out.mod.YPred = yPred;
    out.mod.YProb = yProb;

    if strcmp(predictionType, 'continuous')
        [out.rObsVsPred(l), out.pObsVsPred(l)] = corr(Y,yPred);
    elseif strcmp(predictionType, 'discrete')
        accuracy= Y==yPred;
        [out.meanAccuracyOverClasses(l)] = mean([mean(accuracy(Y==1)) mean(accuracy(Y==-1))]);
    end
end
end