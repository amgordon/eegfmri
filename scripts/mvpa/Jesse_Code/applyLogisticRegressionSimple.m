%
% Apply a classifier learned with trainLogisticRegressionSimple to a new dataset
% 
% In:
% - examples - #examples x #voxels
% - model    - the classifier returned by trainLogisticRegressionSimple
%
% Out:
% - predictedLabels - #examples x 1 - the labels predicted
%
% History
% 2009 June - created - fpereira@princeton.edu
%

function [predictedLabels] = applyLogisticRegressionSimple( examples, model )

[nExamples,nFeatures] = size(examples);

sortedLabelValues = model.sortedLabelValues;

decisions = examples * model.weights + repmat(model.biases,nExamples,1);
[maxval,maxdec] = max(decisions,[],2);

predictedLabels = sortedLabelValues(maxdec);
