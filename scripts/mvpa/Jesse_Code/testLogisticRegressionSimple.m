% 
% In:
% - examples
% - labels
% - scratchpad
%
% Out:
% - scratchpad
% - acts
%

function [acts scratchpad] = testLogisticRegressionSimple( trainpats, traintargs, scratchpad )

examples = trainpats'; nExamples = size(examples,1);
acts = exp( examples * scratchpad.weights + repmat(scratchpad.biases,nExamples,1) )';

