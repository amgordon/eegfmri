% Train a K class logistic regression classifier, with optional regularization via L2 norm of weight vector(s)
%
% In:
% - examples - #examples x #voxels
% - labels   - #examples x 1
% - optional (see examples below for usage guidelines)
%   - 'lambda',<value> - scalar (if 0, no regularization, default is 1)
%   - 'lambda','crossvalidation',<vector of lambda values>
%   - 'labelsGroup',<labelsGroup> (#examples x 1, useful together with 'lambda','crossvalidation')
%  
%
% Out:
% - a model structure, with fields
%   - W - a 1+#features x #classes matrix (first row contains bias terms for each class)
%   - weights - #features x #classes matrix (same as W(2:end,:))
%   - biases  - 1 x #classes vector (same as W(1,:))

% Dependencies:
% - Carl Rasmussen's minimize.m function (an adapted version bundled below)
% - classifierLogisticRegressionSimple_gradients.m
%   (contains the functions computing gradients, necessary because minimize requires a function
%    to produce the objective value and gradients for the current weight guess)
%
% Notes:
% - learns weights W to maximize the log likelihood of the data
% - the L2 penalty is -0.5 * lambda * sum_#classes_c norm(W(2:end,c))^2
%   (the square of the L2 norm of the weight vector for each class, added over classes)
%
% History:
% - 2009 June - created from previous code - fpereira@princeton.edu 
%
% Examples:
%
%   [model] = classifierLogisticRegressionMulti(examples,labels);
%   [model] = classifierLogisticRegressionMulti(examples,labels,'lambda',10); % L2 lambda = 10
%   [model] = classifierLogisticRegressionMulti(examples,labels,'lambda',0);  % no regularization
%
% Cross-validates within the training set (leave-one-example-of-each-class-out) to decide on
% which lambda value to use (slow)
%   [model] = classifierLogisticRegressionMulti(examples,labels,'lambda','crossvalidation',[1,0.1,10,0.01,100,0.001,1000]
%
% Cross-validates within the training set as before but uses the group label of every example to
% do leave-one-group-out cross-validation within the training set (faster). The groups could be
% runs, for instance, or any kind of division of the data, as long as examples that are very
% close together in time (i.e. the hemodynamic responses overlap, they are in the same block)
% are kept in the same group.
%
%   [model] = classifierLogisticRegressionMulti(examples,labels,'lambda','crossvalidation',[1,0.1,10,0.01,100,0.001,1000],'labelsGroup',labelsGroup);


function [scratchpad] = trainLogisticRegressionSimple(trainpats,traintargs,class_args,cv_args)

%
% process arguments
%

this = 'trainLogisticRegressionSimple';

%% process arguments
                         
% default parameters
defaults.regularization = 'L2';
defaults.lambda         = 1;
defaults.optimization   = 'minimize';
defaults.labelsGroup    = []; 


if 1
  % normal operation
  
  class_args = mergestructs(class_args, defaults);

  regularization = class_args.regularization;
  lambda         = class_args.lambda;
  optimization   = class_args.optimization;
else
  % debug
  regularization = defaults.regularization;
  lambda         = defaults.lambda;
  optimization   = defaults.optimization;
end
  
%% call the classifier function

model = classifierLogisticRegressionSimple( trainpats', traintargs', 'lambda', lambda );

%% pack the results

scratchpad.W = model.W;
scratchpad.weights = model.weights;
scratchpad.biases  = model.biases;

scratchpad.lambda         = lambda;
scratchpad.regularization = regularization;
scratchpad.optimization   = optimization;
