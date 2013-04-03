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
%

function [model] = classifierLogisticRegressionSimple( varargin )

%
% process arguments
%

this = 'classifierLogisticRegressionSimple';

if nargin < 2; eval(sprintf('help %s;',this)); return; else

  examples = varargin{1};
  labels   = varargin{2};
  
  % defaults for options
  regularization = 'L2'; lambda = 1;
  runSilent = 0;
  labelsGroup = [];
  optimizationMethod = 'minimize';
  
  if nargin > 2
    % there are additional arguments to process
  
    idx = 3;
    while idx <= nargin
      argName = varargin{idx}; idx = idx + 1;
      switch argName
       case {'lambda'}
        lambda = varargin{idx}; idx = idx + 1;
        switch lambda
         case {'crossvalidation'}
          % extra arguments are needed
          lambdaRange = varargin{idx}; idx = idx + 1;
         otherwise
          % we're done
        end
       case {'labelsGroup'}
        labelsGroup = varargin{idx}; idx = idx + 1;
       case {'runSilent'}
        runSilent = 1;
       otherwise
	fprintf('%s: unknown argument %s\n',this,argName); pause
      end 
    end
  end

end

%% compute information about labels

sortedLabelValues     = sort(unique(labels)); conditions = sortedLabelValues;
nClasses              = length(sortedLabelValues);
[nExamples,nFeatures] = size(examples);

% find examples in each class and create a binary indicator mask
delta = zeros(nExamples,nClasses);

for c = 1:nClasses
  label        = sortedLabelValues(c);  
  indices{c}   = find(labels == label);
  nperclass(c) = length(indices{c});

  delta(indices{c},c) = 1;
end


%% if requested, run a cross validation to set the value of lambda

if isequal(lambda,'crossvalidation')
  fprintf('%s: using crossvalidation to find lambda\n',this);

  if isempty(labelsGroup)

    % we have no group labels, create fake group labels to run leave-one-example-out
    nGroups = min(nperclass);
    labelsGroup = zeros(nExamples,1);

    idx = 1;
    for c = 1:nClasses
      nc  = nperclass(c);
      tmp = zeros(nc,1);
      idx = 1;
      while nc >= nGroups
        tmp(idx:(idx+nGroups-1)) = 1:nGroups;
        nc  = nc - nGroups;
        idx = idx + nGroups;
      end
      if nc > 0
        tmp(idx:(idx+nc-1)) = 1:nc;
      end
      labelsGroup(indices{c}) = tmp;
    end
  end

  %% set up a leave-one-group-out cross-validation
  
  nl = length(lambdaRange);
  groupValues = unique(labelsGroup); nGroups = length(groupValues);

  for g = 1:nGroups
    group = groupValues(g); mask = (labelsGroup == group);
    indicesTest{g}  = find( mask); nTest(g)  = length(indicesTest{g});
    indicesTrain{g} = find(~mask); nTrain(g) = length(indicesTrain{g});
  end

  %% run it for every value of lamba in the range given
  
  cvresults = zeros(1,nl);
  
  for r = 1:nl
    lambda = lambdaRange(r);
    
    fprintf('\ttesting lambda=%s (%d folds)\n',num2str(lambda),nGroups);
    
    predictedLabels = zeros(nExamples,1);

    for g = 1:nGroups
      examplesTest  = examples(indicesTest{g},:);
      labelsTest    = labels(  indicesTest{g},:);
      examplesTrain = examples(indicesTrain{g},:);
      labelsTrain   = labels(  indicesTrain{g},:);
    
      cvmodel   = classifierLogisticRegressionSimple( examplesTrain,labelsTrain,'lambda',lambda );
      decisions = examplesTest * cvmodel.weights + repmat(cvmodel.biases,nTest(g),1);
      [discard,maxpos] = max(decisions,[],2);
      
      predictedLabels(indicesTest{g}) = sortedLabelValues(maxpos);
    end

    cvresults(r) = sum(predictedLabels == labels)/nExamples;
  end
  
  [maxval,maxpos] = max(cvresults);

  lambda = lambdaRange(maxpos);

  fprintf('%s: crossvalidation estimated lambda is %s\n',this,num2str(lambda));
  
  clear examplesTrain examplesTest labelsTrain labelsTest;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Train classifier
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial guess for the weight matrix
W = randn(nFeatures+1,nClasses) * 0.1;

switch optimizationMethod
 case {'minimize'}
  % use carl rasmussen's minimize.m (renamed to optimize.m to avoid a clash with CVX, appended below)
  [Wvector,Wvalue] = optimize(W(:),'classifierLogisticRegressionMulti_gradients',1000,examples,labels,lambda,indices,delta,regularization,'minimize');    
end
  
% get the objective at the end
[f,df,fpart,lambdapart] = classifierLogisticRegressionSimple_gradients(Wvector,examples,labels,lambda,indices,delta,regularization,'minimize');
if ~runSilent; fprintf('minimize: fpart = %s lambdapart = %s\n',num2str(fpart),num2str(lambdapart));  end
  
W = reshape(Wvector,size(W));
tmp = W(2:end,:);

if ~runSilent; fprintf('L1(W)=%s\tL2(W)=%s\n',num2str(sqrt(sum(tmp(:).^2))),num2str(sum(abs(tmp(:)))));end
obj = -Wvalue; % it returns -1*objective

%% Pack the results into a struct

model.W       = W;
model.biases  = W(1,:);     % to help people
model.weights = W(2:end,:); % decipher the output
model.sortedLabelValues = sortedLabelValues;

% training Set information
trainingSetInfo.nExamples         = nExamples;
trainingSetInfo.nFeatures         = nFeatures;
trainingSetInfo.nLabels           = nClasses;
trainingSetInfo.nClasses          = nClasses; % same thing as labels
trainingSetInfo.sortedLabelValues = sortedLabelValues;
trainingSetInfo.classPriors       = nperclass / sum(nperclass);
model.trainingSetInfo = trainingSetInfo;



% Carl Edward Rasmussen's minimize modified to 
% - not print anything
% - use less stringent convergence criteria
% to be used in inner loops of other code

function [X, fX, i] = optimize(X, f, length, varargin)

% Minimize a differentiable multivariate function. 
%
% Usage: [X, fX, i] = minimize(X, f, length, P1, P2, P3, ... )
%
% where the starting point is given by "X" (D by 1), and the function named in
% the string "f", must return a function value and a vector of partial
% derivatives of f wrt X, the "length" gives the length of the run: if it is
% positive, it gives the maximum number of line searches, if negative its
% absolute gives the maximum allowed number of function evaluations. You can
% (optionally) give "length" a second component, which will indicate the
% reduction in function value to be expected in the first line-search (defaults
% to 1.0). The parameters P1, P2, P3, ... are passed on to the function f.
%
% The function returns when either its length is up, or if no further progress
% can be made (ie, we are at a (local) minimum, or so close that due to
% numerical problems, we cannot get any closer). NOTE: If the function
% terminates within a few iterations, it could be an indication that the
% function values and derivatives are not consistent (ie, there may be a bug in
% the implementation of your "f" function). The function returns the found
% solution "X", a vector of function values "fX" indicating the progress made
% and "i" the number of iterations (line searches or function evaluations,
% depending on the sign of "length") used.
%
% The Polack-Ribiere flavour of conjugate gradients is used to compute search
% directions, and a line search using quadratic and cubic polynomial
% approximations and the Wolfe-Powell stopping criteria is used together with
% the slope ratio method for guessing initial step sizes. Additionally a bunch
% of checks are made to make sure that exploration is taking place and that
% extrapolation will not be unboundedly large.
%
% See also: checkgrad 
%
% Copyright (C) 2001 - 2006 by Carl Edward Rasmussen (2006-09-08).

INT = 0.1;    % don't reevaluate within 0.1 of the limit of the current bracket
EXT = 3.0;                  % extrapolate maximum 3 times the current step-size
MAX = 20;                         % max 20 function evaluations per line search
RATIO = 10;                                       % maximum allowed slope ratio
SIG = 0.1; RHO = SIG/2; % SIG and RHO are the constants controlling the Wolfe-
% Powell conditions. SIG is the maximum allowed absolute ratio between
% previous and new slopes (derivatives in the search direction), thus setting
% SIG to low (positive) values forces higher precision in the line-searches.
% RHO is the minimum allowed fraction of the expected (from the slope at the
% initial point in the linesearch). Constants must satisfy 0 < RHO < SIG < 1.
% Tuning of SIG (depending on the nature of the function to be optimized) may
% speed up the minimization; it is probably not worth playing much with RHO.

% The code falls naturally into 3 parts, after the initial line search is
% started in the direction of steepest descent. 1) we first enter a while loop
% which uses point 1 (p1) and (p2) to compute an extrapolation (p3), until we
% have extrapolated far enough (Wolfe-Powell conditions). 2) if necessary, we
% enter the second loop which takes p2, p3 and p4 chooses the subinterval
% containing a (local) minimum, and interpolates it, unil an acceptable point
% is found (Wolfe-Powell conditions). Note, that points are always maintained
% in order p0 <= p1 <= p2 < p3 < p4. 3) compute a new search direction using
% conjugate gradients (Polack-Ribiere flavour), or revert to steepest if there
% was a problem in the previous line-search. Return the best value so far, if
% two consecutive line-searches fail, or whenever we run out of function
% evaluations or line-searches. During extrapolation, the "f" function may fail
% either with an error or returning Nan or Inf, and minimize should handle this
% gracefully.

if max(size(length)) == 2, red=length(2); length=length(1); else red=1; end
if length>0, S='Linesearch'; else S='Function evaluation'; end 

i = 0;                                            % zero the run length counter
ls_failed = 0;                             % no previous line search has failed
[f0 df0] = feval(f, X, varargin{:});          % get function value and gradient
fX = f0;
i = i + (length<0);                                            % count epochs?!
s = -df0; d0 = -s'*s;           % initial search direction (steepest) and slope
x3 = red/(1-d0);                                  % initial step is red/(|s|+1)
f0prev = f0; % FP

while i < abs(length)                                      % while not finished
  i = i + (length>0);                                      % count iterations?!

  X0 = X; F0 = f0; dF0 = df0;                   % make a copy of current values
  if length>0, M = MAX; else M = min(MAX, -length-i); end

  while 1                             % keep extrapolating as long as necessary
    x2 = 0; f2 = f0; d2 = d0; f3 = f0; df3 = df0;
    success = 0;
    while ~success && M > 0
      try
        M = M - 1; i = i + (length<0);                         % count epochs?!
        [f3 df3] = feval(f, X+x3*s, varargin{:});
        if isnan(f3) || isinf(f3) || any(isnan(df3)+isinf(df3)), error(''), end
        success = 1;
      catch                                % catch any error which occured in f
        x3 = (x2+x3)/2;                                  % bisect and try again
      end
    end
    if f3 < F0, X0 = X+x3*s; F0 = f3; dF0 = df3; end         % keep best values
    d3 = df3'*s;                                                    % new slope
    if d3 > SIG*d0 || f3 > f0+x3*RHO*d0 || M == 0  % are we done extrapolating?
      break
    end
    x1 = x2; f1 = f2; d1 = d2;                        % move point 2 to point 1
    x2 = x3; f2 = f3; d2 = d3;                        % move point 3 to point 2
    A = 6*(f1-f2)+3*(d2+d1)*(x2-x1);                 % make cubic extrapolation
    B = 3*(f2-f1)-(2*d1+d2)*(x2-x1);
    x3 = x1-d1*(x2-x1)^2/(B+sqrt(B*B-A*d1*(x2-x1))); % num. error possible, ok!
    if ~isreal(x3) || isnan(x3) || isinf(x3) || x3 < 0 % num prob | wrong sign?
      x3 = x2*EXT;                                 % extrapolate maximum amount
    elseif x3 > x2*EXT                  % new point beyond extrapolation limit?
      x3 = x2*EXT;                                 % extrapolate maximum amount
    elseif x3 < x2+INT*(x2-x1)         % new point too close to previous point?
      x3 = x2+INT*(x2-x1);
    end
  end                                                       % end extrapolation

  while (abs(d3) > -SIG*d0 || f3 > f0+x3*RHO*d0) && M > 0  % keep interpolating
    if d3 > 0 || f3 > f0+x3*RHO*d0                         % choose subinterval
      x4 = x3; f4 = f3; d4 = d3;                      % move point 3 to point 4
    else
      x2 = x3; f2 = f3; d2 = d3;                      % move point 3 to point 2
    end
    if f4 > f0           
      x3 = x2-(0.5*d2*(x4-x2)^2)/(f4-f2-d2*(x4-x2));  % quadratic interpolation
    else
      A = 6*(f2-f4)/(x4-x2)+3*(d4+d2);                    % cubic interpolation
      B = 3*(f4-f2)-(2*d2+d4)*(x4-x2);
      x3 = x2+(sqrt(B*B-A*d2*(x4-x2)^2)-B)/A;        % num. error possible, ok!
    end
    if isnan(x3) || isinf(x3)
      x3 = (x2+x4)/2;               % if we had a numerical problem then bisect
    end
    x3 = max(min(x3, x4-INT*(x4-x2)),x2+INT*(x4-x2));  % don't accept too close
    [f3 df3] = feval(f, X+x3*s, varargin{:});
    if f3 < F0, X0 = X+x3*s; F0 = f3; dF0 = df3; end         % keep best values
    M = M - 1; i = i + (length<0);                             % count epochs?!
    d3 = df3'*s;                                                    % new slope
  end                                                       % end interpolation

  if abs(d3) < -SIG*d0 && f3 < f0+x3*RHO*d0          % if line search succeeded
    X = X+x3*s; f0 = f3; fX = [fX' f0]';                     % update variables

    change = abs((f0prev-f0)/(f0prev)); % FP
    f0prev = f0; % FP

    if ((i<=10) | (~rem(i,10))) fprintf('%s %6i;  Value %4.10e %f\r', S, i, f0,change); end % FP
    %    fprintf('%s %6i;  Value %4.6e\r', S, i, f0);
    s = (df3'*df3-df0'*df3)/(df0'*df0)*s - df3;   % Polack-Ribiere CG direction
    df0 = df3;                                               % swap derivatives
    d3 = d0; d0 = df0'*s;
    if d0 > 0                                      % new slope must be negative
      s = -df0; d0 = -s'*s;                  % otherwise use steepest direction
    end
    x3 = x3 * min(RATIO, d3/(d0-realmin));          % slope ratio but max RATIO
    ls_failed = 0;                              % this line search did not fail

       
    if change < 0.0001; fprintf('converged with objective percent change %f\n',change); break; end
  else
    X = X0; f0 = F0; df0 = dF0;                     % restore best point so far
    if ls_failed || i > abs(length)         % line search failed twice in a row
      break;                             % or we ran out of time, so we give up
    end
    s = -df0; d0 = -s'*s;                                        % try steepest
    x3 = 1/(1-d0);                     
    ls_failed = 1;                                    % this line search failed
  end

end
