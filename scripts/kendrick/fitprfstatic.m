function [params,paramsse,r,modelfit,numiters,responseB,rTRAIN,modelfitTRAIN,responseTRAIN] = ...
  fitprfstatic(stimulus,response,model,basismode,maxiter,wantresample,xvaltype,tol,outputfcn,metric,wantix,extraopt,wanttrain)

% function [params,paramsse,r,modelfit,numiters,responseB,rTRAIN,modelfitTRAIN,responseTRAIN] = ...
%   fitprfstatic(stimulus,response,model,basismode,maxiter,wantresample,xvaltype,tol,outputfcn,metric,wantix,extraopt,wanttrain)
%
% <stimulus> is the stimulus with dimensions A x B x C where A indicates different data points,
%   B indicates the components of a stimulus (e.g. pixels), and C indicates different stimulus frames
%   whose responses are to be averaged together.  note that C can equal 1.  note also that individual stimulus
%   frames (i.e. matrix slices with dimensions 1 x B x 1) can consist of all NaNs; we handle these cases
%   intelligently (i.e. we ignore these frames), giving us the ability to associate different data points 
%   with different numbers of stimulus frames.  <stimulus> can also be a function that returns the stimulus
%   when evaluated.  <stimulus> can also be 0, which means to use eye(A) where A is the size of <response>
%   along the first dimension.
% <response> is the data with dimensions A x D where A indicates different data points and D indicates
%   different trials.  note that D can equal 1.  all values should be finite and double-format.
% <model> is
%   {X Y Z W} where
%     X is the initial seed (1 x p)
%     Y are the bounds (2 x p).  NaNs in the first row indicate parameters to fix (and not change).
%     Z is a function handle that accepts two arguments, a set of
%       parameters (1 x p) and stimuli (N x B), and outputs
%       predicted responses (N x 1).  if <wantix>, then the stimuli
%       passed to the function handle is actually a cell vector {F G} where 
%       F is like before (N x B) and G is a vector of indices 
%       referring to the first dimension of <stimulus> (N x 1).
%     W (optional) is a function handle that transforms <stimulus> (N x B)
%       (possibly modifying the B dimension) into a new form, prior to model evaluation.
%       if you don't want to use this, omit it or pass in [].
%   in this case, we use lsqcurvefit.m optimization to fit the model.
%  OR
%   {M1 M2 M3 ...} where M1 is of the form {X Y Z W} as described above,
%   and the remaining Mi are of the form {F G H I} where
%     F is a function that takes fitted parameters (1 x p) from the previous model
%       and outputs an initial seed (1 x pnew) for the current model
%     G are the bounds (2 x pnew).  NaNs in the first row indicate parameters to fix (and not change).
%     H is a function that takes fitted parameters (1 x p) from the previous model
%       and outputs a function handle that accepts two arguments, a set of
%       parameters (1 x pnew) and stimuli (N x B), and outputs
%       predicted responses (N x 1).  if <wantix>, then the stimuli
%       passed to the function handle is actually a cell vector {F G} where
%       F is like before (N x B) and G is a vector of indices
%       referring to the first dimension of <stimulus> (N x 1).
%     I (optional) is a function that takes fitted parameters (1 x p) from the previous model
%       and outputs a function handle that transforms <stimulus> (N x B) (possibly modifying
%       the B dimension) into a new form, prior to model evaluation.  if you don't want to use
%       this, omit it or pass in [].
%   the numbers of parameters for the different models do not have to be the same.
%   also, there must be at least one model (i >= 1).
%  OR
%   a matrix of dimensions B x F with different basis functions
%   along the second dimension.  F must be >= 1.
%  OR
%   0 is a special case that means act as if we have delta basis functions,
%   e.g., eye(B).  the point of this is to enable us to speed-up the code.
% <basismode> (optional) is
%   0 means use ordinary least-squares estimation
%   {J K L} means use regularized regression as implemented in 
%     gradientdescent_wrapper.m.  J, K, and L should be the <fittype>, 
%     <randomtype>, and <stepsize> inputs to gradientdescent_wrapper.m
%     (we use defaults for the other inputs, except for <nodc> which is always
%     set to 1). 
%   the <basismode> input matters only when <model> is not the cell vector
%   case.  (in the <model> cell vector case, we always use lsqcurvefit.m 
%   optimization.)  default: {0 [3 0 0] []}.
% <maxiter> (optional) matters only when <model> is the cell vector case, 
%   and is the maximum number of iterations.  can be [M1 M2 M3 ...] where
%   Mi is the number to use for each successive model.  default: Inf.
% <wantresample> (optional) is
%   0 means do not bootstrap nor cross-validate (and just take the mean across trials).
%   j means treat each trial separately.
%   N or {N SEED} means number of bootstraps to take and the random number seed to use.
%     if the seed is [] or not provided, we use sum(100*clock).
%   [A B C ...] specifies the cross-validation scheme.  this should be a vector
%     of length A or D, depending on whether you want to cross-validate data points or 
%     trials (see <xvaltype>), and should have at least two elements.  elements of 
%     the vector should be non-negative integers and should cover either the 
%     range 0 to J for some positive J or the range 1 to J for some positive J.  
%     on the first cross-validation run, we train on all data not labeled 1 and 
%     then validate on data labeled 1.  on the second cross-validation run, we 
%     train on all data not labeled 2 and then validate on data labeled 2.  and so forth.
%     special case is -1 which means to use [1:S] where S is the total number of data points
%     or trials.
%  -M specifies the cross-validation scheme.  (this scheme subsumes the previous scheme
%     and provides additional functionality.)  M should be a matrix with dimensions
%     A x B where A >= 1 indicates different cross-validation cases and B >= 2 matches
%     the number of data points or trials, depending on whether you want to cross-validate
%     data points or trials (see <xvaltype>).  elements of M should be -1, 0, or 1.  
%     1 indicates items to use for training; -1 indicates items to use for testing.  
%     every cross-validation case must involve at least one training item and one
%     testing item.
%   default: 0.
% <xvaltype> (optional) matters only when <wantresample> is not 0 nor j, and should be 
%   either 0, which means to bootstrap or cross-validate data points, or 1, which
%   means to bootstrap or cross-validate trials.  default: 0.
% <tol> (optional) is the tolerance to use for TolFun and TolX.  default: 1e-6.
%   can be [A B] where A is the tolerance on TolFun and B is the tolerance on TolX.
%   can also be {T1 T2 T3 ...} where Ti is the tolerance to use for each successive model.
%   matters only when <model> is the cell vector case.
% <outputfcn> (optional) is a function suitable for use as an 'OutputFcn'.
%   <outputfcn> matters only when <model> is the cell vector case.  default: @alwayszero.
%   special case is {FN} where FN is a function that is like an 'OutputFcn' but is 
%   expecting two additional arguments:
%     data, a column vector with the data being fitted
%     datase, a column vector with the standard error on the data
%   note that data and datase will not reflect the original data but rather will reflect
%   division by the std dev of the original data.
% <metric> (optional) is a function that accepts two column vectors (the first is the
%   model, the second is the data) and outputs a value.  default: @(x,y) calccod(x,y,[],[],0).
% <wantix> (optional) matters only when <model> is the cell vector case and is whether to 
%   pass not only stimuli but also indices to the model (see <model>).  the point of <wantix>
%   is to allow models that have parameters whose influence is dependent on the specific
%   index of the data point being predicted.  default: 0.
% <extraopt> (optional) is a cell vector of extra parameter/value pairs to use
%   in the lsqcurvefit.m optimization options (e.g. {'Display' 'final'}).  default: {}.
% <wanttrain> (optional) is whether to set the <rTRAIN>, <modelfitTRAIN>, and <responseTRAIN>
%   output variables.  this input matters only when <wantresample> is a cross-validation case.
%   default: 0.
%
% estimate static PRF model.
%
% return:
%  <params> is 1 x p with the parameters.
%    when <wantresample>, we return <params> as b x p where the rows correspond to different resamples.
%  <paramsse> is 1 x p with the standard deviation across resamples.
%    NaNs are returned if isequal(<wantresample>,0).
%  <r> is the <metric> between <modelfit> and <responseB>.
%    when <wantresample> is 0 or j or N or {N SEED}, this is straightforward:
%     <modelfit> is A x 1 with the model fit.
%       if <wantresample> is j or N or {N SEED}, we calculate the mean model fit across cases.
%     <responseB> is A x 1 with the data averaged across trials.
%    when <wantresample> is [A B C ...] or -M, this is a bit tricky:
%     <modelfit> is V x 1 with the cross-validated model predictions.
%     <responseB> is V x 1 with the cross-validation data.
%  <modelfit> is as described above.
%  <numiters> is 1 x 1 with the number of iterations taken.  when <wantresample>, we return <numiters>
%    as 1 x b where each element corresponds to a different resampling.  in the OLS case, <numiters>
%    is returned as [].
%  <responseB> is as described above.
%  <rTRAIN> is the <metric> between <modelfitTRAIN> and <responseTRAIN>.  these variables get 
%    computed only when <wantresample> is [A B C ...] or -M and when <wanttrain>.
%    <modelfitTRAIN> is T x 1 with the model fits to the training data.
%    <responseTRAIN> is T x 1 with the training data.
%
% history:
% 2011/11/26 - add additional case for outputfcn
% 2011/10/09 - outputfcn no longer has all those inputs.  now it is a usual OutputFcn.  also, we automatically
%              enforce a sanity check via outputfcnsanitycheck.m.
% 2011/09/22 - add <wanttrain> and the corresponding output variables
% 2011/07/02 - oops, fix minor bug (would have crashed)
% 2011/06/27 - previously, we would crash when the data had no variance (e.g. all zeros).  now we don't crash.
% 2011/06/09 - add <extraopt> input
% 2010/01/12 - make <maxiter> more general
% 2010/12/06 - introduce nanreplace into the lsqcurvefit call to prevent crazy parametric model return values from causing chaos.
% 2010/12/02 - implement {N SEED} for <wantresample>
% 2010/11/13 - change name of <gdparams> to <basismode>; implement <basismode> 0 case; implement <stimulus> 0 case
% 2010/11/02 - implement -1 case for <wantresample>
% 2010/09/28 - fix bug (it would have crashed)
% 2010/09/27 - implement j case for <wantresample>
% 2010/09/09 - add -M case for <wantresample>
% 2010/08/25 - too many changes to document!  let's just start again fresh.
% 2010/08/16 - implemenet NaN case for the bounds (to fix parameters)
% 2010/08/15 - allow cell vector case for <tol>
% 2010/08/14 - add funvalcheck; change model to require @(ss) before 3 and 4th entries in <model>; 
% 2010/08/13 - add reporting of seeds in the multiple model case; add reporting of paramteers estimated
% 2010/08/12 - allow NaN cases in <stimulus>; implement the W parameter in <model>; add responseB output
% 2010/08/11 - change default for <metric> and <extraopt>.
% 2010/08/10 - first version

% TODO: perhaps we should save the final model function handle?
% TODO: also get the bootstrap "GROUP" input from fitprf.m?

% input
if ~exist('basismode','var') || isempty(basismode)
  basismode = {0 [3 0 0] []};
end
if ~exist('maxiter','var') || isempty(maxiter)
  maxiter = Inf;
end
if ~exist('wantresample','var') || isempty(wantresample)
  wantresample = 0;
end
if ~exist('xvaltype','var') || isempty(xvaltype)
  xvaltype = 0;
end
if ~exist('tol','var') || isempty(tol)
  tol = 1e-6;
end
if ~exist('outputfcn','var') || isempty(outputfcn)
  outputfcn = @alwayszero;
end
if ~exist('metric','var') || isempty(metric)
  metric = @(x,y) calccod(x,y,[],[],0);
end
if ~exist('wantix','var') || isempty(wantix)
  wantix = 0;
end
if ~exist('extraopt','var') || isempty(extraopt)
  extraopt = {};
end
if ~exist('wanttrain','var') || isempty(wanttrain)
  wanttrain = 0;
end
if iscell(model)
  if ~iscell(model{1})
    model = {model};
  end
  for zz=1:length(model)
    if length(model{zz}) < 4 || isempty(model{zz}{4})
      if zz==1
        model{zz}{4} = @identity;
      else
        model{zz}{4} = @(ss) @identity;
      end
    end
  end
end
if length(maxiter) ~= length(model)
  assert(length(maxiter)==1);
  maxiter = repmat(maxiter,[1 length(model)]);
end
if ~iscell(tol)
  tol = {tol};
end
for zz=1:length(tol)
  if length(tol{zz}) == 1
    tol{zz} = repmat(tol{zz},[1 2]);
  end
end
if length(tol) ~= length(model)
  assert(length(tol)==1);
  tol = repmat(tol,[1 length(model)]);
end
if iscell(outputfcn)
  outputfcn = outputfcn{1};
  outputfcnextra = 1;
else
  outputfcnextra = 0;
end
  
% calc stuff
if isequal(stimulus,0)
  stimulus = eye(size(response,1));
end
if isa(stimulus,'function_handle')
  stimulus = feval(stimulus);
end
numdata = size(stimulus,1);
numtrials = size(response,2);
numxval = choose(xvaltype==0,numdata,numtrials);

% calc stuff for NaN handling
valid = permute(~isnan(stimulus(:,1,:)),[3 1 2]);  % logical matrix C x A indicating which frames are valid

% convert wantresample -1 case
if isequal(wantresample,-1)
  wantresample = 1:numxval;
end

% convert wantresample case [A B C ...] to -M format
if ~iscell(wantresample) && numel(wantresample) > 1 && all(wantresample(:) >= 0)
  wantresampleNEW = [];
  for p=1:max(wantresample)
    wantresampleNEW(p,:) = zeros(1,numxval);
    wantresampleNEW(p,wantresample~=p) = 1;  
    wantresampleNEW(p,wantresample==p) = -1;
  end
  wantresample = -wantresampleNEW; clear wantresampleNEW;
end

% handle fitting in the cross-validation cases
isxval = ~iscell(wantresample) && numel(wantresample) > 1;
if isxval
  
  % do the fitting
  params = []; numiters = []; responseB = []; modelfit = []; responseTRAIN = []; modelfitTRAIN = [];
  for cnt=1:size(wantresample,1)
  
    % figure out train and test indices
    test = find(-wantresample(cnt,:) < 0); assert(~isempty(test));
    train = find(-wantresample(cnt,:) > 0); assert(~isempty(train));
    
    % precompute
    if xvaltype==0
      data = mean(response(train,:),2);
      datase = std(response(train,:),[],2)/sqrt(numtrials);
      trainx = train;
      testx = test;
    else
      data = mean(response(:,train),2);
      datase = std(response(:,train),[],2)/sqrt(length(train));
      trainx = 1:numdata;
      testx = 1:numdata;
    end
    datastd = choose(std(data)==0,1,std(data));

    % fit the training data
    if ~iscell(model)
      if model==0
        X = nanmean(stimulus(trainx,:,:),3);
        Xt = nanmean(stimulus(testx,:,:),3);
      else
        X = nanmean(stimulus(trainx,:,:),3)*model;
        Xt = nanmean(stimulus(testx,:,:),3)*model;
      end
      if iscell(basismode)
        [h,dc,numiters0] = gradientdescent_wrapper(data,X,basismode{1},[],basismode{2},basismode{3},[],[],[],[],[],[],[],[],[],1);
        params(cnt,:) = full(h');
        numiters(cnt) = numiters0;
      else
        params0 = olsmatrix(X) * data;
        params(cnt,:) = params0';
      end
      fprintf('estimated params is %s.\n',mat2str(params(cnt,:),5));
    else
      for zz=1:length(model)
        if zz==1
          seed = model{zz}{1};
          mmm = model{zz}{3};
          tran = model{zz}{4};
        else
          seed = feval(model{zz}{1},params0);
          mmm = feval(model{zz}{3},params0);
          tran = feval(model{zz}{4},params0);
        end
        lb = model{zz}{2}(1,:); assert(length(seed)==length(lb));
        ub = model{zz}{2}(2,:);
        ix = ~isnan(lb);  % indices of free parameters
        fprintf('for model %d of %d, the seed is ',zz,length(model)); fprintf('%.3f ',seed); fprintf('\n');
        stim = subscript(squish(permute(feval(tran,stimulus(trainx,:,:)),[3 1 2]),2),{vflatten(valid(:,trainx)) ':'});
        stimindices = subscript(squish(permute(repmat(trainx(:),[1 1 size(stimulus,3)]),[3 1 2]),2),{vflatten(valid(:,trainx)) ':'});
        fun = @(pp) chunkfun(feval(mmm,copymatrix(seed,ix,pp),choose(wantix,{stim stimindices},stim)),sum(valid(:,trainx),1),@(y)nanmean(y,1))' / datastd;
        outputfcnextras = choose(outputfcnextra,{data/datastd datase/datastd},{});
        options = optimset('Display','iter','FunValCheck','on','MaxFunEvals',Inf,'MaxIter',maxiter(zz), ...
                           'TolFun',tol{zz}(1),'TolX',tol{zz}(2), ...
                           'OutputFcn',@(a,b,c) outputfcnsanitycheck(a,b,c,tol{zz}(1),10) | feval(outputfcn,a,b,c,outputfcnextras{:}), ...
                           extraopt{:});
        [params0,d,d,exitflag,output] = lsqcurvefit(@(x,y)nanreplace(fun(x),0,2),seed(ix),[],data/datastd,lb(ix),ub(ix),options);
        assert(exitflag >= -1);
        params0 = copymatrix(seed,ix,params0);
        fprintf('estimated params is %s.\n',mat2str(params0,5));
      end
      params(cnt,:) = params0;
      numiters(cnt) = output.iterations;
    end

    % construct the data and model fit for the testing data
    if xvaltype==0
      responseB = cat(1,responseB,mean(response(test,:),2));
    else
      responseB = cat(1,responseB,mean(response(:,test),2));
    end
    if ~iscell(model)
      modelfit = cat(1,modelfit,Xt*params(cnt,:)');
    else
      stim = subscript(squish(permute(feval(tran,stimulus(testx,:,:)),[3 1 2]),2),{vflatten(valid(:,testx)) ':'});
      stimindices = subscript(squish(permute(repmat(testx(:),[1 1 size(stimulus,3)]),[3 1 2]),2),{vflatten(valid(:,testx)) ':'});
      modelfit = cat(1,modelfit,chunkfun(feval(mmm,params(cnt,:),choose(wantix,{stim stimindices},stim)),sum(valid(:,testx),1),@(y)nanmean(y,1))');
    end
    
    % do the same for the training data (if desired)
    if wanttrain
      if xvaltype==0
        responseTRAIN = cat(1,responseTRAIN,mean(response(train,:),2));
      else
        responseTRAIN = cat(1,responseTRAIN,mean(response(:,train),2));
      end
      if ~iscell(model)
        modelfitTRAIN = cat(1,modelfitTRAIN,X*params(cnt,:)');
      else
        stim = subscript(squish(permute(feval(tran,stimulus(trainx,:,:)),[3 1 2]),2),{vflatten(valid(:,trainx)) ':'});
        stimindices = subscript(squish(permute(repmat(trainx(:),[1 1 size(stimulus,3)]),[3 1 2]),2),{vflatten(valid(:,trainx)) ':'});
        modelfitTRAIN = cat(1,modelfitTRAIN,chunkfun(feval(mmm,params(cnt,:),choose(wantix,{stim stimindices},stim)),sum(valid(:,trainx),1),@(y)nanmean(y,1))');
      end
    end

  end

% handle fitting in the other cases
else

  % figure out bootfun (for data) and bootfun2 (for stimulus)
  if ~isequal(wantresample,0)
    bootfun = {}; bootfun2 = {};
    if isequal(wantresample,j)
      for xx=1:numtrials
        bootfun{xx} = @(x) x(:,xx,:);
        bootfun2{xx} = @identity;
        trainx = 1:numdata;
      end
    else
      if iscell(wantresample)
        wantresample0 = wantresample{1};
        if length(wantresample) < 2 || isempty(wantresample{2})
          setrandstate;
        else
          setrandstate(wantresample{2});
        end
      else
        wantresample0 = wantresample;
        setrandstate;
      end
      if xvaltype==0
        for xx=1:wantresample0
          ix = ceil(rand(1,numdata)*numdata);
          bootfun{xx} = @(x) x(ix,:,:);
          bootfun2{xx} = @(x) x(ix,:,:);
          trainx = ix;
        end
      else
        for xx=1:wantresample0
          ix = ceil(rand(1,numtrials)*numtrials);
          bootfun{xx} = @(x) x(:,ix,:);
          bootfun2{xx} = @identity;
          trainx = 1:numdata;
        end
      end
    end
  else
    bootfun = {@identity};  % this does nothing
    bootfun2 = {@identity};  % this does nothing
    trainx = 1:numdata;
  end

  % precompute
  data = response; datase = std(response,[],2)/sqrt(numtrials); datastd = choose(std(data(:))==0,1,std(data(:)));
  if ~iscell(model)
    if model==0
      X = nanmean(stimulus,3);
    else
      X = nanmean(stimulus,3)*model;
    end
  end

  % do the fitting
  params = []; numiters = [];
  for cnt=1:length(bootfun)
    
    % do it
    if ~iscell(model)
      if iscell(basismode)
        [h,dc,numiters0] = gradientdescent_wrapper(mean(feval(bootfun{cnt},data),2),feval(bootfun2{cnt},X),basismode{1},[],basismode{2},basismode{3},[],[],[],[],[],[],[],[],[],1);
        params(cnt,:) = full(h');
        numiters(cnt) = numiters0;
      else
        params0 = olsmatrix(feval(bootfun2{cnt},X)) * mean(feval(bootfun{cnt},data),2);
        params(cnt,:) = params0';
      end
      fprintf('estimated params is %s.\n',mat2str(params(cnt,:),5));
    else
      for zz=1:length(model)
        if zz==1
          seed = model{zz}{1};
          mmm = model{zz}{3};
          tran = model{zz}{4};
        else
          seed = feval(model{zz}{1},params0);
          mmm = feval(model{zz}{3},params0);
          tran = feval(model{zz}{4},params0);
        end
        lb = model{zz}{2}(1,:); assert(length(seed)==length(lb));
        ub = model{zz}{2}(2,:);
        ix = ~isnan(lb);  % indices of free parameters
        fprintf('for model %d of %d, the seed is ',zz,length(model)); fprintf('%.3f ',seed); fprintf('\n');
        stim = subscript(squish(permute(feval(tran,feval(bootfun2{cnt},stimulus)),[3 1 2]),2),{vflatten(feval(bootfun2{cnt},valid')') ':'});
        stimindices = subscript(squish(permute(repmat(trainx(:),[1 1 size(stimulus,3)]),[3 1 2]),2),{vflatten(feval(bootfun2{cnt},valid')') ':'});
        fun = @(pp) chunkfun(feval(mmm,copymatrix(seed,ix,pp),choose(wantix,{stim stimindices},stim)),sum(feval(bootfun2{cnt},valid')',1),@(y)nanmean(y,1))' / datastd;
        outputfcnextras = choose(outputfcnextra,{mean(feval(bootfun{cnt},data),2)/datastd feval(bootfun2{cnt},datase)/datastd},{});
        options = optimset('Display','iter','FunValCheck','on','MaxFunEvals',Inf,'MaxIter',maxiter(zz), ...
                           'TolFun',tol{zz}(1),'TolX',tol{zz}(2), ...
                           'OutputFcn',@(a,b,c) outputfcnsanitycheck(a,b,c,tol{zz}(1),10) | feval(outputfcn,a,b,c,outputfcnextras{:}), ...
                           extraopt{:});
        [params0,d,d,exitflag,output] = lsqcurvefit(@(x,y)nanreplace(fun(x),0,2),seed(ix),[],mean(feval(bootfun{cnt},data),2)/datastd,lb(ix),ub(ix),options);
        assert(exitflag >= -1);
        params0 = copymatrix(seed,ix,params0);
        fprintf('estimated params is %s.\n',mat2str(params0,5));
      end
      params(cnt,:) = params0;
      numiters(cnt) = output.iterations;
    end

  end

  % construct the data and model fit
  responseB = mean(response,2);
  if ~iscell(model)
    modelfit = mean(rowfun(params,@(x) X*x',2),2);
  else
    stim = subscript(squish(permute(feval(tran,stimulus),[3 1 2]),2),{vflatten(valid) ':'});
    stimindices = subscript(squish(permute(repmat((1:numdata)',[1 1 size(stimulus,3)]),[3 1 2]),2),{vflatten(valid) ':'});
    modelfit = mean(rowfun(params,@(x) chunkfun(feval(mmm,x,choose(wantix,{stim stimindices},stim)),sum(valid,1),@(y)nanmean(y,1)),1),1)';
  end
  
  % these are always []
  responseTRAIN = []; modelfitTRAIN = [];

end

% how well did we do?
r = feval(metric,modelfit,responseB);
if isxval && wanttrain
  rTRAIN = feval(metric,modelfitTRAIN,responseTRAIN);
else
  rTRAIN = [];
end

% deal with standard error
if ~isequal(wantresample,0)
  paramsse = std(params,[],1);
else
  paramsse = NaN*zeros(1,size(params,2));
end








% accepts:
% %     params, optimValues, state, data, datase, trainx, seed, ix
% %         trainx is a vector of indices referring to the original data, indicating which
% %           data points are currently being fitted
% %         seed is a vector with the initial seed
% %         ix is a logical vector indicating the parameters to change
