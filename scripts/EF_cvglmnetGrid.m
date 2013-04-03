function CVerr = cvglmnet(x,y,nfolds,foldid,type,family,options,verbous)
% Do crossvalidation of glmnet model. The coordinate descent algorithm 
% chooses a set of lambda to optimize over based on the input given in
% the options struct.Parameter tuning should be done in such a way that for
% a fixed alpha, a set of lambda values are evaluated. Basically it does
% not matter if lambda corresponds across alpha values, as each cv-plot
% should inspected seperatly.
% So e.g. to find optimal tuning parameters, fit a total of 10 alpha
% values, beeing alpha = 0:0.1:1.lambdas will then be chosen according to 
% the specific alpha. 
% Call: CVerr = cvglmnet(x,y,nfolds,foldid,type,family,options,verbous)
% Example:
% x=randn(100,2000);
% y=randn(100,1);
% g2=randsample(2,100,true);
% CVerr=cvglmnet(x,y,100,[],'response','gaussian',glmnetSet,1);
% CVerr=cvglmnet(x,g2,100,[],'response','binomial',glmnetSet,1);
% x         : Covariates
% y         : Response (For now only elastic net = continous data supported
% nfolds    : How many folds to evaluate. nfolds = size(x,1) for LOOCV
% foldid    : Possibility for supplying own folding series. [] for nothing
% type      : Used in the glmnetPredict routine. (Now only "response" works)
% family    : Used in glmnet routine. (Now only "gaussian" and "binomial" work)
% options   : See function glmnetSet()
% verbous   : Print model plot
% 
% Written by Bjørn Skovlund Dissing (27-02-2010)
% Editted by Alan Gordon (30-05-2012)

glmnet_object = glmnet(x, y, family,options);
options.lambda = glmnet_object.lambda;
options.nlambda = length(options.lambda);
N = size(x,1);
if (isempty(foldid))
    foldid = randsample([repmat(1:nfolds,1,floor(N/nfolds)) 1:mod(N,nfolds)],N);
else
    nfolds = max(foldid)
end
%predmat = glmnetPredict(glmnet_object, type,x, options.lambda);
predmat = nan(N,length(options.lambda), length(options.alpha_set));
for j=1:length(options.alpha_set);
    options.alpha = options.alpha_set(j);
    for i=1:nfolds
        which=foldid==i;
        if verbous, disp(['Fitting fold # ' num2str(i) ' of ' num2str(nfolds)]);end
        cvfit = glmnet(x(~which,:), y(~which),family, options);
        predmat(which,:,j) = glmnetPredict(cvfit, type,x(which,:),options.lambda);
    end
end
%%yy=repmat(y,length(options.lambda));
yy=repmat(y,[1,length(options.lambda),length(options.alpha_set)]);

if strcmp(family,'gaussian')
    cvraw=(yy-predmat).^2;
elseif strcmp(family,'binomial')
    if     strcmp(type,'response')
        cvraw=-2*((yy==2).*log(predmat)+(yy==1).*log(1-predmat));
    elseif strcmp(type,'class')
        cvraw=double(yy~=predmat);
    end
elseif strcmp(family,'multinomial')
    error('Not implemented yet')
end
CVerr.cvm=squeeze(mean(cvraw));
CVerr.stderr=sqrt(squeeze(var(cvraw)/N));
CVerr.cvlo=CVerr.cvm-CVerr.stderr;
CVerr.cvup=CVerr.cvm+CVerr.stderr;

minError = min(min(CVerr.cvm));
[idxMinLambda,idxMinAlpha] = find(CVerr.cvm==(minError));

CVerr.lambda_min = options.lambda(idxMinLambda(1));
CVerr.alpha_min = options.alpha_set(idxMinAlpha(1));

% if there are several minima, choose largest lambda of the smallest cvm
%CVerr.lambda_min=max(options.lambda(CVerr.cvm<=min(CVerr.cvm)));
%Find stderr for lambda(min(sterr))
semin=CVerr.cvup(idxMinLambda(1),idxMinAlpha(1));
% find largest lambda which has a smaller mse than the stderr belonging to
% the largest of the lambda belonging to the smallest mse
% In other words, this defines the uncertainty of the min-cv, and the min
% cv-err could in essence be any model in this interval.
CVerrThisAlpha=CVerr.cvm(:,idxMinAlpha(1));
CVerr.lambda_1se=max(options.lambda(CVerrThisAlpha<semin));
CVerr.glmnetOptions=options;
CVerr.glmnet_object = glmnet_object;
if verbous, cvglmnetPlot(CVerr);end
end