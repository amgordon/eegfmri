
function eegClassifyROI(subjects, dur, shufl, contrasts)

% this function sets up the necessary variables to perform binary
% classification on eeg data
%
% it is setup to perform classification on ROI's
% inputs:
% expt:     experiment string, either 'PC1' for posner or 'SS2' for memory
% subjects: a cell array of strings containing the subject numbers
% dur:      2 element vector that delimits the start and end of the
%           classification window
% shulf:    flag specifying that the class labes

addpath(genpath('/biac4/wagner/biac3/wagner7/ecog/scripts/alex/lib'))

for s = subjects
	S=[];
		
	S.subj = s{1};
	S.binsize = 50; % bin size in ms
	S.sldwin = 50; % sliding window
	S.rois = channel_rois; % channels by roi
	S.classify_win_def = dur; % classification window default
	S.shufl = shufl; % shuffle labels flag
	S.nfolds = []; % number of cross validation folds
	S.nbands = 7; % number of bands
	S.trial_reject_method = 'var'; % method of trial rejection {'var', 'corr'}
	S.trial_rejection_threh = 1.5; % trial rejection threshold
	S.trial_rejection_channels = 'LPS'; % use these channels to detect bad trials
	toolboxes = {'liblinear', 'libsvm','glmnet'}; % classification toolboxes
	S.toolbox = toolboxes{3};
	S.perfmetric = 'bac'; % optimize for balanced accuracy
	
	% bands
	%S.classificationbands = {'amp','delta','theta','alpha','beta','gamma','hgamma'};
	% contrasts
	S.contrasts = contrasts;
	
	% options for each toolbox
	if strcmp(S.toolbox, 'liblinear')
		addpath(genpath('/biac4/wagner/biac3/wagner7/ecog/scripts/alex/classification_scripts/liblinear-1.8/matlab/'))
		S.CparamSet = logspace(-2,3,10);
		S.solver_type = '1'; % solver for liblinear
		S.trainOpts_orig = ['-q -s ' S.solver_type];
		S.predictOpts = [];%'-b 1';
		nit = length(S.CparamSet);
		
	elseif strcmp(S.toolbox, 'libsvm')
		addpath(genpath('/biac4/wagner/biac3/wagner7/ecog/scripts/alex/classification_scripts/libsvm-3.11/matlab/'))
		
		S.GammaSet = logspace(-3,-1,3);%[1/105 1/200];
		S.CparamSet = logspace(-3,1,5);
		nit = length(S.CparamSet);
		
		S.solver_type = '0';
		S.trainOpts_orig = ['-q -s ' S.solver_type '-t 2'];
		S.predictOpts = [];
	elseif strcmp(S.toolbox, 'glmnet')
		addpath(genpath('/biac4/wagner/biac3/wagner7/ecog/scripts/alex/classification_scripts/glmnet_matlab'))
		S.opts = glmnetSet; % glmnet
		S.opts.alpha = 0.9; % L2-L1 parameter
		S.opts.nlambda = 150; % number of lambda
		S.opts.maxit = 150; % max number of iterations
		S.opts.perfmetric = S.perfmetric;
		S.solver_type = ['alpha' strrep(num2str(S.opts.alpha),'.','p')];
		S.lambdafolds = 5;  % lambda parameter selection folds.
	end
	
	% load the data and set up directories
	S.mainpath = '/biac4/wagner/biac3/wagner5/alan/eegfmri/fmri_data/'; % returns experimental directory
	datapath = [ S.mainpath S.subj '/erpData/'];
	temp = load([datapath 'trial_data_downsamp_LP30Hz_zscored.mat']);
	S.data = temp.S; temp =[];
	
	savepath = [datapath 'classificationResults/'];
	mkdir(savepath)
	
	tc1 = strrep(num2str((S.classify_win_def(1))),'.','p');
	tc2 = strrep(num2str(S.classify_win_def(2)),'.','p');
	
	% load behavioral data
	load([datapath '/results/behdata.mat']);
	S.behavioral = behdata;
	S.binstoclassify = ((S.data.bins(:,1) >= 1000*S.classify_win_def(1)) ...
		& (S.data.bins(:,2) <= 1000*S.classify_win_def(2)));
	
	
	% zscore and find bad trials based on amplitude data
	channels = S.rois.(S.trial_rejection_channels);
	thr = S.trial_rejection_threh;
	X = S.data.bin_trials(channels,:,:);
	[nch ntr nbins]=size(X);
	
	% zscore the channels
	XZ = zscore(X(:,:),[],2);
	XZ = reshape(XZ,[nch ntr nbins]);
	
	% collapse across channels
	X = permute(XZ,[2 1 3]);
	
	if strcmp(S.trial_reject_method,'var')
		% calculate the trial covariane matrix
		trial_var = var(X(:,:),0,2);
		
		% calculate an upper bound on the variance for each trial
		high_var_bound = mean(trial_var)+thr*std(trial_var);
		
		% deterime good trials
		S.good_trials = (trial_var < high_var_bound);
		S.good_trials = S.good_trials(:);
		
	elseif strcmp(S.trial_reject_method,'corr')
		% mean trial
		mt = squeeze(mean(X(:,:),1));
		
		% correlation of each trial to the mean trial
		mtc = corr(mt',squeeze(X(:,:))');
		
		lower_corr_bound = mean(mtc)-thr*std(mtc);
		
		% good trials
		S.good_trials = (mtc > lower_corr_bound);
		S.good_trials = S.good_trials(:);
		
	end
	
	XZ =[]; X = [];
	
	for cont = 1:size(S.contrasts,1);
		S.conds = S.contrasts(cont,:); % classification conditions
        nevents = numel(S.behavioral.(S.conds{1}));
		S.cond_trials = [S.behavioral.(S.conds{1})&S.good_trials(1:nevents) , ...
			S.behavioral.(S.conds{2})&S.good_trials(1:nevents)];
		S.indicator_trials = logical(sum(S.cond_trials,2));
		S.Y = [S.indicator_trials(S.cond_trials(:,1)==1); ...
			-1*S.indicator_trials(S.cond_trials(:,2)==1)];
		S.Yorig = S.Y;
		
		for ro = 1:numel(S.rois.roinames)
            
            
			param_model=[];
			S.whichroi = S.rois.roinames{ro};
            display(['Classifying ROI ' S.whichroi])
			S.channels = S.rois.(S.whichroi); % classify across channels
			S.nchan = numel(S.channels);
			
			filename = [savepath S.subj '_' ...
				S.conds{1} '_' S.conds{2} '_' S.toolbox 'solver' S.solver_type ...
				'_time' tc1 '_' tc2 '_' S.whichroi '_binsize' num2str(S.binsize)];
				
			for band = [1]
				
				S.band_to_classify = band;
				% number of predictors
				S.nPredictors = sum(S.binstoclassify)*numel(S.channels)*numel(band);
				
				% for each band calculate the mean sd of each channel
				X = S.data.bin_trials(S.channels,:,:);
				[nch ntr nbins]=size(X);
				mX = mean(X(:,:),2);
				sdX = std(X(:,:),[],2);
				
				% standarize the data
				X = bsxfun(@minus,X(:,:),mX);
				X = bsxfun(@rdivide,X,sdX);
				X = reshape(X,[nch ntr nbins]);
				
				% data for condition 1
				c1 = S.cond_trials(:,1)==1;
				S.nc1 = sum(c1); % number of trials for condition 1
				Xc1 = X(:,c1 ,S.binstoclassify);
				Xc1 = permute(Xc1,[2 1 3]);
				Xc1 = reshape(Xc1,S.nc1,S.nPredictors);
				% data for condition 2
				c2	= S.cond_trials(:,2)==1;
				S.nc2 = sum(c2);
				Xc2 = X(:,c2,S.binstoclassify);
				Xc2 = permute(Xc2,[2 1 3]);
				Xc2 = reshape(Xc2,S.nc2,S.nPredictors);
				
				% final predictor matrix
				S.X = [Xc1;Xc2]; X =[]; Xc1 =[]; Xc2 =[];
				
				% total number of observations
				S.N = size(S.X,1);
				S.nfolds = 10; % number of cross validation folds
				
				% actual classification....
				if strcmp(S.toolbox, 'liblinear')
					Sout = [];
					for c = 1:nit
						% select a c
						if S.shufl
							
							S.Y = shuffle(S.Yorig);
							shuflstr = 'ShuffledLabels';
							temp = load([filename1 '.mat']);
							S.Cparam = temp.param_model.optC_bac(band);
							temp=[];
						else
							shuflstr = '';
							S.Cparam = S.CparamSet(c);
						end
						
						fprintf('\nClassifying with C = %f\n',S.Cparam)
						% set training parameters
						S.trainOpts = [S.trainOpts_orig ' -c ' num2str(S.Cparam)];
						% classify
						Sout = ecogSVMClassifyROI2(S);
						param_model.Cparams(c) = S.Cparam;
						% performance by rois for a given c
						param_model.Fscores(band,c,:) = Sout.perf.F_Score;
						param_model.Acc(band,c,:) = Sout.perf.Accuracy;
						param_model.dPrime(band,c,:) = Sout.perf.dPrime;
						param_model.acc_class1(band,c,:) = Sout.perf.acc_class1;
						param_model.acc_class2(band,c,:) = Sout.perf.acc_class2;
						param_model.bac(band,c,:) = Sout.perf.bac;
						param_model.mcc(band,c,:) = Sout.perf.mcc;
						
						param_model.perf(band,c,:) = Sout.perf;
						param_model.weights{band,c} = Sout.weights;
						param_model.predictions(band,c,:) =Sout.predictions;
						param_model.actual(band,c,:) =Sout.actual_outs;
						param_model.probs(band,c,:) =Sout.prob_est;
						
					end
					
				elseif strcmp(S.toolbox, 'libsvm')
					Sout = [];
					
					for c = 1:length(S.CparamSet)
						for g = 1:length(S.GammaSet)
							% select a c
							S.Cparam = S.CparamSet(c);
							% select gamma
							S.Gparam = S.GammaSet(g);
							fprintf('\nClassifying with C = %f, gamma = %f\n',S.Cparam,S.Gparam)
							% set training parameters
							S.trainOpts = [S.trainOpts_orig ' -c ' num2str(S.Cparam) ' -g ' num2str(S.Gparam)];
							% classify
							Sout = ecogSVMClassifyROI2(S,allchdata,Sout);
							param_model.Cparams(c) = S.Cparam;
							param_model.Gamma_params(g) = S.Gparam;
							% performance by roi for a given g and c
							param_model.Fscores(band,g,c,:) = Sout.perf.F_Score;
							param_model.Acc(band,g,c,:) = Sout.perf.Accuracy;
							param_model.dPrime(band,g,c,:) = Sout.perf.dPrime;
							param_model.acc_class1(band,g,c,:) = Sout.perf.acc_class1;
							param_model.acc_class2(band,g,c,:) = Sout.perf.acc_class2;
							param_model.bac(band,g,c,:) = Sout.perf.bac;
							param_model.mcc(band,g,c,:) = Sout.perf.mcc;
							
							param_model.perf(band,g,c,:) = Sout.perf;
							param_model.predictions(band,g,c,:) =Sout.predictions;
							param_model.actual(band,g,c,:) =Sout.actual_outs;
							
							param_model.perf(band,g,c) =Sout.perf;
							
						end
					end
					
				elseif strcmp(S.toolbox, 'glmnet')
					
					if S.shufl
						shuflstr = 'ShuffledLabels';
					else
						shuflstr = '';
					end
					
					it_dist = shuffle(1:length(S.Y));
					iters = ceil(it_dist/(max(it_dist)/S.nfolds));
					
					X = double(S.X); %data must in double format.
					
					% scale data from 0 to 1
					X = (X - repmat(min(X,[],1),size(X,1),1))* ...
						diag(1./(max(X,[],1)-min(X,[],1)));
					
					Y = S.Y;
					
					%initialize hypothesis vector
					prob_est_cell = cell(S.nfolds,1);
					prob_est = nan(S.N,1);
					betas = zeros(S.nfolds,S.nPredictors);
					a0 = zeros(S.nfolds,1);
					minlambdas = zeros(S.nfolds,1);
					lambdafolds = S.lambdafolds;
					opts = S.opts;
					
					parfor fold = 1:S.nfolds
						display(sprintf('fold # %d',fold))
						trainIdx = iters~=fold;
						testIdx = iters==fold;
						
						TrainX =  X(trainIdx,:);
						TrainY =  Y(trainIdx);
						
						TestX = X(testIdx,:);
						
						m = cvglmnet2(TrainX,(TrainY==1)+1,lambdafolds,[],'class','binomial',opts,0);
						minlambdas(fold)=m.lambda_min;
						betas(fold,:) = m.glmnet_object.beta(:,m.glmnet_object.lambda== m.lambda_min);
						a0(fold) = m.glmnet_object.a0(m.glmnet_object.lambda== m.lambda_min);
						prob_est_cell{fold} = glmnetPredict(m.glmnet_object, 'response', TestX, minlambdas(fold));
					end
					
					for fold =1:S.nfolds
						prob_est(iters==fold) = prob_est_cell{fold};
					end

					
					predictions = prob_est > 0.5;
					perf = class_perf_metrics(Y==1,predictions);
					param_model.maxperf_Fs(band)			= perf.F1score;
					param_model.maxperf_dP(band)			= perf.dP;
					param_model.maxperf_Ac(band)			= perf.acc;
					param_model.maxperf_BAC(band)			= perf.bac;
					param_model.maxperf_AUC(band)			= perf.AUC;
					param_model.maxperf_mcc(band)			= perf.mcc;
					param_model.betas(band,1:S.nfolds,:)	= betas;
					param_model.a0(band,1:S.nfolds)			= a0;
					param_model.minlambdas(band,1:S.nfolds)	= minlambdas;
					param_model.prob_est(band,:)			= prob_est;
					param_model.predictions(band,:)			= predictions;
					param_model.actual(band,:)				= Y;
				end
				
				%param_model.Sout = Sout;
				
				save_now=1;
				
				if save_now
					%savepath = [Sout.mainpath '/oldNewClassification/data/'
					%'subj' Sout.subjnum '/classify_roi/' Sout.expt '/'];
					save([filename shuflstr '.mat'],'param_model')
				end
				
			end
			
			% 			% summary of results
			% 			[param_model.maxperf_Fs param_model.max_Fs] = max(param_model.Fscores,[],2);
			% 			[param_model.maxperf_dP param_model.max_dP] = max(param_model.dPrime,[],2);
			% 			[param_model.maxperf_Ac param_model.max_Ac] = max(param_model.Acc,[],2);
			%
			% 			[param_model.maxperf_BAC i] = max(param_model.bac,[],2);
			% 			param_model.optC_bac = S.CparamSet(i);
			% 			[~,j] = max(param_model.maxperf_BAC);
			% 			[param_model.maxClassiferWeights{1:7}] = deal(param_model.weights{sub2ind([7 10],1:7,i')});
			%
			% 			[param_model.maxperf_mcc i] = max(param_model.mcc,[],2);
			% 			param_model.optC_mcc = S.CparamSet(i);
			
		end
	end
end


% liblinear options

% -s type : set type of solver (default 1)
% 	0 -- L2-regularized logistic regression (primal)
% 	1 -- L2-regularized L2-loss support vector classification (dual)
% 	2 -- L2-regularized L2-loss support vector classification (primal)
% 	3 -- L2-regularized L1-loss support vector classification (dual)
% 	4 -- multi-class support vector classification by Crammer and Singer
% 	5 -- L1-regularized L2-loss support vector classification
% 	6 -- L1-regularized logistic regression
% 	7 -- L2-regularized logistic regression (dual)
% -c cost : set the parameter C (default 1)
% -e epsilon : set tolerance of termination criterion
% 	-s 0 and 2
% 		|f'(w)|_2 <= eps*min(pos,neg)/l*|f'(w0)|_2,
% 		where f is the primal function and pos/neg are # of
% 		positive/negative data (default 0.01)
% 	-s 1, 3, 4 and 7
% 		Dual maximal violation <= eps; similar to libsvm (default 0.1)
% 	-s 5 and 6
% 		|f'(w)|_1 <= eps*min(pos,neg)/l*|f'(w0)|_1,
% 		where f is the primal function (default 0.01)
% -B bias : if bias >= 0, instance x becomes [x; bias]; if < 0, no bias term added (default -1)
% -wi weight: weights adjust the parameter C of different classes (see README for details)
% -v n: n-fold cross validation mode
% -q : quiet mode (no outputs)
% col:
% 	if 'col' is setted, training_instance_matrix is parsed in column
% 	format, otherwise is in row format
% libsvm options
%  "Usage: model = svmtrain(training_label_vector, training_instance_matrix, 'libsvm_options');\n"
% 	"libsvm_options:\n"
% 	"-s svm_type : set type of SVM (default 0)\n"
% 	"	0 -- C-SVC\n"
% 	"	1 -- nu-SVC\n"
% 	"	2 -- one-class SVM\n"
% 	"	3 -- epsilon-SVR\n"
% 	"	4 -- nu-SVR\n"
% 	"-t kernel_type : set type of kernel function (default 2)\n"
% 	"	0 -- linear: u'*v\n"
% 	"	1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
% 	"	2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
% 	"	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
% 	"	4 -- precomputed kernel (kernel values in training_instance_matrix)\n"
% 	"-d degree : set degree in kernel function (default 3)\n"
% 	"-g gamma : set gamma in kernel function (default 1/num_features)\n"
% 	"-r coef0 : set coef0 in kernel function (default 0)\n"
% 	"-c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"
% 	"-n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
% 	"-p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
% 	"-m cachesize : set cache memory size in MB (default 100)\n"
% 	"-e epsilon : set tolerance of termination criterion (default 0.001)\n"
% 	"-h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)\n"
% 	"-b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)\n"
% 	"-wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)\n"
% 	"-v n : n-fold cross validation mode\n"
% 	"-q : quiet mode (no outputs)\n"
