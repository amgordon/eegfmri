function [subj] = AG_interpret_weights_linreg(subj,results,results_IW,S, varargin)

% USAGE : [subj] = interpret_weights(subj,results,varargin)
%
% This algorithm is used to determine the level of influence
% of a given voxel to a given category. The value produced
% by the algorithm is the 'importance value' of the
% voxel. In the context of a neural network classifier, an
% input unit contains signal from a single voxel. The
% activation value of an input unit can be positive or
% negative.
%
%     imp_{i,j} = w_{i,j} * a_{i,j}      [Latex formatting]
%
% This eqn describes the calculation of this importance
% value. The subscript i is used to step across input units
% in the network. The subscript j is used to step across
% output units in the network. Thus, a single input unit has
% an importance value for each category. This importance
% value is defined as the product of two values, the weight
% between input unit i and output unit j, and the average
% activity of input unit i for the training patterns drawn
% from category j. The importance value is meant to reflect
% how effectively this input unit can alter the value of the
% corresponding output unit. The sign of the importance
% value reflects whether activity in this input unit tends
% to drive the activation value of the corresponding output
% unit up or down. If both wij and aij are positive (or both
% are negative), signal in the input unit causes the output
% unit to turn on; it has a positive importance. If the two
% have opposite sign, then the input unit¡Çs signal acts to
% turn off this output unit, and it has a negative
% importance value. Network classifiers use both positive
% and negative evidence to reduce error on the training set
% (Polyn et al., 2005). For example, if the signal in a
% voxel increases when a face is viewed, then a strong
% positive weight to the face output unit will reduce
% error. A strong negative weight from this input unit to
% the location output unit will also reduce error, by
% turning off the location output when face is the correct
% category.  (Explanation taken from
% http://polyn.com/struct/PolynDissertation.pdf)
%
% NOTE : This algorithm is designed to operate on a cell
% array, RESULTS, containing a separate individual 'results'
% structure. So if you've run the same cross-validation
% multiple times with a backprop classifier, each time
% initialized randomly, then you can average over these
% multiple attempts to get a better estimate of each voxel's
% contribution. If you feed in a single 'results' structure,
% it will be turned into a 1-cell cell array.
%
% In the end we can get two type importance
% maps one with both positive and negative canonicals and
% the other with only positive canonicals, depending on the
% 'type_canon' argument.
%
% VARARGIN :
% ---------
% TYPE_CANON:
% You can give an option to what type of importance maps you
% want. The default is both positive and negative canonical,
% but if you need only the positive then you can give this
% option in varargin as 'TYPE_CANON' as 'pos'.
%
% IMPMAP_NAME:
% An importance map a is cell array with either both the
% importance maps (all and positives) or only positives.This
% cell array will be added to the subj structure, as new
% patterns using the default group name : IMPMAP_NAME. So by
% default this names are going to be IMPMAP_ALL_1,
% IMPMAP_ALL_2.,etc.. based on the no. of iterations of N-1
% you have in your subj structure. But you can your use own
% name by using varargin 'impmap_name'.
%
% WHAT DOES THE STRUCTURE LOOK LIKE:
%
% The importance map will be an average over all the times
% you run your cross validation script. Thus if you have 8
% runs and 3 conditions and you run the cross validation 50
% times in the end you will end up averaging over all the 50
% iterations, and you will get 8 new maps holding a matrix
% of nVoxels (per N-1 iteration since it different for each
% iteration due to the no peeking anova.) by nConds (3 in
% this case).


% setting defaults
defaults.type_canon = 'all';
defaults.impmap_name = 'impmap';
args = propval(varargin,defaults);

% if RESULTS is not a cell array, then convert into a 1-cell
% cell array
if ~iscell(results)
    temp = results;
    results=[];
    results{1} = temp;
end

% number of times the cross_validation function has been run
nTimes = length(results);
% number of iterations of the N-1
nRuns = length(results{1}.iterations);

sanity_check(results,nTimes,nRuns,args);

% now getting the patterns to average
for i=1:nRuns
    % using RESULTS{1} on the assumption that all of the
    % RESULTS cells will be using the same object names
    pats_name = results{1}.iterations(i).created.patname;
    % since i need this for masked by field
    mask_name{i} = results{1}.iterations(i).created.maskname;
    sel_name  = results{1}.iterations(i).created.selname;
    reg_name  = results{1}.iterations(i).created.regsname;
    
    % really useful fuction Greg!
    masked_pats = get_masked_pattern(subj,pats_name,mask_name{i});
    selectors   = get_mat(subj,'selector',sel_name);
    regressors  = get_mat(subj,'regressors',reg_name);
    
    % this picks out the only the training patterns.
    %selectors ==2 <-- changed from 1 to select testing patterns only by AG
    %10/13/09
    
    pats_to_avg  = masked_pats(:,(selectors == 2));
    targs_to_avg = regressors(:,(selectors == 2));
    nConds       = size(regressors,1);
    
    % average the patterns per category.
    for j=1:nConds
        
        
        curData  = pats_to_avg;
        nVox_mean = mean(curData,2);
        
        %         if strcmp(args.type_canon,'pos')
        %             nVox_mean(find(nVox_mean < 0)) = 0;
        %         end
        iteration{i}.canonical(:,j) = nVox_mean;
        
        
        for it = 1:length(S.impType)
            %for cds = 1:length(S.condnames)
            importance.(S.impType{it}).(S.condnames{j}).map= zeros(size(iteration{1}.canonical(:,j)));
            
            %     importance{i}.pos.map= zeros(size(iteration{i}.canonical));
            %     importance{i}.neg.map= zeros(size(iteration{i}.canonical));
            %     importance{i}.both.map= zeros(size(iteration{i}.canonical));
            %end
        end
    end
end


cds = S.condnames{1};
% this step collapses across every condition per iteration
% per nTimes
for i=1:nTimes
    %for j=1:nRuns
    
        
        
        
        
        for j=1:nRuns
            curWts_h(:,j) = (results_IW{j}.iterations(1).scratchpad.net.IW{1}');
            curCanonical_h(:,j) = iteration{j}.canonical;
        end
       
        if nRuns==1
            curWts = curWts_h;
            curCanonical = curCanonical_h;
        else
            curWts = mean(curWts_h')';
            curCanonical = mean(curCanonical_h')';
        end
        
        
        % these are the weights of the voxels per iterations.
        %for k=1:nConds
        % for all the canonical
        % these are the weights per nTimes per nRuns per cat
        %curWts=results_IW{j}.iterations(1).scratchpad.net.IW{1}(cds,:)'; % JR: get weights for output unit 1
        %curWts{2}=results_IW{i}.iterations(j).scratchpad.net.IW{1}(2,:)'; % JR: get weights for output unit 2
        
        
        
        %curCanonical = iteration{j}.canonical(:,cds); % JR: get mean activity for output unit 1
        %curCanonical{2} = iteration{j}.canonical(:,2); % JR: get mean activity for output unit 2
        
        
        
%         Vox.pos.(cds).act = find(curCanonical>0);
%         Vox.pos.(cds).wts = find(curWts>0);
%         Vox.pos.(cds).sameSign = intersect(Vox.pos.(cds).wts , Vox.pos.(cds).act);
%         
%         Vox.neg.(cds).act = find(curCanonical<0);
%         Vox.neg.(cds).wts = find(curWts<0);
%         Vox.neg.(cds).sameSign = intersect(Vox.neg.(cds).wts , Vox.neg.(cds).act);
%         
%         %Vox.both.(cds).act = find(curCanonical);
%         %Vox.both.(cds).wts = find(curWts);
%         %Vox.both.(cds).sameSign = intersect(Vox.both.(cds).wts , Vox.both.(cds).act);
%         
%         Vox.both.(cds).act = 1: length(curCanonical);
%         Vox.both.(cds).wts = 1: length(curWts);
%         Vox.both.(cds).sameSign = intersect(Vox.both.(cds).wts , Vox.both.(cds).act);
%         
%         curImp.pos.(cds) = zeros(size(curWts));
%         curImp.neg.(cds) = zeros(size(curWts));
%         curImp.both.(cds) = zeros(size(curWts));
%         
%         curImp.pos.(cds)(Vox.pos.(cds).sameSign)  = curWts(Vox.pos.(cds).sameSign) .* curCanonical(Vox.pos.(cds).sameSign);
%         curImp.neg.(cds)(Vox.neg.(cds).sameSign) = curWts(Vox.neg.(cds).sameSign) .* curCanonical(Vox.neg.(cds).sameSign);
%         curImp.both.(cds)(Vox.both.(cds).sameSign) = curImp.pos.(cds) + curImp.neg.(cds);
%         
%         curImp.raw.(cds) = curWts.*curCanonical;
%                            
%         %Vox.(S.impType{it}).(cds).sameSign = intersect( Vox.(S.impType{it}).(cds).act, Vox.(S.impType{it}).(cds).wts);
%       
%         
%         importance.pos.(cds).map = importance.pos.(cds).map + curImp.pos.(cds);
%         importance.neg.(cds).map =  importance.neg.(cds).map + curImp.neg.(cds);
%         importance.both.(cds).map = importance.both.(cds).map + curImp.both.(cds);
%         
%         importance.raw.(cds).map = importance.raw.(cds).map + curImp.raw.(cds);
        %importance.raw(S.condnames.map = importance.raacds.map
        
        %
        %                 importance{j}.pos.map(:,1) = importance{j}.pos.map(:,1) + FinalImp.pos{1};
        %                 importance{j}.pos.map(:,2) =  importance{j}.pos.map(:,2) + FinalImp.pos{2};
        %                 importance{j}.neg.map(:,1) = importance{j}.neg.map(:,1) + FinalImp.neg{1};
        %                 importance{j}.neg.map(:,2) = importance{j}.neg.map(:,2) + FinalImp.neg{2};
        %                 importance{j}.both.map(:,1) = importance{j}.both.map(:,1) + FinalImp.both{1};
        %                 importance{j}.both.map(:,2) = importance{j}.both.map(:,2) + FinalImp.both{2};
    
end
%end


% we divide by the no. of times we add each block to
% get an average importance map per iteration of n-1 per condition
% and put everything into the subj structure.
for it = 1%:length(S.impType) %changed to write out only beta values b/c of linReg analysis.  AG 7/6/10.

        %for j=1:nRuns
%         importance.(S.impType{it}).(cds).map= importance.(S.impType{it}).(cds).map;
        
        patname = strcat(args.impmap_name,  cds );  %changed to write out only beta values b/c of linReg analysis.  AG 7/6/10.
        subj = init_object(subj,'pattern',patname);
        subj = set_mat(subj,'pattern',patname,curWts); %changed to write out only beta values b/c of linReg analysis.  AG 7/6/10.
        subj = set_objfield(subj,'pattern', patname ,'masked_by', mask_name{1});
        subj = set_objfield(subj,'pattern', patname ,'group_name',args.impmap_name);
        created.type_canonicals = args.type_canon;
        subj = add_created(subj,'pattern',patname,created);
    
end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = sanity_check(results,nTimes,nRuns,args)

% checking if all the names are ok.
pats_name = results{1}.iterations(1).created.patname;
mask_name = results{1}.iterations(1).created.maskname;
sel_name  = results{1}.iterations(1).created.selname;
reg_name  = results{1}.iterations(1).created.regsname;

if nTimes > 1
    for i=2:nTimes
        newpats_name = results{i}.iterations(1).created.patname;
        newmask_name = results{i}.iterations(1).created.maskname;
        newsel_name  = results{i}.iterations(1).created.selname;
        newreg_name  = results{i}.iterations(1).created.regsname;
        
        if ~strcmp( pats_name,newpats_name) || ~strcmp( mask_name,newmask_name)  ||  ~strcmp( sel_name,newsel_name)|| ~strcmp(reg_name,newreg_name)
            error('The names of your patterns, masks , selectors and regressors matrices must be the same for iterations');
        end
    end
end

for i=2:nTimes
    if ~isequal(nRuns,length(results{i}.iterations))
        error('The no. of runs in each of the cross validations results structures should be the same');
    end
end

% this is to check if it's a backpropagation classifier and
% if so then check if there are hidden units;
% for i=1:nTimes
%     if ~strcmp(results{i}.iterations(1).scratchpad.class_args.train_funct_name,'train_bp')
%         error('This only works with the backpropagation classifier');
%     end
%
%     if ~isequal(results{i}.iterations(1).scratchpad.class_args.nHidden,0)
%         error('No. of nHidden units must be 0');
%     end
% end

if ~(strcmp(args.type_canon,'pos') || strcmp(args.type_canon,'all'))
    error('type_canon can be only ''all'' or ''pos''');
end
%end
