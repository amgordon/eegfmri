function [subj, S] = AG_balanceTrainPats(S, subj, reg_name)

%Get the regressor if interest
Reg = get_mat(subj,'regressors',reg_name);
regLength = length(Reg);

for i = 1:length(unique(Reg));
    init_cond{i} = shuffle(find(Reg(i,:)==1)); 
    remaining_cond{i} = init_cond{i};
end
% init_cond1 = shuffle(find(Reg(1,:)==1));
% init_cond2 = shuffle(find(Reg(2,:)==1));
init_total = [init_cond{:}];

% remaining_cond1 = init_cond1;
% remaining_cond2 = init_cond2;
remaining_total = init_total;

%%
if S.preserveRunsStruct
    ItLevels  = get_mat(subj, 'selector', 'meta_runs_condensed');
    S.runsSelector = 'meta_runs_condensed_balanced';
    
    subj = init_object(subj,'selector',S.runsSelector);
    subj = set_mat(subj,'selector',S.runsSelector, ItLevels);
elseif S.leaveOneTrialOut
    
    ItLevels = 1:regLength;
 
    S.runsSelector = 'leaveOneTrialOut';
    
    subj = init_object(subj,'selector',S.runsSelector);
    subj = set_mat(subj,'selector',S.runsSelector, ItLevels);

else
    %this block ensures that all iterations have equal(ish) ratios of each
    %class present in each test/train sets.  This maximizes the amount of
    %balanced training patterns we'll be allowed to include, and thus optimizes
    %our classificaiton power.
    for i = 1:S.xValIterations
        
        %number of total patterns in each iterations
        sz = round(i*regLength/S.xValIterations) - (round((i-1)*regLength/S.xValIterations));
        
        % proportion of cond1 patterns present in the remaining set
        cond_ratio = length(remaining_cond{1})/(length(remaining_total));
        
        NCondToTake{1} = round(cond_ratio*sz);
        NCondToTake{2} = sz - NCondToTake{1};
        
        for j = 1:2
            %number of cond1 and cond2 patterns to select
            
            
            %randomly chosen cond1 and cond2 patterns
            condTaken{j} = remaining_cond{j}(1:NCondToTake{j});
            
            %update the sets of patterns that have yet to be selected
            remaining_cond{j} = shuffle(setxor(remaining_cond{j}, condTaken{j}));
            %remaining_cond2 = shuffle(setxor(remaining_cond2, cond2Taken));
        end
        
        remaining_total = [remaining_cond{1}, remaining_cond{2}];
        
        %store the trials selected for the given iteration
        itTrials{i} = [condTaken{1}, condTaken{2}];
        
        %the corresponding iteration
        w{i} = i*ones(1,sz);
        
    end
    
    
    %concatenate
    Iters = [itTrials{:}];
    
    %create index linking chronilogical order to iteration
    [~,ix_iters] = sort(Iters);
    
    %concatenate
    wConcat = [w{:}];
    
    %which iteration each pattern is assigned to.  patterns remain in
    %chronoligcal order.
    ItLevels = wConcat(ix_iters);
    S.runsSelector = 'shuffledRuns';
    
    subj = init_object(subj,'selector',S.runsSelector);
    subj = set_mat(subj,'selector',S.runsSelector, ItLevels);
end

% save mat into the subj selector "shuffled runs"



%%
%this block of code creates selectors for each iteration.  These selectors
%insure that there are equal numbers of each condition represented in the
%training set.
k = 0;
for i = S.ItLevelConstrain;
    for j = 1:S.numTrainIts;
        
        k = k+1;
        
        actives = ones(size(ItLevels));
        
        if S.ItLevelConstrain==0
            levelNum = j;
        else
            levelNum = i;
        end
        
        %all train patterns for a given iteration
        train_trials_w = (ItLevels~=levelNum);
        
        trainReg1 = train_trials_w.*Reg(1,:);
        trainReg2 = train_trials_w.*Reg(2,:);
        
        size_cond1 = sum(trainReg1==1);
        size_cond2 = sum(trainReg2==1);
        
        %the number of patterns in the minority class
        regs_size = min([size_cond1 size_cond2]);
        
        v1 = shuffle(find(trainReg1));
        v2 = shuffle(find(trainReg2));
        
        %cut patterns fromthe majority class
        toCut = [v1(1:(size_cond1 - regs_size)), v2(1:(size_cond2 - regs_size))];
        
        actives(toCut) = 0;
        
        %cell array of iteration-specific actives selectors
        itActives{k} = actives;
    end
end


%create subj selectors

for i =1:length(itActives)
    
    if S.onlyOneIt
        unBalancedMat = get_mat(subj, 'selector', 'meta_runs_condensed');
    else
        unBalancedMat = get_mat(subj, 'selector', ['meta_runs_condensed_xval_' num2str(i)]);
    end
    
    
    thisActiveIter = ['activeIter_' num2str(i)];
    
    subj = init_object(subj,'selector',thisActiveIter);
    subj = set_mat(subj,'selector',thisActiveIter, unBalancedMat.*itActives{i});
    
    subj = set_objfield(subj,'selector',thisActiveIter,'group_name', 'activeIter');
    
end

S.thisSelector = [S.runsSelector '_xval'];
    
    