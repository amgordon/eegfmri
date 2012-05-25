subjs = { 'ef_040412'  'ef_040512' 'ef_040712' 'ef_040712_2'  'ef_041112' 'ef_042912' 'ef_050112'};


for s=1:length(subjs)
    y = load(['/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/' subjs{s} '/erpData/results/trial_data_LP30Hz_RT.mat']);
    
    badch = [];
    
    z = y.data.trialdata(:,:);
    
    stdThresh = 3;
    channelIdx = 1:256;
    x = z*z';
    
    variances = diag(abs(x));
    orig_variances = variances;
    %iteratively remove channels that exceed the std threshold.
    while true
        theseBadCh=[];
        
        mx = mean(variances);
        
        sx = std(variances);
        
        theseBadCh=find(variances>=(mx+stdThresh*sx))';
        
        if isempty(find(theseBadCh, 1))
            break;
        end
        
        badch = [badch channelIdx(theseBadCh)];
        
        variances(theseBadCh) = [];
        channelIdx(theseBadCh) = [];
    end
    
    figure; semilogy(orig_variances);
    hold on;  
    semilogy(1:256,stdThresh*sx*ones(256,1)+mx*ones(256,1),'r')
    
    badChannel{s} = badch;
end


for s=1:length(subjs)
    y = load(['/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/' subjs{s} '/erpData/results/trial_data_LP30Hz_RT.mat']);
    
    badtr = [];
    
    z = permute(y.data.trialdata,[2 1 3]);
    cz = z(:,:);
    
    stdThresh = 5;

    x = cz*cz';
    
    variances = diag(abs(x));
    nTrials = length(variances);
    trialIdx = 1:nTrials;
        
    orig_variances = variances;
    
    badtr = 80:80:nTrials;
    %iteratively remove channels that exceed the std threshold.
    while true
        theseBadTr=[];
        
        mx = mean(variances);
        
        sx = std(variances);
        
        theseBadTr=[find(variances>=(mx+stdThresh*sx))'  find(variances<=(mx-stdThresh*sx))'];
        
        if isempty(find(theseBadTr, 1))
            break;
        end
        
        badtr = [badtr trialIdx(theseBadTr)];
        
        variances(theseBadTr) = [];
        trialIdx(theseBadTr) = [];
    end
    
 
    figure; semilogy(orig_variances);
    ylim([10^8 10^12]);
    hold on;  
    semilogy(1:nTrials,stdThresh*sx*ones(nTrials,1)+mx*ones(nTrials,1),'r')
    
    badTrial{s} = badtr;
end





% imagesc(x)
% x
% abs(x)
% log10(abs(x))
% imagesc(log10(abs(x)))
% plot(diag(log10(abs(x))))
%badch = [89,229,256];
% 
% z = y.trialdata(:,:);
% z = y.data.trialdata(:,:);
% z(badch,:) = 0;
% x = z*z';
% plot(diag(log10(abs(x))))
% mx = mean(diag(log10(abs(x))));
% sx = std(diag(log10(abs(x))));
% hold on;
% plot(1:256,[mx mx])
% plot(1:256,mx*ones(256,1))
% mx
% x(badch,badch)=0;
% mx = mean(diag(log10(abs(x))));
% mx
% % x
% 
% %mx
% % logplot(diag(x))
% % logyplot(diag(x))
% % help plot
% semilogy(diag(x))
% semilogy(1:256,mx*ones(256,1))
% sx = std(diag(log10(abs(x))));
% sx = std(diag((abs(x))));
% semilogy(1:256,sx*ones(256,1),'r')
% semilogy(1:256,2*sx*ones(256,1),'r')
% semilogy(1:256,2*sx*ones(256,1)+mx*ones(256,1),'r')
% abs(x)>=2*sx*ones(256,1)+mx*ones(256,1)
% diag(abs(x))>=2*sx*ones(256,1)+mx*ones(256,1)
% find(diag(abs(x))>=2*sx*ones(256,1)+mx*ones(256,1))
% dummy=find(diag(abs(x))>=2*sx*ones(256,1)+mx*ones(256,1))';
% badch
% badch=[badch dummy];
% badch
% z = y.data.trialdata(:,:);
% size(z)
% z = perm(y.data.trialdata,[2 1 3]);
% z = permutate(y.data.trialdata,[2 1 3]);
% z = permute(y.data.trialdata,[2 1 3]);
% size(z)
% cz = z(:,:);
% size(cz)
% x = cz*cz';
% size(x)
% imagesc(x)
% imagesc(log10(abs(x)))
% plot(diag(log10(abs(x))))
% cz(80:80:320,:)=0;
% x = cz*cz';
% imagesc(log10(abs(x)))
% plot(diag(log10(abs(x))))
% badtrials = 80:80:320;
% 
% badch{21} = [89   229   256     1    26    27    54    91   226   231   233   236   238   239   240]