function [out] = EF_FFTConfidenceAnalysis(par)

NIts = 1;
fs = 1/5.5195; %this is the presentation frequency

% read in retrieval data
[res , idx] = EF_BehAnalyzer(par);

% study data
studyDataStr = dir(fullfile(par.behavdir, '*Study*'));
studyData = load(fullfile(par.behavdir, studyDataStr.name));
studyItems = studyData.respSelData.item;

% map retrieval data to study data
studyConf = nan(size(studyItems));
for i = 1:length(studyItems)
   [~ , ix] = ismember(studyItems{i}, idx.item);
   
   % if no retrieval response exists for the item, label it '0'
   if ix==0
       studyConf(i) = 0;
   else
       studyConf(i) = idx.conf(ix);
   end
end

% replace '0' values with the mean of non-zero values
studyConf(studyConf==0) = mean(studyConf(studyConf~=0));
studyConf = zscore(studyConf);

out.studyConf = studyConf;
%% Fourier Transform analysis
n=length(studyConf);
nEff = 2^nextpow2(n);

Y = fft(studyConf, nEff); %fft the data

nyquist = 1/2;
out.fftpower = abs(Y(1:floor(nEff/2))).^2; %power in data
out.fftfreq = (fs/2)*(1:nEff/2)/(nEff/2)*nyquist; % frequency


for i=1:NIts
    Ys(i,:) = fft(shuffle(studyConf), nEff); % fft shuffled data
end

fftPowerShuff = abs(Ys(:,1:floor(nEff/2))).^2;
out.fftshuffPower = median(fftPowerShuff,1)'; % power in shuffled data

%permutation testing to determine where on the shuffled distribution the
%fft power lies
for k = 1:length(out.fftpower)
    pF(k) =  mean(out.fftpower(k)>fftPowerShuff(:,k));
end
pF(pF==1) = 1 - 1/NIts;
pF(pF==0) = 1/NIts;
out.fftZStat = norminv(pF);
%figure; plot(out.fftfreq, out.fftpower, 'r'); hold on;
%plot(out.fftfreq, out.fftshuffPower, 'b');


%% Power Analysis with specified bands
N = 10; %related to the number of neighboring points to use for lp filtering
Nhp = 16; %related to the number of neighboring points to use for hp filtering

freqconv  = 2/fs; 
%freqbounds = logspace(log10(.005), log10(fs/2), 11); %define frequency bands
freqbounds = logspace(log10(.005), log10(.04), 11); %define frequency bands

for i=1:length(freqbounds) % cycle through frequency bands
    
    x = studyConf;
        
    if i==length(freqbounds)
        lowpasscut = .04;
        highpasscut = .005;
    else
        lowpasscut = freqbounds(i+1);
        highpasscut = freqbounds(i);
    end
        
    %make lp filter
    f0 = lowpasscut*freqconv;
    df = (f0/20)*freqconv; 
    lp=firls(N,[0 f0-df/2 f0+df/2 1],[1 1 0 0]);
    
    %make hp filter
    f0 = highpasscut*freqconv;
    df = (f0/20)*freqconv; 
    w = [1e-5 1e3];
    hp=firls(Nhp,[0 f0-df/2 f0+df/2 1],[0 0 1 1],w);
    
    %filter real data
    x = filtfilt(hp,1,x);
    x = filtfilt(lp,1,x);
    
    %extract power for real data
    pwr_h = hilbert(x')';
    out.pwr(i,:) = abs(pwr_h).^2;
    out.pwrFreq(i) = mean([highpasscut lowpasscut]);
    
    for j=1:NIts
        xShuff = shuffle(studyConf);
        
        %filter shuffled data
        xShuff = filtfilt(hp,1,xShuff);
        xShuff = filtfilt(lp,1,xShuff);
        
        %extract power for shuffled data
        pwrShff_h(j,:) = hilbert(xShuff')';
    end
    pwrShff_h = abs(pwrShff_h).^2; %extract power;
    medianShuffledPwr = median(pwrShff_h,1);
    out.shuffledPwr(i,:) = medianShuffledPwr;
    
    % using permutation testing, find out where on the shuffled
    % distribution the power at each time point lies.
    for k = 1:size(out.pwr,2)
       p(k) =  mean(out.pwr(i,k)>pwrShff_h(:,k));
    end
    p(p==1) = 1 - 1/NIts;
    p(p==0) = 1/NIts;
    out.pwrZStat(i,:) = norminv(p);
end

%figure;
%imagesc(out.pwrZStat(:,21:140));


