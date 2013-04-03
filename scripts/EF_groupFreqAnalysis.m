function [freqOut] = EF_groupFreqAnalysis(sa, varargin)

timePointsToDiscard = 20;
fftPlotSmoothKernel = 20;

if (nargin < 1)
      sa = {'ef_092111' 'ef_092211'  ...
        'ef_092911' 'ef_100511' 'ef_101411' 'ef_040512'  ...
        'ef_041112' 'ef_042912' 'ef_111112_1'  'ef_111112_3'  'ef_111412_2' ...
        'ef_111112_2'   'ef_111512'};
end
 

for s = 1:length(sa)
    
    par = EF_Params(sa(s));    
    par = propval(varargin,par);

    freqOut(s) = EF_FFTConfidenceAnalysis(par); %analyze each sujbect
end


%% trialwise power in specified bands
pwr_h = cat(3,freqOut(:).pwr); %power for time, frequency, and subjects
pwr_clipped = pwr_h(:,(timePointsToDiscard+1):(end-timePointsToDiscard),:); %remove first N and last N timepoints

pwrShuffled_h = cat(3,freqOut(:).shuffledPwr); %power for time, frequency, and subjects
pwrShuffled_clipped = pwrShuffled_h(:,(timePointsToDiscard+1):(end-timePointsToDiscard),:); %remove first N and last N timepoints

%t-test of each power value for each time/frequency point, to see whether
%real and shuffled data differ
for i=1:size(pwr_clipped,1)
    for j=1:size(pwr_clipped,2)
        [~, ~, ~, a4] = ttest(squeeze(pwr_clipped(i,j,:)), squeeze(pwrShuffled_clipped(i,j,:)));       
        tMap1_h(i,j) = a4.tstat;        
    end 
    tMap1(i,:) = smooth(tMap1_h(i,:),160,'loess');
end
%figure;
%imagesc(tMap1);


%% mean power across trials in specified bands
%mean power across time
meanPwrAcrossTime = squeeze(mean(pwr_clipped,2));
meanPwrShuffledAcrossTime = squeeze(mean(pwrShuffled_clipped,2));

%t-test to see whether the mean power across time differs for the real vs.
%shuffled data
for i=1:size(pwr_clipped,1)
    [~,~,~,a4] = ttest(meanPwrAcrossTime(i,:), meanPwrShuffledAcrossTime(i,:));
    tMap2(i) = a4.tstat;
end
figure; 
imagesc(tMap2');

figure; plot(freqOut(1).pwrFreq(1:(end-1)), tMap2(1:(end-1)));
ylim([0 3])

%% FFT analyses
% FFT power for real and shuffled data
allFFT = [freqOut(:).fftpower];
allFFTShuffled = [freqOut(:).fftshuffPower];

freqs  = freqOut(1).fftfreq;
% test whether FFT power differs for real and shuffled data
for i=1:size(allFFT,1)
    [~, ~, ~, a4] = ttest(allFFT(i,:), allFFTShuffled(i,:));
    tMap3(i) = a4.tstat;
end

sw = sliding_window(128, 16, 8);

for i=1:size(sw, 2);
    theseFreqs = sw(:,i);
    freqLabels(i) = mean(freqs(theseFreqs));
    [~, ~, ~, a4] = ttest(mean(allFFT(theseFreqs,:)), mean(allFFTShuffled(theseFreqs,:)));
    tMap4(i) = a4.tstat;
end

figure; plot(freqLabels, tMap4)

%mean power across subjects
meanFFT = mean([freqOut(:).fftpower],2);
meanFFTShuffled = mean([freqOut(:).fftshuffPower],2);

% plot fft power for real and shuffled data
figure; 
plot(freqs, smooth(meanFFT, fftPlotSmoothKernel ), 'r');
hold on
plot(freqs, smooth(meanFFTShuffled, fftPlotSmoothKernel), 'b');



doc sprintf('end')
end


