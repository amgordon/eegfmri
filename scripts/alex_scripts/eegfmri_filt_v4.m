function dataout = eegfmri_filt_v4(data,fs,lowpasscut,highpasscut,notch,par,nworkers)

% eegfmri_filt provides highly efficient High Pass, Lowpass and Notch
% filters, applied in that order. It uses FIR zero phase filtering to avoid phase
% lags. All filterss are performed by a least squares FIR.
%
% These filters were optimized for:
% lowpasscut = 125;
% highpasscut = 1;
% notch = 60;
% fs = 500;
%
% par is a flag that indicates if you want run the process in parallel
% default is 1, i.e. it will create 3 matlab workers to run the processs in
% parallel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eegfmri filt
% Version 4.0
% Created by Alex Gonzalez
% Stanford Memory Lab
% Feb 13, 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%conversion to normalized frequency
freqconv  = 2/fs;
% filter order
N = 200;
Nhp = 500;

if ~exist('par','var')
    par =1;
    nworkers = 3;
end
%% create filter for high pass
% uses interpolated FIR filter
if exist('highpasscut','var')
    if ~isempty(highpasscut)
        f0 = highpasscut*freqconv;
        df = f0*freqconv;
        w = [1e-5 1e3];
        hp=firls(Nhp,[0 f0-df/2 f0+df/2 1],[0 0 1 1],w);
    else
        hp =[];
    end
else
    hp =[];
end
%% create filter for low pass
% uses FIR least squares filter
if exist('lowpasscut','var')
    if ~isempty(lowpasscut)
        f0 = lowpasscut*freqconv;
        df = (f0/10)*freqconv;
        lp=firls(N,[0 f0-df/2 f0+df/2 1],[1 1 0 0]);
    else
        lp =[];
    end
else
    lp =[];
    
end
%% create filter for notch filter
% uses FIR least squares filter
if exist('notch','var')
    if ~isempty(notch)
        f0 = notch*freqconv;
        df = (f0/5)*freqconv;
        h=firls(N,[0 f0-df/2 f0+df/2 1],[1 1 0 0]);
        nf = 2*h - ((0:N)==(N/2));
    else
        nf = [];
    end
else
    nf = [];
end
%% use filtfilt with fir filters for now phase effects

% assume data is by row
n = size(data,1);
dataout = zeros(size(data));
if par
    isOpen = matlabpool('size') > 0;
    
    if isOpen
        matlabpool close
    end
    
 %   t = tic;
    matlabpool('open', nworkers)
    parfor i = 1:n
        %display(sprintf('Filtering Channel %i',i))
        x = data(i,:);
        %highpass
        if ~isempty(hp)
            x = filtfilt(hp,1,x);
        end
        % lowpass
        if ~isempty(lp)
            x = filtfilt(lp,1,x);
        end
        % notch
        if ~isempty(nf)
            x = filtfilt(nf,1,x);
        end
        dataout(i,:) = x;
        
    end
    matlabpool close
    
else
    for i = 1:n
        %display(sprintf('Filtering Channel %i',i))
        x = data(i,:);
        %highpass
        if ~isempty(hp)
            x = filtfilt(hp,1,x);
        end
        % lowpass
        if ~isempty(lp)
            x = filtfilt(lp,1,x);
        end
        % notch
        if ~isempty(nf)
            x = filtfilt(nf,1,x);
        end
        dataout(i,:) = x;
        
    end
end
%display(sprintf('Time to filter data %g seconds',toc(t)))
return


