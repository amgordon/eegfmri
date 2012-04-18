function dataout = eegfmri_multibandpass_v3(data,filters,comp)
% eegfmri_multibandpass does multi bandpass filtering and compression
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eegfmri bandpass
% Version 3.0 (uses zero-phase IIR to compress)
% Created by Alex Gonzalez
% Stanford Memory Lab
% Jan 11, 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% use filtfilt with fir filters for now phase effects

% bandnames
bands = fieldnames(filters);
% number of bands
n = numel(bands);
% dataout
nsamps = ceil(length(data)/comp);
dataout = zeros(n-1,nsamps,'single');
for i = 2:n
    x = filtfilt(filters.(bands{i}).Numerator,1,double(data)); 
    dataout(i-1,:) = single(decimate(double(x),comp));
end
