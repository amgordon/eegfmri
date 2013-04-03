function hrf = getcanonicalhrf(duration,tr)

% function hrf = getcanonicalhrf(duration,tr)
%
% <duration> is the duration of the stimulus in seconds.
%   must be a multiple of 0.1.
% <tr> is the TR in seconds.
%
% average empirical HRFs from various datasets (where the HRFs were
% in response to a 3-s stimulus) and then do the requisite 
% deconvolution and convolution to obtain a predicted HRF to a 
% stimulus of duration <duration>, with data sampled at a TR of <tr>.
%
% the resulting HRF is a row vector whose first point is 
% coincident with stimulus onset.  the HRF is normalized such 
% that the maximum value is one.  note that if <duration> is 
% small, the resulting HRF will be quite noisy.
%
% example:
% hrf = getcanonicalhrf(4,1);
% figure; plot(0:length(hrf)-1,hrf,'ro-');

% load HRFs from five datasets and then take the average.
% these were the empirical response to a 3-s stimulus, TR 1.323751 s
hrf = mean(catcell(2,getsamplehrf([9 10 11 12 14],1)),2)';  % 1 x time
trorig = 1.323751;

% resample to 0.1-s resolution
trnew = 0.1;
hrf = interp1((0:length(hrf)-1)*trorig,hrf,0:trnew:(length(hrf)-1)*trorig,'cubic');

% deconvolve to get the predicted response to 0.1-s stimulus
hrf = deconvolvevectors(hrf,ones(1,3/trnew));

% convolve to get the predicted response to the desired stimulus duration
hrf = conv(hrf,ones(1,duration/trnew));

% resample to desired TR
hrf = interp1((0:length(hrf)-1)*trnew,hrf,0:tr:(length(hrf)-1)*trnew,'cubic');

% make the peak equal to one
hrf = hrf / max(hrf);
