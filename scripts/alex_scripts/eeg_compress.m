function newdata = eeg_compress(data,factor,par)
% this function compresses data by factor using decimation
% inputs:
%           data:   matrix chan X samps, double format
%           factor: integer indicating compression factor
%           par:    flag indicatin if it should run in parallel
%
% outputs: 
%           newdata: matrix chan X nsamps/factor, created by default
%           decimation filter parameters, i.e. IIR, 8 order

% Alex Gonzalez July 27, 2012

if nargin <2
    error('Not enough inputs');
elseif nargin==2
    par =0; % default not to run in parallel
end

[nch,nsamps]= size(data);


newdata = zeros(nch,floor(nsamps/factor));
if ~par
    for c=1:nch
        newdata(c,:) = decimate(data(c,:),factor);
    end
else
    parfor c=1:nch
        newdata(c,:) = decimate(data(c,:),factor);
    end
end

