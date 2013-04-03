
% this script downsamples the data in order to make it more manageble for
% classification

% alex g. June 19 2012

clear all
close all

datapath = '/biac4/wagner/biac3/wagner5/alan/eegfmri/fmri_data/';
subjs = { 'ef_040412'  'ef_040512' 'ef_040712' 'ef_040712_2'  'ef_041112' 'ef_042912' 'ef_050112'};
ext = 'LP30Hz_zscored.mat';

for s=1:length(subjs)
    load([datapath subjs{s} '/erpData/trial_data_' ext]);
    
    S = [];
    S.subj = subjs{s};
    S.binsize = 50;
    S.sldwin = 50;
    [S.nchan, S.nevents, S.nsamps] = size(data.trialdata);
    
    % corresponding time points to the data
    time_points = round(linspace(data.dur(1)*1000,data.dur(2)*1000-1,S.nsamps));
    
    binpoints = time_points(1):S.sldwin:time_points(end)+1;
    % start point of each bin
    S.binStartsamps = binpoints(binpoints<=(time_points(end)+1- S.binsize));
    
    % end point of each bin
    S.binEndsamps = binpoints(binpoints>=(time_points(1) + S.binsize));
    
    % for reference, this are the processed bins are stored
    S.bins = [S.binStartsamps' S.binEndsamps'];
    
    % number of bins in the classification window
    S.nbins = size(S.bins,1);
    
    [~,S.timebin_allocation]=histc(time_points,binpoints);
    
    % create matrix with reduced data size
    S.bin_trials = zeros(S.nchan,S.nevents,S.nbins,'single');
    
    for bin = 1:S.nbins
        bin_logical = S.timebin_allocation==bin |S.timebin_allocation==bin+1;
        S.bin_trials(:,:,bin) = nanmean(data.ztrialdata(:,:,bin_logical),3);
    end
    
    save([datapath subjs{s} '/erpData/trial_data_downsamp_' ext], 'S')
end
