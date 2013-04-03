% collect behavioral data into single file

clear all
close all

datapath = '/biac4/wagner/biac3/wagner5/alan/eegfmri/fmri_data/';

subjs = { 'ef_040412'  'ef_040512' 'ef_040712' 'ef_040712_2'  'ef_041112' 'ef_042912' 'ef_050112'};

fields = {'old','new','hits','cr','Remhits','HChits','LChits','LCcr','HCcr','FA','miss','dP','RT'};
oldfields = {'old','new','hits','cr','Rem_hits','HC_hits','LC_hits','LC_cr','HC_cr','FA','misses','dP','RT'};

nruns = 5;

for s = 1:numel(subjs)
    behdata =  [];
    behdata.subject = subjs{s};
    
    for f = 1:numel(fields)
        behdata.(fields{f}) =[];
    end
    
    for r = 1:nruns
        load([datapath subjs{s} '/erpData/r' num2str(r) '/behdata.mat']);
        for f = 1:numel(fields)
            behdata.(fields{f}) = [behdata.(fields{f}) ; data.(oldfields{f})];
        end 
    end
    
   save([datapath subjs{s} '/erpData/results/behdata.mat'], 'behdata') 
    
end