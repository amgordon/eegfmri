% wrapper for ecogClassifyROI2, this scripts calls the function and changes
% the time parameter.


clearvars
close all

shufl = 0; % shuffle y vector (labels) to find out mean
%durations_vect = [zeros(10,1) (0.1:0.1:1)'];
%durations_vect = [0 1; 0 0.5; 0.5 1; 0.25 0.75; 0 0.25; 0.25 0.5; 0.5 0.75; 0.75 1];
durations_vect = [0 1];

subjs = { 'ef_040412'  'ef_040512' 'ef_040712' 'ef_040712_2'  'ef_041112' 'ef_042912' 'ef_050112'};

for d = 1:size(durations_vect,1)
    dur = durations_vect(d,:);
    eegClassifyROI(subjs,dur,shufl,{'hits','cr'});
end
