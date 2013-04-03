function [freqOut] = EF_batchshebang(sa, varargin)



if (nargin < 1)
    sa = {'ef_091211' 'ef_091511' 'ef_092111' 'ef_092211'  'ef_092711' ...
        'ef_092911' 'ef_100511' 'ef_101411' 'ef_040512' 'ef_040712' ...
        'ef_041112' 'ef_042912'};
end
    
for s = 1:length(sa)
    
    
    par = EF_Params(sa(s));
    
    par = propval(varargin,par);
%    [dat.res(s) dat.idx(s)] = EF_BehAnalyzer(par);
    if par.goodSub && (~isempty(dir(fullfile(par.behavdir, '*Study*'))));
       %EF_wholeshebang(par, 'j');
       %EF_MakeBetaMaps(sa{s});
       freqOut(s) = EF_FFTConfidenceAnalysis(par);
    end
    %hasHires = exist(par.hiresimg);
    %
    %[~, r{s}] = EF_EEG_ResponseFunction(sa{s});
    %par = EF_Params(sa{s});
    %EF_wholeshebang(par, 'h');
    %EF_makeAnas(par)
end
