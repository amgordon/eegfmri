function [ output_args ] = EF_MakeBetaMaps(S, sess)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% TO DO:  HP FILTER THE DESIGN MATRIX.  

[S par] = EF_EEG_ResponseParams(subj_id);

if (S.use_premade_workspace&&existWorkspace) % if we are supposed to use premade workspace, and one with the correct name exists
    load(S.workspace, 'subj')
else
    [subj] = EF_mvpa_load_and_preprocess_raw_data(S);
end

xSPMOrig = load(fullfile(S.univar_dir, 'SPM.mat'));
xSPMOrig = xSPMOrig.SPM;

%create new U struct for trial-specific regressor
xSPMOrig.Sess.U(2:(end+1)) = xSPMOrig.Sess.U(1:(end));
xSPMOrig.Sess.U(1).P.name = 'none';
xSPMOrig.Sess.U(1).dur = 0;

k = xSPMOrig.nscan;
fMRI_T     = xSPMOrig.xBF.T;
fMRI_T0    = xSPMOrig.xBF.T0;
bf = xSPMOrig.xBF.bf;
V = xSPMOrig.xBF.Volterra;

pat = subj.patterns{end}.mat; %MAKE THIS MORE GENERAL
sess = xSPMOrig.Sess.C.C;
constant_col = ones(size(sess,1),1);

idx.sessTR = subj.selectors{1}.mat;

for a = 1:length(idx.allTrials);
    
    xSPM = xSPMOrig;    
    thisOnset = idx.allTrials(a);
    thisOnsetName = [S.betaMapPrefix prepend(num2str(a),4)];
    
    %remove thisOnset from its original location
    for u=2:length(xSPM.Sess.U);
        idxToDelete = ismember(xSPM.Sess.U(u).ons, thisOnset);
        xSPM.Sess.U(u).ons(idxToDelete) = [];
        xSPM.Sess.U(u).dur(idxToDelete) = [];
    end
    
    xSPM.Sess.U(1).ons = thisOnset;
    xSPM.Sess.U(1).dur = 0;
    xSPM.Sess.U(1).name = {thisOnsetName};
    
    %remove conditions with no onsets.
    idxToRemove = [];
    for u=2:length(xSPM.Sess.U);
        if isempty(xSPM.Sess.U(u).ons)
            idxToRemove = [u idxToRemove];
        end
    end
    xSPM.Sess.U(idxToRemove)=[];
    
    % Get inputs, neuronal causes or stimulus functions U
    U = spm_get_ons(xSPM,1);
    
    % Convolve stimulus functions with basis functions
    X = spm_Volterra(U,bf,V);
    
    % Resample regressors at acquisition times (32 bin offset)
    X = X([0:(k - 1)]*fMRI_T + fMRI_T0 + 32,:);
    
    betaMaps = {'1'};
    
    %replace the first regressor with a set of FIR regs
    if strcmp(S.bf, 'FIR')
        X(:,1) = [];
        
        xFIR = zeros(size(X,1),S.num_FIR_bins);
        
        for h=1:S.num_FIR_bins
            thisTR = h + floor(thisOnset/S.TR);
            
            % only add the regressor if its value is less
            % than the number of volumes in a session.
            if thisTR<=size(X,1)
                xFIR(thisTR, h)=1;
                betaMaps{h} = ['TR' num2str(h)];
            end
        end
        X = [xFIR X];   
    end
    
    % only put data and model from the session to which the
    % given trial belongs
    if S.modelSubSample==1
        
        idx.thisSess = (idx.sessTR==idx.sess(a));
        
        XMat = [X(idx.thisSess,:) constant_col(idx.thisSess)];
        
        %remove all-zero columns
        idx.notAllZero = sum(abs(XMat))~=0;
        XMat = XMat(:,idx.notAllZero);
        
        patMat = pat(:,idx.thisSess);
        betaMaps(~idx.notAllZero(1:length(betaMaps))) = [];
    else
        XMat = [X sess constant_col ];
    end
    
    % preallocate betas
    b = nan(size(XMat ,2), size(patMat,1));
    
    % do the glms
    for i = 1:size(patMat,1)
        thisPat = patMat(i,:);
        b(:,i) = regress(thisPat', XMat);
    end
     
    % write out the beta maps
    vol_info = S.vol_info;
    voxel_inds = subj.masks{end}.mat;
    
    for m = 1:length(betaMaps)     
        outMap=  zeros(vol_info.dim);
        outMap(voxel_inds==1) = b(m,:);
        
        vol_info.dir = S.patRegDir;
        vol_info.fname = [ vol_info.dir '/' thisOnsetName '_' betaMaps{m} '.img'];
        
        S.pat{i}.map{m}.fname = vol_info.fname;
        
        if isempty(dir([vol_info.dir]))
            mkdir(vol_info.dir);
        end
        
        spm_write_vol(vol_info,outMap);
        fprintf('\n wrote out %s \n', vol_info.fname);
    end    
end

end



