yC = load('LPS_withBOLD_CRs');
for i=1:7, for j = 1:8, for k=1:12, zC(i,j,k) =  yC.res{k}.EEGWithBOLD.BOLDPat(i).EEGBin(j).rObsVsPred; end, end, end
zC = shiftdim(zC,2);
[~,~,~,zC2] = ttest(zC);
zC2 = squeeze(zC2.tstat);

yH = load('LPS_withParietalBOLD_hits.mat');
for i=1:7, for j = 1:8, for k=1:12, zH(i,j,k) =  yH.res{k}.EEGWithBOLD.BOLDPat(i).EEGBin(j).rObsVsPred; end, end, end
zH = shiftdim(zH,2);
[~,~,~,zH2] = ttest(zH);
zH2 = squeeze(zH2.tstat);


yA = load('LPSLeftBOLDAcrossERP.mat');
for i=1:7, for j = 1:8, for k=1:11, zA(i,j,k) =  yA.res{k}.EEGWithBOLD.BOLDPat(i).EEGBin(j).rObsVsPred; end, end, end
zA = shiftdim(zA,2);
[~,~,~,zA2] = ttest(zA);
zA2 = squeeze(zA2.tstat);




yC = load('LPI_CRsWithLeftBOLDAcrossERP.mat');
for i=1:7, for j = 1:8, for k=1:12, zC(i,j,k) =  yC.res{k}.EEGWithBOLD.BOLDPat(i).EEGBin(j).rObsVsPred; end, end, end
zC = shiftdim(zC,2);
[~,~,~,zC2] = ttest(zC);
zC2 = squeeze(zC2.tstat);

yH = load('LPI_HitsWithLeftBOLDAcrossERP.mat');
for i=1:7, for j = 1:8, for k=1:12, zH(i,j,k) =  yH.res{k}.EEGWithBOLD.BOLDPat(i).EEGBin(j).rObsVsPred; end, end, end
zH = shiftdim(zH,2);
[~,~,~,zH2] = ttest(zH);
zH2 = squeeze(zH2.tstat);


yA = load('LPI_allTrialsWithLeftBOLDAcrossERP.mat');
for i=1:7, for j = 1:8, for k=1:11, zA(i,j,k) =  yA.res{k}.EEGWithBOLD.BOLDPat(i).EEGBin(j).rObsVsPred; end, end, end
zA = shiftdim(zA,2);
[~,~,~,zA2] = ttest(zA);
zA2 = squeeze(zA2.tstat);
