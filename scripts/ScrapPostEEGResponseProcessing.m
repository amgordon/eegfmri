y1 = load('canonicalHRF_s1to5.mat');
y2 = load('canonicalHRF_s8.mat');
y3 = load('canonicalHRF_s9.mat');
y4 = load('canonicalHRF_s10.mat');
y5 = load('canonicalHRF_s11.mat');
y6 = load('canonicalHRF_s12.mat');

y1.res{8} = y2.res{1};
y1.res{9} = y3.res{1};
y1.res{10} = y4.res{1};
y1.res{11} = y5.res{1};
y1.res{12} = y6.res{1};

y1.S{8} = y2.S{1};
y1.S{9} = y3.S{1};
y1.S{10} = y4.S{1};
y1.S{11} = y5.S{1};
y1.S{12} = y6.S{1};

for i=1:11, 
    for j = 1:8
        w(i,j) = y.res{i}.EEGWithBOLD.BOLDPat.EEGBin(j).rObsVsPred;
    end
end

[a1 a2 a3 a4] = ttest(w);
figure; plot(a4.tstat)
set(gca,'XTickLabel',{ -125 125 375 625 875 1125 1375 1625})