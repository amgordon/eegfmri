function chnls = channel_rois()

chnls = [];

% occipital
chnls.LOS = [107 108 115:117 123 124];
chnls.ROS = [139 149:151 158:160];

% parietal
chnls.LPI = [75 76 84 85 96 97];
chnls.LPS = [77,78,86:88,98,99];
chnls.PM = [100 101 110 118 119 127:129];
chnls.RPS = [141 142 152:154 162 163];
chnls.RPI = [161 170:172 179 180];

% central
chnls.LCI = [51 58 59 65 66 72];
chnls.LCS = [17 43 44 52 53 60];
chnls.CM = [9 45 81 132 186];
chnls.RCS = [144 155 184 185 197 198];
chnls.RCI = [ 164 173 182 183 195 196];

% frontal
chnls.LFI = [40 48:50 56 57];
chnls.LFS = [23 24 29 30 36 41 42];
chnls.FM = [7 14:16 22];
chnls.RFS = [5 6 26 207 214 215 224];
chnls.RFI = [204 205 212 213 222 223];

% temporal
chnls.LTPI = [92 93 102 103 111 112 120];
chnls.RTPI = [187 188 199:201 208 209];

% polar
chnls.LFP = [34 37:39 46 47];
chnls.FPM = [19 20 25:27 32 33];
chnls.RFP = [2:3 10:12 18];

chnls.allroinames = fieldnames(chnls);

% classify these ones
%chnls.roinames ={'LPI','LPS','PM','RPS','RPI','LFI','LFS','FM','RFS','RFI'};
chnls.roinames = {'allchnls'};
chnls.allchnls = 1:256;
