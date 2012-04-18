

chnls = [];
chnls.all = 1:256;

% occipital
chnls.occ.LOS = [107 108 115:117 123 124];
chnls.occ.ROS = [139 149:151 158:160];

% parietal
chnls.parietal.LPI = [75 76 84 85 96 97];
chnls.parietal.LPS = [77,78,86:88,98,99];
chnls.parietal.PM = [100 101 110 118 119 127:129];
chnls.parietal.RPS = [141 142 152:154 162 163];
chnls.parietal.RPI = [161 170:172 179 180];

% central
chnls.central.LCI = [51 58 59 65 66 72];
chnls.central.LCS = [17 43 44 52 53 60];
chnls.central.CM = [9 45 81 132 186];
chnls.central.RCS = [144 155 184 185 197 198];
chnls.central.RCI = [ 164 173 182 183 195 196];

% frontal
chnls.frontal.LFI = [40 48:50 56 57];
chnls.frontal.LFS = [23 24 29 30 36 41 42];
chnls.frontal.FM = [7 14:16 22];
chnls.frontal.RFS = [5 6 26 207 214 215 224];
chnls.frontal.RFI = [204 205 212 213 222 223];

% temporal
chnls.temporal.LTPI = [92 93 102 103 111 112 120];
chnls.temporal.RTPI = [187 188 199:201 208 209];

% polar
chnls.polar.LFP = [34 37:39 46 47];
chnls.polar.FPM = [19 20 25:27 32 33];
chnls.polar.RFP = [2:3 10:12 18];

