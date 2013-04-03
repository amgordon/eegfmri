function llbands = makellbands(smoothingSzInMM, x, y, sliceThickness, time)
%
% llbands = llbandsMM2Voxel(smoothingSzInMM, x, [y], [sliceThickness], [time])
% 
% Set smoothing sizes for 4-D local linear regression in units of voxels
% and scans. Output arg, llbands, has 4 columns: 3 spatial vallues and one
% time value. It has as many rows as there are smoothing sizes (equal in
% length to first input arg).
%
% <smoothingSzInMM> : a vector of smoothing sizes in millimeters, assumed
%                       to be isotropic in 3d
% <x>               : pixel size in mm (x direction)
% <y>               : pixel size in mm (y direction; usu same as x)
% <sliceThickness>  : slice thickness in mm
% <time>            : a scalar or vector for smoothing in time (across
%                     scans)
%
% Example 1:
%
%  smoothingSzInMM   = 15 * 2.^([-2:3 Inf]);
%  x                 = 2.5;
%  llbands = makellbands(smoothingSzInMM, x);
%
% Example 2:
%
%  smoothingSzInMM   = 15 * 2.^(-2:3);
%  x                 = 3;
%  y                 = 3;
%  sliceThickness    = 4;
%  llbands = makellbands(smoothingSzInMM, x,y, sliceThickness);
%
% Example 3:
%
%  smoothingSzInMM   = 15 * 2.^(-2:3);
%  x                 = 3;
%  y                 = 3;
%  sliceThickness    = 4;
%  time              = 2;
%  llbands = makellbands(smoothingSzInMM, x,y, sliceThickness, time);
%
% This function was written by Jon Winawer.

% check input args; make defaults if missing
if nargin < 5, time = inf(1, length(smoothingSzInMM)); end
if nargin < 4, sliceThickness = x; end
if nargin < 3, y = x; end

% if smoothing in time is a scalar
if length(time) == 1, time = zeros(1, length(smoothingSzInMM)) + time; end
    
llbands = inf(length(smoothingSzInMM), 4);

% time
llbands(:,4) = time;

% space
for ii = 1:length(smoothingSzInMM)
    llbands(ii,1:3) = smoothingSzInMM(ii) ./ ([x y sliceThickness]);
end

end