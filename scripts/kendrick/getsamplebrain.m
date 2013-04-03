function f = getsamplebrain(n)

% function f = getsamplebrain(n)
%
% <n> (optional) is
%   1 means 64 x 64 x 19 (2.5-mm isotropic) with values in [0,1]
%   2 means 256 x 256 x 16 (0.75 mm x 0.75 mm x 3 mm) with values in [0,3]
%   3 means 240 x 276 x 240 (1-mm isotropic) with values in [0,500]
%   4 means 64 x 64 x 20 x 50 (2.5-mm isotropic, 1337.702 ms TR) with values in [0,3000]
%   default: 1.
%
% return a sample brain(s).

% input
if ~exist('n','var') || isempty(n)
  n = 1;
end

% do it
f = double(loadmulti(strrep(which('getsamplebrain'),'getsamplebrain.m',sprintf('getsamplebrain%d.mat',n)),'samplebrain'));
