function f = checkmemoryworkspace

% function f = checkmemoryworkspace
% 
% figure out the number of megabytes used in the workspace of
% the caller of this function and return that number.
%
% example:
% a = zeros(10000,1000);
% checkmemoryworkspace

a = evalin('caller','whos');
f = sum(cat(2,a.bytes))/1024/1024;
