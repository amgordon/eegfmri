function str = prepend(input, numChars, filler, append)
% Written by J. Drucker 09/06/07
% str = prepend(input, numChars, filler, append)
% Pads the front of string "input" with char "filler" until it is of length
% "numChars". If boolean "append" is flagged, the padding will be appended
% instead of prepended.
%
% Defaults:
% numChars == 2
% filler == '0'
% append == 0

switch(nargin)
    case 1,
        numChars = 2;
        filler = '0';
        append = 0;
    case 2,
        filler = '0';
        append = 0;
    case 3,
        append = 0;
end


% Error if input is longer than numChars
if numel(input) > numChars
    error('Input string is too long');
end

padstr = zeros(1, numChars - numel(input)) + filler;

if append
    str = [input padstr];
else
    str = [padstr input];
end