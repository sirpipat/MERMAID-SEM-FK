function s = float2filestr(x, numsig, sep)
% s = FLOAT2FILESTR(x, numsig, sep)
%
% Converts a float to a string for file naming. By default, it replaces the
% dot to letter p to prevent confusion as decimal points may appear like
% file extensions.
%
% INPUT:
% x             value
% numsig        number of decimal points
% sep           separation string
%
% OUTPUT:
% s             string
%
% EXAMPLE:
% % simple call: equivalent to replace(sprintf('%.2f', pi), '.', 'p')
% % which returns '3p14'
% s = float2filename(pi);
%
% % You can change the number of decimal points and separation
% % This returns '3_1416'
% s = float2filename(pi, 4, '_');
%
% Last modified by sirawich-at-princeton.edu, 04/23/2025

defval('numsig', 2)
defval('sep', 'p')

s = replace(sprintf(sprintf('%%.%df', numsig), x), '.', sep);
end