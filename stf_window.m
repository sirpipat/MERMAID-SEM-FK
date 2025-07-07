function w = stf_window(t)
% w = STF_WINDOW(t)
%
% Makes a tapering window for preparing source-time function to be injected
% to a FK-SPECFEM3D run.
%
% INPUT:
% t         time
%
% OUTPUT:
% w         window, value is between 0 and 1
%                  t < -12,     w = 0
%           -12 <= t < -10,     w = cosine taper: 0->1
%           -10 <= t <  80,     w = 1
%            80 <= t < 100,     w = cosine taper: 1->0
%           100 <= t      ,     w = 0
%
%           This window works for twindow = 100 in FKMODEL file
%
% SEE ALSO:
% MAKEFKMODEL
%
% Last modified by sirawich-at-princeton.edu, 07/07/2025

w = and(t >= -12, t < -10) .* (1 + cos(pi * (t - (-10)) / ...
    ((-10) - (-12)))) / 2 + and(t >= -10, t < 80) + ...
    and(t >= 80, t < 100) .* (1 + cos(pi * (t - 80) / (100 - 80))) / 2;
end