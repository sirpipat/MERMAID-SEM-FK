function [Ypk,Xpk,Wpk,Ppk] = findpeakstopbottom(Yin,varargin)
% [Ypk,Xpk,Wpk,Ppk] = FINDPEAKSTOPBOTOM(Yin,varargin)
%
% Finds the locations and values at both upright and upside-down local peaks.
%
% INPUT:
% Yin           input data vector
%
% OUTPUT:
% Ypk           the peak values
% Xpk           the peak locations: either indices or the corresponding X
%               when Xin is specified
% Wpk           the peak widths
% Ppk           the peak prominence
% SEE ALSO:
% FINDPEAKS
%
% Last modified by sirawich-at-princeton.edu, 05/02/2025

% top peaks
[Ypk_top,Xpk_top,Wpk_top,Ppk_top] = findpeaks(Yin,varargin{:});

% bottom peaks
[Ypk_bot,Xpk_bot,Wpk_bot,Ppk_bot] = findpeaks(-Yin,varargin{:});

% merge
Ypk = [Ypk_top; -Ypk_bot];
Xpk = [Xpk_top; Xpk_bot];
Wpk = [Wpk_top; Wpk_bot];
Ppk = [Ppk_top; -Ppk_bot];

% sort by X-value
[Xpk, iXpk] = sort(Xpk);
Ypk = Ypk(iXpk);
Wpk = Wpk(iXpk);
Ppk = Ppk(iXpk);
end