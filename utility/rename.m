function rename(varargin)
% RENAME - Rename a variable in the workspace
% Usage: 
%   rename('oldvar', 'newvar')  % function syntax
%   rename oldvar newvar        % command syntax
%
% Last modified by sirawich-at-princeton.edu, 07/15/2025

% Handle different input formats
if nargin == 1
    % Command syntax: single string with space-separated arguments
    args = strsplit(varargin{1});
    if length(args) ~= 2
        error('Command syntax requires exactly two arguments: rename oldvar newvar');
    end
    oldvar = args{1};
    newvar = args{2};
elseif nargin == 2
    % Function syntax: two separate arguments
    oldvar = varargin{1};
    newvar = varargin{2};
else
    error('Usage: rename(''oldvar'', ''newvar'') or rename oldvar newvar');
end

% Check if old variable exists in caller workspace
if ~evalin('caller', sprintf('exist(''%s'', ''var'')', oldvar))
    error('Variable ''%s'' does not exist in the workspace', oldvar);
end

% Check if new variable name is valid
if ~isvarname(newvar)
    error('''%s'' is not a valid variable name', newvar);
end

% Get the value from the old variable
value = evalin('caller', oldvar);

% Assign to new variable in caller workspace
assignin('caller', newvar, value);

% Clear the old variable from caller workspace
evalin('caller', sprintf('clear %s', oldvar));

fprintf('Variable ''%s'' renamed to ''%s''\n', oldvar, newvar);
end