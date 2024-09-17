function fkmodel = loadfkmodel(fname)
% fkmodel = LOADFKMODEL(fname)
%
% Loads a FKMODEL file of SPECFEM3D_Cartesian.
%
% INPUT:
% fname         name of the FKMODEL file
%
% OUTPUT:
% fkmodel       struct containing elastic properties of the media in
%               FK-domain and the properties of the incoming wave
%
% SEE ALSO:
% MAKEFKMODEL, WRITEFKMODEL
%
% Last modified by sirawich-at-princeton.edu, 09/16/2024

fid = fopen(fname, 'r');

line = fgetl(fid);
% skip the comments
while strcmp(line(1), '#')
    line = fgetl(fid);
end
fkmodel.nlayer = sscanf(line(7:end), '%d', 1);
fkmodel.layers = cell(fkmodel.nlayer, 1);
for ii = 1:fkmodel.nlayer
    line = fgetl(fid);
    numbers = sscanf(line(6:end), '%f');
    fkmodel.layers{ii}.rho = numbers(2);
    fkmodel.layers{ii}.vp = numbers(3);
    fkmodel.layers{ii}.vs = numbers(4);
    fkmodel.layers{ii}.ztop = numbers(5);
end

line = fgetl(fid);
% skip the comments
while strcmp(line(1), '#')
    line = fgetl(fid);
end
fkmodel.wave = cindeks(split(line), 2);

line = fgetl(fid);
% skip the comments
while strcmp(line(1), '#')
    line = fgetl(fid);
end
fkmodel.baz = str2double(cindeks(split(line), 2));
line = fgetl(fid);
fkmodel.theta = str2double(cindeks(split(line), 2));

line = fgetl(fid);
% skip the comments
while strcmp(line(1), '#')
    line = fgetl(fid);
end
fkmodel.fmax = str2double(cindeks(split(line), 2));

line = fgetl(fid);
% skip the comments
while strcmp(line(1), '#')
    line = fgetl(fid);
end
fkmodel.twindow = str2double(cindeks(split(line), 2));

% optional
fkmodel.origin_wavefront = [nan nan nan];
fkmodel.origin_time = nan;
for ii = 1:2
    line = fgetl(fid);
    % skip the comments
    while strcmp(line(1), '#')
        line = fgetl(fid);
        if isempty(line) || ~ischar(line)
            return
        end
    end
    if isempty(line) || ~ischar(line)
        return
    end
    words = split(strip(line));
    if strcmp(words{1}, 'ORIGIN_WAVEFRONT')
        fkmodel.origin_wavefront = [str2double(words{2}) ...
            str2double(words{3}) str2double(words{4})];
    elseif strcmp(words{1}, 'ORIGIN_TIME')
        fkmodel.origin_time = str2double(words{2});
    end
end

fclose(fid);
end