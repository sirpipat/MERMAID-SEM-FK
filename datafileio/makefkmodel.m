function fkmodel = makefkmodel
% fkmodel = MAKEFKMODEL
%
% Makes a sensible fkmodel for FK-SPECFEM3D_Cartesian.
%
% OUTPUT:
% fkmodel           model for embeded FK modelling
%
% SEE ALSO:
% LOADFKMODEL, WRITEFKMODEL
%
% Last modified by sirawich-at-princeton.edu, 02/26/2025

% model description
fkmodel.nlayers = 3;
rho = [1020 2677.5 3300];
vp = [1500 5250 8031.15];
vs = [0 2835 4576.15];
ztop = [0 -5000 -7500];
for ii = 1:fkmodel.nlayers
    fkmodel.layers{ii,1} = struct('rho', rho(ii), 'vp', vp(ii), ...
        'vs', vs(ii), 'ztop', ztop(ii));
end

% incident wave
fkmodel.wave = 'p';

% angle of incoming wave
fkmodel.baz = 0;
fkmodel.theta = 10;

% frequency max
fkmodel.fmax = 1;

% sampling frequency
fkmodel.fs = 10;

% time window
fkmodel.twindow = 50;

% amplitude
fkmodel.amplitude = 1;

% optional
fkmodel.origin_wavefront = nan(1,3);
fkmodel.origin_time = nan;

fkmodel.stf_type = nan;
fkmodel.stf_file = "";
end