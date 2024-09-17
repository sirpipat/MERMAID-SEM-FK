function writefkmodel(fkmodel, fname)
% WRITEFKMODEL(fkmodel, fname)
%
% Writes a FKMODEL file of SPECFEM3D_Cartesian.
%
% INPUT:
% fkmodel       struct containing elastic properties of the media in
%               FK-domain and the properties of the incoming wave
% fname         name of the FKMODEL file
%
% SEE ALSO:
% LOADFKMODEL, MAKEFKMODEL
%
% Last modified by sirawich-at-princeton.edu, 09/17/2024

defval('fkmodel', makefkmodel)
defval('fname', [])

if isempty(fname)
    % standard output aka console output
    fid = 1;
else
    fid = fopen(fname, 'w');
end

fprintf(fid, '#\n');
fprintf(fid, '#  input file for embedded FK modelling\n');
fprintf(fid, '#\n');
fprintf(fid, '#  for each layer we give :\n');
fprintf(fid, '#  LAYER ilayer rho vp vs ztop\n');
fprintf(fid, '#  the last layer is the homogeneous half space\n');
fprintf(fid, '#\n');
fprintf(fid, '#\n');

fprintf(fid, '# model description  ---------------------\n');
fprintf(fid, 'NLAYER             %d\n', fkmodel.nlayers);
for ii = 1:fkmodel.nlayers
    fprintf(fid, 'LAYER %d %.2f %.2f %.2f %.2f\n', ii, ...
        fkmodel.layers{ii}.rho, fkmodel.layers{ii}.vp, ...
        fkmodel.layers{ii}.vs, fkmodel.layers{ii}.ztop);
end

fprintf(fid, '#----------------------------------------\n');
fprintf(fid, '# incident wave p or sv\n');
fprintf(fid, 'INCIDENT_WAVE     %s\n', fkmodel.wave);

fprintf(fid, '#----------------------------------------\n');
fprintf(fid, '# anlges of incomming wave\n');
fprintf(fid, 'BACK_AZIMUTH           %.2f\n', fkmodel.baz);
fprintf(fid, 'TAKE_OFF               %.2f\n', fkmodel.theta);

fprintf(fid, '#----------------------------------------\n');
fprintf(fid, 'FREQUENCY_MAX      %.2f\n', fkmodel.fmax);

fprintf(fid, '#----------------------------------------\n');
fprintf(fid, 'TIME_WINDOW       %.2f\n', fkmodel.twindow);

fprintf(fid, '#----------------------------------------\n');
fprintf(fid, '# optional\n');
if isnan(fkmodel.origin_wavefront)
    fprintf(fid, '#ORIGIN_WAVEFRONT  0 0 -40000\n');
else
    fprintf(fid, 'ORIGIN_WAVEFRONT  %.2f %.2f %.2f\n', ...
        fkmodel.origin_wavefront(1), fkmodel.origin_wavefront(2), ...
        fkmodel.origin_wavefront(3));
end
if isnan(fkmodel.origin_time)
    fprintf(fid, '#ORIGIN_TIME 0.');
else
    fprintf(fid, '#ORIGIN_TIME %.2f', fkmodel.origin_time);
end

if fid >= 3
    fclose(fid);
end
end