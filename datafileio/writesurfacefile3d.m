function writesurfacefile3d(fname, elev)
% WRITESURFACEFILE3D(fname, elev)
%
% Writes a surface file containing elevations at grid points on the
% interface used by an interfacefile.
%
% INPUT:
% fname         name of the surface file
% elev          elevation grid of a surface
%
% SEE ALSO:
% LOADSURFACEFILE3D, WRITEINTERFACEFILE3D
%
% Last modified by sirawich-at-princeton.edu, 06/12/2025

fid = fopen(fname, 'w');
if size(elev, 2) > 1
    elev = reshape(elev', numel(elev), 1);
end
fprintf(fid, '%g\n', elev);
fclose(fid);
end