function writecmtsolution3d(cmt, fname)
% WRITECMTSOLUTION3D(cmt, fname)
%
% Writes a CMTSOLUTION file of SPECFEM3D_Cartesian.
%
% INPUT:
% cmt           CMT solutions, array of struct(s) with following fields
%               - HDR
%                   - TYPE
%                   - YEAR
%                   - MONTH
%                   - DAY
%                   - HOUR
%                   - MINUTE
%                   - SECOND
%                   - LAT
%                   - LON
%                   - DEPTH
%                   - mb
%                   - Ms
%                   - NAME
%               - NAME              event name
%               - TSHIFT            time shift
%               - HALFDUR           half duration
%               - LAT               latitude
%               - LON               longitude
%               - DEPTH             depth
%               - Mrr
%               - Mtt
%               - Mpp
%               - Mrt
%               - Mrp
%               - Mtp
%               [default: makecmtsolution3d()]
% fname         name of the Par_file
%
% The CMT solution format can be found at
% https://specfem3d.readthedocs.io/en/latest/05_running_the_solver/
%
% SEE ALSO:
% LOADCMTSOLUTION3D, MAKECMTSOLUTION3D
%
% Last modified by Sirawich Pipatprathanporn, 09/30/2024

defval('fname', [])
defval('cmt', makecmtsolution3d)

% open the file
if isempty(fname)
    % standard output aka console output
    fid = 1;
else
    fid = fopen(fname, 'w');
end

for ii = 1:length(cmt)
    hdr = cmt{ii}.HDR;
    fprintf(fid, '%4s%d %d %d %d %d %f %f %f %f %f %f %s\n', ...
        hdr.TYPE, hdr.YEAR, hdr.MONTH, hdr.DAY, hdr.HOUR, hdr.MINUTE, ...
        hdr.SECOND, hdr.LAT, hdr.LON, hdr.DEPTH, hdr.mb, hdr.Ms, hdr.NAME);
    fprintf(fid, '%s:\t%s\n', 'event name', cmt{ii}.NAME);
    fprintf(fid, '%s:\t%f\n', 'time shift', cmt{ii}.TSHIFT);
    fprintf(fid, '%s:\t%f\n', 'half duration', cmt{ii}.HALFDUR);
    fprintf(fid, '%s:\t%f\n', 'latitude', cmt{ii}.LAT);
    fprintf(fid, '%s:\t%f\n', 'longitude', cmt{ii}.LON);
    fprintf(fid, '%s:\t%f\n', 'depth', cmt{ii}.DEPTH);
    fprintf(fid, '%s:\t%e\n', 'Mrr', cmt{ii}.Mrr);
    fprintf(fid, '%s:\t%e\n', 'Mtt', cmt{ii}.Mtt);
    fprintf(fid, '%s:\t%e\n', 'Mpp', cmt{ii}.Mpp);
    fprintf(fid, '%s:\t%e\n', 'Mrt', cmt{ii}.Mrt);
    fprintf(fid, '%s:\t%e\n', 'Mrp', cmt{ii}.Mrp);
    fprintf(fid, '%s:\t%e\n\n', 'Mtp', cmt{ii}.Mtp);
end

% close the file
if fid >= 3
    fclose(fid);
end
end