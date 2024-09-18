function writemeshparfile3d(meshparams, fname)
% WRITEMESHPARFILE3D(meshparams, fname)
%
% Writes parameters to a Mesh_Par_file of SPECFEM3D_Cartesian.
%
% DISCLAIMER: This is not the official way to read/write Mesh_Par_file. I 
% just go through comments and parameters in an instant of Mesh_Par_file 
% and read/write accordingly.
%
% INPUT:
% meshparams        parameters
% fname             name of the Mesh_Par_file
%
% SEE ALSO:
% LOADMESHPARFILE3D, MAKEMESHPARAMS3D
%
% Last modified by Sirawich Pipatprathanporn, 09/18/2024

defval('params', makemeshparams3d())
defval('fname', []);

%% List of Helper functions (see the end of this file)
% writetitle(fid, title)                    -- write a title section
% writeblank(fid)                           -- write a blank line
% writecomment(fid, comment)                -- write a line/block comment
% writebool(fid, name, value, comment)      -- write a boolean variable
% writeint(fid, name, value, comment)       -- write an integer variable
% writefloat(fid, name, value, comment, numsigfig)-- write a float variable
% writestring(fid, name, value, comment)    -- write a string variable

%% open the file
if isempty(fname)
    % standard output aka console output
    fid = 1;
else
    fid = fopen(fname, 'w');
end

%% Meshing input parameters
writetitle(fid, ' Meshing input parameters');
writeblank(fid);

writecomment(fid, ['# coordinates of mesh block in latitude/longitude ' ...
    'and depth in km']);
writefloat(fid, 'LATITUDE_MIN', meshparams.LATITUDE_MIN);
writefloat(fid, 'LATITUDE_MAX', meshparams.LATITUDE_MAX);
writefloat(fid, 'LONGITUDE_MIN', meshparams.LONGITUDE_MIN);
writefloat(fid, 'LONGITUDE_MAX', meshparams.LONGITUDE_MAX);
writefloat(fid, 'DEPTH_BLOCK_KM', meshparams.DEPTH_BLOCK_KM);
writeint(fid, 'UTM_PROJECTION_ZONE', meshparams.UTM_PROJECTION_ZONE);
writebool(fid, 'SUPPRESS_UTM_PROJECTION', ...
    meshparams.SUPPRESS_UTM_PROJECTION);
writeblank(fid);

writecomment(fid, ['# file that contains the interfaces of the model ' ...
    '/ mesh']);
writestring(fid, 'INTERFACES_FILE', meshparams.INTERFACES_FILE);
writeblank(fid);

writecomment(fid, '# file that contains the cavity');
writestring(fid, 'CAVITY_FILE', meshparams.CAVITY_FILE);
writeblank(fid);

writecomment(fid, ['# number of elements at the surface along edges ' ...
    'of the mesh at the surface']);
writecomment(fid, ['# (must be 8 * multiple of NPROC below if mesh is ' ...
    'not regular and contains mesh doublings)']);
writecomment(fid, ['# (must be multiple of NPROC below if mesh is ' ...
    'regular)']);
writeint(fid, 'NEX_XI', meshparams.NEX_XI);
writeint(fid, 'NEX_ETA', meshparams.NEX_ETA);
writeblank(fid);

writecomment(fid, ['# number of MPI processors along xi and eta ' ...
    '(can be different)']);
writeint(fid, 'NPROC_XI', meshparams.NPROC_XI);
writeint(fid, 'NPROC_ETA', meshparams.NPROC_ETA);
writeblank(fid);

%% Doubling layers
writetitle(fid, 'Doubling layers');
writeblank(fid);

writecomment(fid, '# Regular/irregular mesh');
writebool(fid, 'USE_REGULAR_MESH', meshparams.USE_REGULAR_MESH);
writecomment(fid, ['# Only for irregular meshes, number of doubling ' ...
    'layers and their position']);
writeint(fid, 'NDOUBLINGS', meshparams.NDOUBLINGS);
writecomment(fid, ['# NZ_DOUBLING_1 is the parameter to set up if ' ...
    'there is only one doubling layer']);
writecomment(fid, ['# (more doubling entries can be added if needed ' ...
    'to match NDOUBLINGS value)']);
for ii = 1:max(2, meshparams.NDOUBLINGS)
    varname = sprintf('NZ_DOUBLING_%d', ii);
    writeint(fid, varname, meshparams.NZ_DOUBLINGS(ii));
end
writeblank(fid);

%% Visualization
writetitle(fid, 'Visualization');
writeblank(fid);

writecomment(fid, ['# create mesh files for visualisation or further ' ...
    'checking']);
writebool(fid, 'CREATE_ABAQUS_FILES', meshparams.CREATE_ABAQUS_FILES);
writebool(fid, 'CREATE_DX_FILES', meshparams.CREATE_DX_FILES);
writebool(fid, 'CREATE_VTK_FILES', meshparams.CREATE_VTK_FILES);
writeblank(fid);

writecomment(fid, ['# stores mesh files as cubit-exported files into ' ...
    'directory MESH/ (for single process run)']);
writebool(fid, 'SAVE_MESH_AS_CUBIT', meshparams.SAVE_MESH_AS_CUBIT);
writeblank(fid);

writecomment(fid, '# path to store the databases files');
writestring(fid, 'LOCAL_PATH', meshparams.LOCAL_PATH);
writeblank(fid);

%% CPML
writetitle(fid, 'CPML');
writeblank(fid);

writecomment(fid, '# CPML perfectly matched absorbing layers');
writefloat(fid, 'THICKNESS_OF_X_PML', meshparams.THICKNESS_OF_X_PML);
writefloat(fid, 'THICKNESS_OF_Y_PML', meshparams.THICKNESS_OF_Y_PML);
writefloat(fid, 'THICKNESS_OF_Z_PML', meshparams.THICKNESS_OF_Z_PML);
writeblank(fid);

%% Domain materials
writetitle(fid, 'Domain materials');
writeblank(fid);

writecomment(fid, '# number of materials');
writeint(fid, 'NMATERIALS', meshparams.NMATERIALS);
block_comment = {...
    '# define the different materials in the model as:', ...
    ['# #material_id  #rho  #vp  #vs  #Q_Kappa  #Q_mu  ' ...
    '#anisotropy_flag  #domain_id'], ...
    '#     Q_Kappa          : Q_Kappa attenuation quality factor', ...
    '#     Q_mu             : Q_mu attenuation quality factor', ...
    ['#     anisotropy_flag  : 0 = no anisotropy / 1,2,... check the ' ...
    'implementation in file aniso_model.f90'], ...
    '#     domain_id        : 1 = acoustic / 2 = elastic'};
writecomment(fid, block_comment);
for ii = 1:meshparams.NMATERIALS
    fprintf(fid, '%d  %.2f %.2f %.2f %.2f %.2f %d %d\n', ii, ...
        meshparams.MATERIALS{ii}.rho, meshparams.MATERIALS{ii}.vp, ...
        meshparams.MATERIALS{ii}.vp, meshparams.MATERIALS{ii}.Q_Kappa, ...
        meshparams.MATERIALS{ii}.Q_mu, ...
        meshparams.MATERIALS{ii}.anisotropy_flag, ...
        meshparams.MATERIALS{ii}.domain_id);
end
writeblank(fid);

%% Domain regions
writetitle(fid, 'Domain regions');
writeblank(fid);

writecomment(fid, '# number of regions');
writeint(fid, 'NREGIONS', meshparams.NREGIONS);
writecomment(fid, '# define the different regions of the model as :');
writecomment(fid, sprintf('%-16s%-16s%-16s%-16s%-16s%-16s%-16s', ...
    '#NEX_XI_BEGIN', '#NEX_XI_END', '#NEX_ETA_BEGIN', '#NEX_ETA_END', ...
    '#NZ_BEGIN', '#NZ_END', '#material_id'));
for ii = 1:meshparams.NREGIONS
    fprintf(fid, '%-16d%-16d%-16d%-16d%-16d%-16d%-16d\n', ...
        meshparams.REGIONS{ii}.NEX_XI_BEGIN, ...
        meshparams.REGIONS{ii}.NEX_XI_END, ...
        meshparams.REGIONS{ii}.NEX_ETA_BEGIN, ...
        meshparams.REGIONS{ii}.NEX_ETA_END, ...
        meshparams.REGIONS{ii}.NZ_BEGIN, ...
        meshparams.REGIONS{ii}.NZ_END, ...
        meshparams.REGIONS{ii}.material_id);
end

%% close the file
if fid >= 3
    fclose(fid);
end
end

%% Helper functions
function writetitle(fid, title)
fprintf(fid, ['#------------------------------------------------------' ...
    '-----\n']);
fprintf(fid, '#\n');
fprintf(fid, '# %s\n', title);
fprintf(fid, '#\n');
fprintf(fid, ['#------------------------------------------------------' ...
    '-----\n']);
end

function writeblank(fid)
fprintf(fid, '\n');
end

function writecomment(fid, comment)
defval('comment', [])
if isempty(comment)
    writeblank(fid);
% write a block comment represented by a cell array of strings
elseif iscell(comment)
    for ii = 1:length(comment)
        writecomment(fid, comment{ii});
    end
% write an individual comment (string)
else
    if ~strcmp(comment(1), '#')
        fprintf(fid, '# %s\n', comment);
    else
        fprintf(fid, '%s\n', comment);
    end
end
end

function writebool(fid, name, value, comment)
defval('comment', [])
if value == 0
    var_string = '.false.';
else
    var_string = '.true.';
end
writestring(fid, name, var_string, comment);
end


function writefloat(fid, name, value, comment, numsigfig)
defval('comment', [])
defval('numsigfig', 3)
if numsigfig < 1
    numsigfig = 1;
end
if value == 0
    var_string = '0';
elseif and(abs(value) >= 1, abs(value) < 100)
    var_string = sprintf(sprintf('%%.%df', round(numsigfig)), value);
else
    var_string = sprintf(sprintf('%%.%de', round(numsigfig)), value);
    var_string = replace(var_string, 'e', 'd');
end
writestring(fid, name, var_string, comment);
end

function writeint(fid, name, value, comment)
defval('comment', [])
if isempty(comment)
    fprintf(fid, '%-31s = %d\n', name, value);
else
    if ~strcmp(comment(1), '#')
        fprintf(fid, '%-31s = %-14d # %s\n', name, value, comment);
    else
        fprintf(fid, '%-31s = %-14d %s\n', name, value, comment);
    end
end
end

function writestring(fid, name, value, comment)
defval('comment', [])
if isempty(comment)
    fprintf(fid, '%-31s = %s\n', name, value);
else
    if ~strcmp(comment(1), '#')
        fprintf(fid, '%-31s = %-14s # %s\n', name, value, comment);
    else
        fprintf(fid, '%-31s = %-14s %s\n', name, value, comment);
    end
end
end