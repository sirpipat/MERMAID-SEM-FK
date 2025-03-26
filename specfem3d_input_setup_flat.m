function ddir = specfem3d_input_setup_flat(ddir, bottom, depth, freq, baz, theta, fs, gpu_mode, stf, origin_time, nsteps)
% ddir = SPECFEM3D_INPUT_SETUP_FLAT(ddir, bottom, depth, freq, baz, theta, fs, gpu_mode, stf, origin_time, nsteps)
%
% Generates Par_file, cmtsolution file, station file, Mesh_Par_file, and 
% interface file for a fluid-solid simulation.
%
% INPUT:
% ddir          directory for the input files
% bottom        depth of the ocean bottom
% depth         depth of the hydrophone float
% freq          characteristic frequency for built-in source-time function
%               in FK code
% baz           backazimuth of the plane wave entering the box
% theta         incidence angle in degrees
% fs            sampling frequency [Default: 10]
% gpu_mode      whether to enable GPU MODE [Default: false]
% stf           source-time function full filename --or--
%               source-time function as [t x]
%               [Default: [], using built-in Gaussian function in FK code]
% origin_time   origin time of the seimogram (not working properly now)
%               [Default: 0]
% nsteps        number of time steps of the simulations [Default: 10000]
%
% OUTPUT:
% ddir          directory for the input files
%
% Last modified by sirawich-at-princeton.edu, 03/26/2025

defval('bottom', 4128)
defval('depth', 1518)
defval('freq', 1)
defval('fs', 10)
defval('baz', 266.235565)
defval('theta', 10.400274)
defval('gpu_mode', false)
defval('stf', nan)
defval('nsteps', 10000)
defval('origin_time', 0)

% create the DATA/ folder for the files generated by this function
system(sprintf('mkdir %s', ddir));
system(sprintf('mkdir %s/DATA/', ddir));
system(sprintf('mkdir %s/DATA/meshfem3D_files/', ddir));

% xspecfem3d parameters
params = makeparams3d;
fparams = fullfile(ddir, 'DATA', 'Par_file');
params.NSTEP = nsteps;
params.GPU_MODE = gpu_mode;
writeparfile3d(params, fparams);

% write source-time function file
if ischar(stf)
    stf_fname = removepath(stf);
    system(sprintf('cp %s %s/%s', stf, ddir, stf_fname));
elseif all(~isnan(stf), 'all')
    stf_fname = 'stf_file.txt';
    writetimeseries(stf(:,1), stf(:,2), fullfile(ddir, stf_fname));
end

% fkmodel parameters
fkmodel = makefkmodel;
ffkmodel = fullfile(ddir, 'DATA', 'FKMODEL');
fkmodel.nlayers = 2;
layer1 = struct('rho', 1020, 'vp', 1500, 'vs', 0, 'ztop', 0);
layer2 = struct('rho', 2500, 'vp', 3400, 'vs', 1963, 'ztop', -bottom);
fkmodel.layers = {layer1; layer2};
fkmodel.baz = baz;
fkmodel.theta = theta;
fkmodel.fmax = freq;
fkmodel.fs = fs;
fkmodel.twindow = params.NSTEP * params.DT;
fkmodel.origin_time = origin_time;
fkmodel.origin_wavefront = [-140000 -14000 -10000];
if any(isnan(stf), 'all')
    fkmodel.stf_type = 1;
    fkmodel.stf_file = "n/a";
else
    fkmodel.stf_type = 4;
    fkmodel.stf_file = stf_fname;
end
writefkmodel(fkmodel, ffkmodel);

% source parameters (unused)
cmt = makecmtsolution3d;
fcmt = fullfile(ddir, 'DATA', 'CMTSOLUTION');
writecmtsolution3d(cmt, fcmt);

% station
stations = struct('name', {{'OBS01', 'P0009'}}, ...
    'network', {{'AA', 'MH'}}, 'lat', [0 0], ...
    'lon', [0 0], 'elev', [-bottom -bottom], 'z', [-bottom -depth]);
fstations = fullfile(ddir, 'DATA', 'STATIONS');
writestations3d(stations, fstations);

% xmeshfem3d parameters
meshparams = makemeshparams3d;
fmesh = fullfile(ddir, 'DATA', 'meshfem3D_files', 'Mesh_Par_file');
% water
material1 = struct('rho', 1020, 'vp', 1500, 'vs', 0, 'Q_Kappa', 9999, ...
    'Q_mu', 9999, 'anisotropy_flag', 0, 'domain_id', 1);
% crust
material2 = struct('rho', 2500, 'vp', 3400, 'vs', 1963, ...
    'Q_Kappa', 9999, 'Q_mu', 9999, 'anisotropy_flag', 0, 'domain_id', 2);
materials = {material1; material2};
meshparams.NMATERIALS = 2;
meshparams.MATERIALS = materials;
% regions
nz_all = 20;
nz_top = round(nz_all*bottom/(meshparams.DEPTH_BLOCK_KM * 1000));
nz_bottom = nz_all - nz_top;
region1 = struct(...
    'NEX_XI_BEGIN'      , 1             , ...
    'NEX_XI_END'        , 56            , ...
    'NEX_ETA_BEGIN'     , 1             , ...
    'NEX_ETA_END'       , 56            , ...
    'NZ_BEGIN'          , nz_bottom+1   , ...
    'NZ_END'            , 20            , ...
    'material_id'       , 1               ...
);

region2 = struct(...
    'NEX_XI_BEGIN'      , 1             , ...
    'NEX_XI_END'        , 56            , ...
    'NEX_ETA_BEGIN'     , 1             , ...
    'NEX_ETA_END'       , 56            , ...
    'NZ_BEGIN'          , 1             , ...
    'NZ_END'            , nz_bottom     , ...
    'material_id'       , 2               ...
);
meshparams.NREGIONS = 2;
meshparams.REGIONS = {region1; region2};
writemeshparfile3d(meshparams, fmesh);

% interface file parameters
finterf = fullfile(ddir, 'DATA', 'meshfem3D_files', ...
    meshparams.INTERFACES_FILE);
% interfaces
itf1 = struct(...
    'SUPPRESS_UTM_PROJECTION'   , true  , ...
    'NXI'                       , 2     , ...
    'NETA'                      , 2     , ...
    'LON_MIN'   , meshparams.LONGITUDE_MIN, ...
    'LAT_MIN'   , meshparams.LATITUDE_MIN, ...
    'SPACING_XI', meshparams.LONGITUDE_MAX - meshparams.LONGITUDE_MIN, ...
    'SPACING_ETA', meshparams.LATITUDE_MAX - meshparams.LATITUDE_MIN, ...
    'FILE', fullfile(ddir, 'DATA', 'meshfem3D_files', 'interf_ocean.dat'), ...
    'Z', ones(2,2) * -bottom ...
);
itf2 = struct(...
    'SUPPRESS_UTM_PROJECTION'   , true  , ...
    'NXI'                       , 2     , ...
    'NETA'                      , 2     , ...
    'LON_MIN'   , meshparams.LONGITUDE_MIN, ...
    'LAT_MIN'   , meshparams.LATITUDE_MIN, ...
    'SPACING_XI', meshparams.LONGITUDE_MAX - meshparams.LONGITUDE_MIN, ...
    'SPACING_ETA', meshparams.LATITUDE_MAX - meshparams.LATITUDE_MIN, ...
    'FILE', fullfile(ddir, 'DATA', 'meshfem3D_files', 'interf_top.dat'), ...
    'Z', zeros(2,2) ...
);
itfs = {itf1; itf2};
layers = [nz_bottom; nz_top];
writeinterfacefiles3d(itfs, layers, finterf);
end