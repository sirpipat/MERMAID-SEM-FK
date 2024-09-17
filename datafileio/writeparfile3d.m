function writeparfile3d(params, fname)
% WRITEPARFILE3D(params, fname)
%
% Writes parameters to a Par_file of SPECFEM3D_Cartesian.
%
% DISCLAIMER: This is not the official way to read/write Par_file. I just
% go through comments and parameters in an instant of Par_file and
% read/write accordingly.
%
% INPUT:
% params        parameters
% fname         name of the Par_file
%
% SEE ALSO:
% LOADPARFILE3D, MAKEPARAMS3D
%
% Last modified by Sirawich Pipatprathanporn, 09/17/2024

defval('params', makeparams3d())

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

%% Simulation input parameters
writetitle(fid, 'Simulation input parameters');
writeblank(fid);

writecomment(fid, '# forward or adjoint simulation');
writecomment(fid, '# 1 = forward, 2 = adjoint, 3 = both simultaneously');
writeint(fid, 'SIMULATION_TYPE', params.SIMULATION_TYPE);
writecomment(fid, ['# 0 = earthquake simulation,  ' ...
    '1/2/3 = three steps in noise simulation']);
writeint(fid, 'NOISE_TOMOGRAPHY', params.NOISE_TOMOGRAPHY);
writebool(fid, 'SAVE_FORWARD', params.SAVE_FORWARD);
writeblank(fid);

writecomment(fid, ['# solve a full FWI inverse problem from a single ' ...
    'calling program with no I/Os, storing everything in memory,']);
writecomment(fid, ['# or run a classical forward or adjoint problem ' ...
    'only and save the seismograms and/or sensitivity kernels to disk ' ...
    '(with costlier I/Os)']);
writebool(fid, 'INVERSE_FWI_FULL_PROBLEM', ...
    params.INVERSE_FWI_FULL_PROBLEM);
writeblank(fid);

writecomment(fid, '# UTM projection parameters');
writecomment(fid, ['# Use a negative zone number for the Southern ' ...
    'hemisphere:']);
writecomment(fid, ['# The Northern hemisphere corresponds to zones ' ...
    '+1 to +60,']);
writecomment(fid, ['# The Southern hemisphere corresponds to zones ' ...
    '-1 to -60.']);
writeint(fid, 'UTM_PROJECTION_ZONE', params.UTM_PROJECTION_ZONE);
writebool(fid, 'SUPPRESS_UTM_PROJECTION', params.SUPPRESS_UTM_PROJECTION);
writeblank(fid);

writecomment(fid, '# number of MPI processors');
writeint(fid, 'NPROC', params.NPROC);
writeblank(fid);

writecomment(fid, '# time step parameters');
writeint(fid, 'NSTEP', params.NSTEP);
writefloat(fid, 'DT', params.DT);
writeblank(fid);

writecomment(fid, '# set to true to use local-time stepping (LTS)');
writebool(fid, 'LTS_MODE', params.LTS_MODE);
writeblank(fid);

writecomment(fid, '# Partitioning algorithm for decompose_mesh');
writecomment(fid, ['# choose partitioner: 1==SCOTCH (default), ' ...
    '2==METIS, 3==PATOH, 4==ROWS_PART']);
writeint(fid, 'PARTITIONING_TYPE', params.PARTITIONING_TYPE);
writeblank(fid);

%% LDDRK time scheme
writetitle(fid, 'LDDRK time scheme');
writebool(fid, 'USE_LDDRK', params.USE_LDDRK);
writebool(fid, 'INCREASE_CFL_FOR_LDDRK', params.INCREASE_CFL_FOR_LDDRK);
writefloat(fid, 'RATIO_BY_WHICH_TO_INCREASE_IT', ...
    params.RATIO_BY_WHICH_TO_INCREASE_IT);
writeblank(fid);

%% Mesh
writetitle(fid, 'Mesh');
writeblank(fid);

block_comment1 = {...
    '# Number of nodes for 2D and 3D shape functions for hexahedra.', ...
    '# We use either 8-node mesh elements (bricks) or 27-node elements.', ...
    ['# If you use our internal mesher, the only option is 8-node ' ...
    'bricks (27-node elements are not supported).']};
writecomment(fid, block_comment1);
writeint(fid, 'NGNOD', params.NGNOD);
writeblank(fid);

block_comment2 = {'# models:', ...
    '# available options are:', ...
    '#   default (model parameters described by mesh properties)', ...
    '# 1D models available are:', ...
    '#   1d_prem,1d_socal,1d_cascadia', ...
    '# 3D models available are:', ...
    '#   aniso,external,gll,salton_trough,tomo,SEP,coupled,...'};
writecomment(fid, block_comment2);
writestring(fid, 'MODEL', params.MODEL);
writeblank(fid);

writecomment(fid, '# path for external tomographic models files');
writestring(fid, 'TOMOGRAPHY_PATH', params.TOMOGRAPHY_PATH);
writecomment(fid, '# if you are using a SEP model (oil-industry format)');
writestring(fid, 'SEP_MODEL_DIRECTORY', params.SEP_MODEL_DIRECTORY);
writeblank(fid);

writecomment(fid, ['#-------------------------------------------------' ...
    '----------']);
writeblank(fid);

writecomment(fid, '# parameters describing the model');
writebool(fid, 'APPROXIMATE_OCEAN_LOAD', params.APPROXIMATE_OCEAN_LOAD, ...
    []);
writebool(fid, 'TOPOGRAPHY', params.TOPOGRAPHY);
writebool(fid, 'ATTENUATION', params.ATTENUATION);
writebool(fid, 'ANISOTROPY', params.ANISOTROPY);
writebool(fid, 'GRAVITY', params.GRAVITY);
writeblank(fid);

writecomment(fid, ['# in case of attenuation, reference frequency in ' ...
    'Hz at which the velocity values in the velocity model are given ' ...
    '(unused otherwise)']);
writefloat(fid, 'ATTENUATION_f0_REFERENCE', ...
    params.ATTENUATION_f0_REFERENCE);
writeblank(fid);

writecomment(fid, ['# attenuation period range over which we try to ' ...
    'mimic a constant Q factor']);
writefloat(fid, 'MIN_ATTENUATION_PERIOD', ...
    params.MIN_ATTENUATION_PERIOD, [], 9);
writefloat(fid, 'MAX_ATTENUATION_PERIOD', ...
    params.MAX_ATTENUATION_PERIOD, [], 9);
writecomment(fid, ['# ignore this range and ask the code to compute ' ...
    'it automatically instead based on the estimated resolution of ' ...
    'the mesh (use this unless you know what you are doing)']);
writebool(fid, 'COMPUTE_FREQ_BAND_AUTOMATIC', ...
    params.COMPUTE_FREQ_BAND_AUTOMATIC);
writeblank(fid);

writecomment(fid, ['# Olsen''s constant for Q_mu = constant * V_s ' ...
    'attenuation rule']);
writebool(fid, 'USE_OLSEN_ATTENUATION', params.USE_OLSEN_ATTENUATION);
writefloat(fid, 'OLSEN_ATTENUATION_RATIO', params.OLSEN_ATTENUATION_RATIO);
writeblank(fid);

%% Absorbing boundary conditions
writetitle(fid, 'Absorbing boundary conditions');
writeblank(fid);

block_comment3 = {...
    '# C-PML boundary conditions for a regional simulation', ...
    ['# (if set to .false., and STACEY_ABSORBING_CONDITIONS is also ' ...
    'set to .false., you get a free surface instead'], ...
    ['# in the case of elastic or viscoelastic mesh elements, and a ' ...
    'rigid surface in the case of acoustic (fluid) elements']};
writecomment(fid, block_comment3);
writebool(fid, 'PML_CONDITIONS', params.PML_CONDITIONS);
writeblank(fid);

writecomment(fid, '# C-PML top surface');
writebool(fid, 'PML_INSTEAD_OF_FREE_SURFACE', ...
    params.PML_INSTEAD_OF_FREE_SURFACE);
writeblank(fid);

writecomment(fid, '# C-PML dominant frequency');
writefloat(fid, 'f0_FOR_PML', params.f0_FOR_PML);
writeblank(fid);

writecomment(fid, ['# parameters used to rotate C-PML boundary ' ...
    'conditions by a given angle (not completed yet)']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% once C-PML boundary is implemented, "uncomment" these two lines below   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
writecomment(fid, '# ROTATE_PML_ACTIVATE           = .false.');
writecomment(fid, '# ROTATE_PML_ANGLE              = 0.');
writeblank(fid);

block_comment4 = {...
    '# absorbing boundary conditions for a regional simulation', ...
    ['# (if set to .false., and PML_CONDITIONS is also set to ' ...
    '.false., you get a free surface instead'], ...
    ['# in the case of elastic or viscoelastic mesh elements, and a ' ...
    'rigid surface in the case of acoustic (fluid) elements']};
writecomment(fid, block_comment4);
writebool(fid, 'STACEY_ABSORBING_CONDITIONS', ...
    params.STACEY_ABSORBING_CONDITIONS);
writeblank(fid);

writecomment(fid, ['# absorbing top surface (defined in mesh as ' ...
    '''free_surface_file'')']);
writebool(fid, 'STACEY_INSTEAD_OF_FREE_SURFACE', ...
    params.STACEY_INSTEAD_OF_FREE_SURFACE);
writeblank(fid);

block_comment5 = {...
    '# When STACEY_ABSORBING_CONDITIONS is set to .true. :', ...
    ['# absorbing conditions are defined in xmin, xmax, ymin, ymax ' ...
    'and zmin'], ...
    '# this option BOTTOM_FREE_SURFACE can be set to .true. to', ...
    '# make zmin free surface instead of absorbing condition'};
writecomment(fid, block_comment5);
writebool(fid, 'BOTTOM_FREE_SURFACE', params.BOTTOM_FREE_SURFACE);
writeblank(fid);

%% undoing attenuation and/or PMLs for sensitivity kernel calculations
writetitle(fid, ['undoing attenuation and/or PMLs for sensitivity ' ...
    'kernel calculations']);
writeblank(fid);

block_comment6 = {...
    ['# to undo attenuation and/or PMLs for sensitivity kernel ' ...
    'calculations or forward runs with SAVE_FORWARD'], ...
    ['# use the flag below. It performs undoing of attenuation and/or ' ...
    'of PMLs in an exact way for sensitivity kernel calculations'], ...
    ['# but requires disk space for temporary storage, and uses a ' ...
    'significant amount of memory used as buffers for temporary ' ...
    'storage.'], ...
    ['# When that option is on the second parameter indicates how ' ...
    'often the code dumps restart files to disk (if in doubt, use ' ...
    'something between 100 and 1000).']};
writecomment(fid, block_comment6);
writebool(fid, 'UNDO_ATTENUATION_AND_OR_PML', ...
    params.UNDO_ATTENUATION_AND_OR_PML);
writeint(fid, 'NT_DUMP_ATTENUATION', params.NT_DUMP_ATTENUATION);
writeblank(fid);

%% visualization
writetitle(fid, 'Visualization');
writeblank(fid);

writecomment(fid, '# save AVS or OpenDX movies');
writecomment(fid, '# MOVIE_TYPE = 1 to show the top surface');
writecomment(fid, ['# MOVIE_TYPE = 2 to show all the external faces ' ...
    'of the mesh']);
writebool(fid, 'CREATE_SHAKEMAP', params.CREATE_SHAKEMAP);
writebool(fid, 'MOVIE_SURFACE', params.MOVIE_SURFACE);
writeint(fid, 'MOVIE_TYPE', params.MOVIE_TYPE);
writebool(fid, 'MOVIE_VOLUME', params.MOVIE_VOLUME);
writebool(fid, 'SAVE_DISPLACEMENT', params.SAVE_DISPLACEMENT);
writebool(fid, 'MOVIE_VOLUME_STRESS', params.MOVIE_VOLUME_STRESS);
writebool(fid, 'USE_HIGHRES_FOR_MOVIES', params.USE_HIGHRES_FOR_MOVIES);
writeint(fid, 'NTSTEP_BETWEEN_FRAMES', params.NTSTEP_BETWEEN_FRAMES);
writefloat(fid, 'HDUR_MOVIE', params.HDUR_MOVIE);
writeblank(fid);

writecomment(fid, '# save AVS or OpenDX mesh files to check the mesh');
writebool(fid, 'SAVE_MESH_FILES', params.SAVE_MESH_FILES);
writeblank(fid);

writecomment(fid, '# path to store the local database file on each node');
writestring(fid, 'LOCAL_PATH', params.LOCAL_PATH);
writeblank(fid);

writecomment(fid, ['# interval at which we output time step info and ' ...
    'max of norm of displacement']);
writeint(fid, 'NTSTEP_BETWEEN_OUTPUT_INFO', ...
    params.NTSTEP_BETWEEN_OUTPUT_INFO);
writeblank(fid);

%% Sources
writetitle(fid, 'Sources');
writeblank(fid);

writecomment(fid, ['# sources and receivers Z coordinates given ' ...
    'directly (i.e. as their true position) instead of as their depth']);
writebool(fid, 'USE_SOURCES_RECEIVERS_Z', params.USE_SOURCES_RECEIVERS_Z);
writeblank(fid);

block_comment7 = {...
    ['# use a (tilted) FORCESOLUTION force point source (or several) ' ...
    'instead of a CMTSOLUTION moment-tensor source.'], ...
    ['# This can be useful e.g. for oil industry foothills ' ...
    'simulations or asteroid simulations'], ...
    ['# in which the source is a vertical force, normal force, tilted ' ...
    'force, impact etc.'], ...
    ['# If this flag is turned on, the FORCESOLUTION file must be ' ...
    'edited by giving:'], ...
    '# - the corresponding time-shift parameter,', ...
    ['# - the half duration (hdur, in s) for Gaussian/Step function, ' ...
    'dominant frequency (f0, in Hz) for Ricker,'], ...
    '# - the coordinates of the source,', ...
    ['# - the source time function type (0=Gaussian function, ' ...
    '1=Ricker wavelet, 2=Step function),'], ...
    '# - the magnitude of the force source,', ...
    ['# - the components of a (non necessarily unitary) direction ' ...
    'vector for the force source in the E/N/Z_UP basis.'], ...
    ['# The direction vector is made unitary internally in the code ' ...
    'and thus only its direction matters here;'], ...
    ['# its norm is ignored and the norm of the force used is the ' ...
    'factor force source times the source time function.']};
writecomment(fid, block_comment7);
writebool(fid, 'USE_FORCE_POINT_SOURCE', params.USE_FORCE_POINT_SOURCE);
writeblank(fid);

writecomment(fid, ['# set to true to use a Ricker source time ' ...
    'function instead of the source time functions set by default']);
writecomment(fid, ['# to represent a (tilted) FORCESOLUTION force ' ...
    'point source or a CMTSOLUTION moment-tensor source.']);
writebool(fid, 'USE_RICKER_TIME_FUNCTION', ...
    params.USE_RICKER_TIME_FUNCTION);
writeblank(fid);

block_comment8 = {'# use an external source time function', ...
    ['# you must add a file with your source time function and the ' ...
    'file name path'], ...
    ['# relative to working directory at the end of CMTSOLUTION or ' ...
    'FORCESOLUTION file'], ...
    '# (with multiple sources, one file per source is required).', ...
    ['# This file must have a single column containing the amplitudes ' ...
    'of the source at all time steps;'], ...
    ['# time step size used must be equal to DT as defined at the ' ...
    'beginning of this Par_file.'], ...
    ['# Be sure when this option is .false. to remove the name of stf ' ...
    'file in CMTSOLUTION or FORCESOLUTION']};
writecomment(fid, block_comment8);
writebool(fid, 'USE_EXTERNAL_SOURCE_FILE', ...
    params.USE_EXTERNAL_SOURCE_FILE);
writeblank(fid);

writecomment(fid, '# print source time function');
writebool(fid, 'PRINT_SOURCE_TIME_FUNCTION', ...
    params.PRINT_SOURCE_TIME_FUNCTION);
writeblank(fid);

block_comment9 = {'# source encoding', ...
    ['# (for acoustic simulations only for now) determines source ' ...
    'encoding factor +/-1 depending on sign of moment tensor'], ...
    ['# (see e.g. Krebs et al., 2009. Fast full-wavefield seismic ' ...
    'inversion using encoded sources, Geophysics, 74 (6), ' ...
    'WCC177-WCC188.)']};
writecomment(fid, block_comment9);
writebool(fid, 'USE_SOURCE_ENCODING', params.USE_SOURCE_ENCODING);
writeblank(fid);

%% Seismograms
writetitle(fid, 'Seismograms');
writeblank(fid);

writecomment(fid, '# interval in time steps for writing of seismograms');
writeint(fid, 'NTSTEP_BETWEEN_OUTPUT_SEISMOS', ...
    params.NTSTEP_BETWEEN_OUTPUT_SEISMOS);
writeblank(fid);

writecomment(fid, ['# set to n to reduce the sampling rate of output ' ...
    'seismograms by a factor of n']);
writecomment(fid, '# defaults to 1, which means no down-sampling');
writeint(fid, 'NTSTEP_BETWEEN_OUTPUT_SAMPLE', ...
    params.NTSTEP_BETWEEN_OUTPUT_SAMPLE);
writeblank(fid);

writecomment(fid, ['# decide if we save displacement, velocity, ' ...
    'acceleration and/or pressure in forward runs (they can be set to ' ...
    'true simultaneously)']);
writecomment(fid, ['# currently pressure seismograms are implemented ' ...
    'in acoustic (i.e. fluid) elements only']);
writebool(fid, 'SAVE_SEISMOGRAMS_DISPLACEMENT', ...
    params.SAVE_SEISMOGRAMS_DISPLACEMENT);
writebool(fid, 'SAVE_SEISMOGRAMS_VELOCITY', ...
    params.SAVE_SEISMOGRAMS_VELOCITY);
writebool(fid, 'SAVE_SEISMOGRAMS_ACCELERATION', ...
    params.SAVE_SEISMOGRAMS_ACCELERATION);
writebool(fid, 'SAVE_SEISMOGRAMS_PRESSURE', ...
    params.SAVE_SEISMOGRAMS_PRESSURE, ['# currently implemented in ' ...
    'acoustic (i.e. fluid) elements only']);
writeblank(fid);

writecomment(fid, '# option to save strain seismograms');
writecomment(fid, '# this option is useful for strain Green''s tensor');
writebool(fid, 'SAVE_SEISMOGRAMS_STRAIN', params.SAVE_SEISMOGRAMS_STRAIN);
writeblank(fid);

writecomment(fid, ['# save seismograms also when running the adjoint ' ...
    'runs for an inverse problem']);
writecomment(fid, ['# (usually they are unused and not very ' ...
    'meaningful, leave this off in almost all cases)']);
writebool(fid, 'SAVE_SEISMOGRAMS_IN_ADJOINT_RUN', ...
    params.SAVE_SEISMOGRAMS_IN_ADJOINT_RUN);
writeblank(fid);

writecomment(fid, ['# save seismograms in binary or ASCII format ' ...
    '(binary is smaller but may not be portable between machines)']);
writebool(fid, 'USE_BINARY_FOR_SEISMOGRAMS', ...
    params.USE_BINARY_FOR_SEISMOGRAMS);
writeblank(fid);

writecomment(fid, ['# output seismograms in Seismic Unix format ' ...
    '(binary with 240-byte-headers)']);
writebool(fid, 'SU_FORMAT', params.SU_FORMAT);
writeblank(fid);

writecomment(fid, '# output seismograms in ASDF (requires asdf-library)');
writebool(fid, 'ASDF_FORMAT', params.ASDF_FORMAT);
writeblank(fid);

writecomment(fid, ['# output seismograms in HDF5 (requires ' ...
    'hdf5-library and WRITE_SEISMOGRAMS_BY_MAIN)']);
writebool(fid, 'HDF5_FORMAT', params.HDF5_FORMAT);
writeblank(fid);

writecomment(fid, ['# decide if main process writes all the ' ...
    'seismograms or if all processes do it in parallel']);
writebool(fid, 'WRITE_SEISMOGRAMS_BY_MAIN', ...
    params.WRITE_SEISMOGRAMS_BY_MAIN);
writeblank(fid);

writecomment(fid, ['# save all seismograms in one large combined file ' ...
    'instead of one file per seismogram']);
writecomment(fid, ['# to avoid overloading shared non-local file ' ...
    'systems such as LUSTRE or GPFS for instance']);
writebool(fid, 'SAVE_ALL_SEISMOS_IN_ONE_FILE', ...
    params.SAVE_ALL_SEISMOS_IN_ONE_FILE);
writeblank(fid);


block_comment10 = {...
    ['# use a trick to increase accuracy of pressure ' ...
    'seismograms in fluid (acoustic) elements:'], ...
    ['# use the second derivative of the source for the source time ' ...
    'function instead of the source itself,'], ...
    ['# and then record -potential_acoustic() as pressure seismograms ' ...
    'instead of -potential_dot_dot_acoustic();'], ...
    ['# this is mathematically equivalent, but numerically ' ...
    'significantly more accurate because in the explicit'], ...
    ['# Newmark time scheme acceleration is accurate at zeroth order ' ...
    'while displacement is accurate at second order,'], ...
    ['# thus in fluid elements potential_dot_dot_acoustic() is ' ...
    'accurate at zeroth order while potential_acoustic()'], ...
    ['# is accurate at second order and thus contains significantly ' ...
    'less numerical noise.']};
writecomment(fid, block_comment10);
writebool(fid, 'USE_TRICK_FOR_BETTER_PRESSURE', ...
    params.USE_TRICK_FOR_BETTER_PRESSURE);
writeblank(fid);

%% Energy calculation
writetitle(fid, 'Energy calculation');
writeblank(fid);

writecomment(fid, ['# to plot energy curves, for instance to monitor ' ...
    'how CPML absorbing layers behave;']);
writecomment(fid, ['# should be turned OFF in most cases because a ' ...
    'bit expensive']);
writebool(fid, 'OUTPUT_ENERGY', params.OUTPUT_ENERGY);
writecomment(fid, ['# every how many time steps we compute energy ' ...
    '(which is a bit expensive to compute)']);
writeint(fid, 'NTSTEP_BETWEEN_OUTPUT_ENERGY', ...
    params.NTSTEP_BETWEEN_OUTPUT_ENERGY);
writeblank(fid);

%% Adjoint kernel outputs
writetitle(fid, 'Adjoint kernel outputs');
writeblank(fid);

writecomment(fid, '# interval in time steps for reading adjoint traces');
writecomment(fid, '# 0 = read the whole adjoint sources at start time');
writeint(fid, 'NTSTEP_BETWEEN_READ_ADJSRC', ...
    params.NTSTEP_BETWEEN_READ_ADJSRC);
writeblank(fid);

writecomment(fid, ['# read adjoint sources using ASDF (requires ' ...
    'asdf-library)']);
writebool(fid, 'READ_ADJSRC_ASDF', params.READ_ADJSRC_ASDF);
writeblank(fid);

block_comment11 = {...
    ['# this parameter must be set to .true. to compute anisotropic ' ...
    'kernels'], ...
    ['# in crust and mantle (related to the 21 Cij in geographical ' ...
    'coordinates)'], ...
    ['# default is .false. to compute isotropic kernels (related to ' ...
    'alpha and beta)']};
writecomment(fid, block_comment11);
writebool(fid, 'ANISOTROPIC_KL', params.ANISOTROPIC_KL);
writeblank(fid);

writecomment(fid, ['# compute transverse isotropic kernels ' ...
    '(alpha_v,alpha_h,beta_v,beta_h,eta,rho)']);
writecomment(fid, ['# rather than fully anisotropic kernels in case ' ...
    'ANISOTROPIC_KL is set to .true.']);
writebool(fid, 'SAVE_TRANSVERSE_KL', params.SAVE_TRANSVERSE_KL);
writeblank(fid);

writecomment(fid, ['# this parameter must be set to .true. to compute ' ...
    'anisotropic kernels for']);
writecomment(fid, ['# cost function using velocity observable rather ' ...
    'than displacement']);
writebool(fid, 'ANISOTROPIC_VELOCITY_KL', params.ANISOTROPIC_VELOCITY_KL);
writeblank(fid);

writecomment(fid, '# outputs approximate Hessian for preconditioning');
writebool(fid, 'APPROXIMATE_HESS_KL', params.APPROXIMATE_HESS_KL);
writeblank(fid);

writecomment(fid, '# save Moho mesh and compute Moho boundary kernels');
writebool(fid, 'SAVE_MOHO_MESH', params.SAVE_MOHO_MESH);
writeblank(fid);

%% Coupling with an injection technique (DSM, AxiSEM, or FK)
writetitle(fid, ['Coupling with an injection technique (DSM, AxiSEM, ' ...
    'or FK)']);
writebool(fid, 'COUPLE_WITH_INJECTION_TECHNIQUE', ...
    params.COUPLE_WITH_INJECTION_TECHNIQUE);
writeint(fid, 'INJECTION_TECHNIQUE_TYPE', ...
    params.INJECTION_TECHNIQUE_TYPE, '# 1 = DSM, 2 = AxiSEM, 3 = FK');
writebool(fid, 'MESH_A_CHUNK_OF_THE_EARTH', ...
    params.MESH_A_CHUNK_OF_THE_EARTH);
writestring(fid, 'TRACTION_PATH', params.TRACTION_PATH);
writestring(fid, 'FKMODEL_FILE', params.FKMODEL_FILE);
writebool(fid, 'RECIPROCITY_AND_KH_INTEGRAL', ...
    params.RECIPROCITY_AND_KH_INTEGRAL, '# does not work yet');
writeblank(fid);

%% Run modes
writetitle(fid, 'Run modes');
writeblank(fid);

block_comment12 = {'# Simultaneous runs', ...
    ['# added the ability to run several calculations (several ' ...
    'earthquakes)'], ...
    ['# in an embarrassingly-parallel fashion from within the same ' ...
    'run;'], ...
    ['# this can be useful when using a very large supercomputer to ' ...
    'compute'], ...
    ['# many earthquakes in a catalog, in which case it can be better ' ...
    'from'], ...
    ['# a batch job submission point of view to start fewer and much ' ...
    'larger jobs,'], ...
    '# each of them computing several earthquakes in parallel.', ...
    ['# To turn that option on, set parameter ' ...
    'NUMBER_OF_SIMULTANEOUS_RUNS to a value greater than 1.'], ...
    ['# To implement that, we create NUMBER_OF_SIMULTANEOUS_RUNS MPI ' ...
    'sub-communicators,'], ...
    ['# each of them being labeled "my_local_mpi_comm_world", and we ' ...
    'use them'], ...
    ['# in all the routines in "src/shared/parallel.f90", except in ' ...
    'MPI_ABORT() because in that case'], ...
    '# we need to kill the entire run.', ...
    ['# When that option is on, of course the number of processor ' ...
    'cores used to start'], ...
    ['# the code in the batch system must be a multiple of ' ...
    'NUMBER_OF_SIMULTANEOUS_RUNS,'], ...
    ['# all the individual runs must use the same number of processor ' ...
    'cores,'], ...
    '# which as usual is NPROC in the Par_file,', ...
    ['# and thus the total number of processor cores to request from ' ...
    'the batch system'], ...
    '# should be NUMBER_OF_SIMULTANEOUS_RUNS * NPROC.', ...
    ['# All the runs to perform must be placed in directories called ' ...
    'run0001, run0002, run0003 and so on'], ...
    '# (with exactly four digits).', ...
    '#', ...
    ['# Imagine you have 10 independent calculations to do, each of ' ...
    'them on 100 cores; you have three options:'], ...
    '#', ...
    '# 1/ submit 10 jobs to the batch system', ...
    '#', ...
    ['# 2/ submit a single job on 1000 cores to the batch, and in ' ...
    'that script create a sub-array of jobs to start 10 jobs,'], ...
    ['# each running on 100 cores (see e.g. ' ...
    'http://www.schedmd.com/slurmdocs/job_array.html )'], ...
    '#', ...
    ['# 3/ submit a single job on 1000 cores to the batch, start ' ...
    'SPECFEM3D on 1000 cores, create 10 sub-communicators,'], ...
    ['# cd into one of 10 subdirectories (called e.g. run0001, ' ...
    'run0002,... run0010) depending on the sub-communicator'], ...
    ['# your MPI rank belongs to, and run normally on 100 cores using ' ...
    'that sub-communicator.'], ...
    '#', ...
    '# The option below implements 3/.', ...
    '#'};
writecomment(fid, block_comment12);
writeint(fid, 'NUMBER_OF_SIMULTANEOUS_RUNS', ...
    params.NUMBER_OF_SIMULTANEOUS_RUNS);
writeblank(fid);

block_comment13 = {...
    ['# if we perform simultaneous runs in parallel, if only the ' ...
    'source and receivers vary between these runs'], ...
    ['# but not the mesh nor the model (velocity and density) then we ' ...
    'can also read the mesh and model files'], ...
    ['# from a single run in the beginning and broadcast them to all ' ...
    'the others; for a large number of simultaneous'], ...
    ['# runs for instance when solving inverse problems iteratively ' ...
    'this can DRASTICALLY reduce I/Os to disk in the solver'], ...
    ['# (by a factor equal to NUMBER_OF_SIMULTANEOUS_RUNS), and ' ...
    'reducing I/Os is crucial in the case of huge runs.'], ...
    ['# Thus, always set this option to .true. if the mesh and the ' ...
    'model are the same for all simultaneous runs.'], ...
    ['# In that case there is no need to duplicate the mesh and model ' ...
    'file database (the content of the DATABASES_MPI'], ...
    ['# directories) in each of the run0001, run0002,... directories, ' ...
    'it is sufficient to have one in run0001'], ...
    '# and the code will broadcast it to the others)'};
writecomment(fid, block_comment13);
writebool(fid, 'BROADCAST_SAME_MESH_AND_MODEL', ...
    params.BROADCAST_SAME_MESH_AND_MODEL);
writeblank(fid);

writecomment(fid, ['#-------------------------------------------------' ...
    '----------']);
writeblank(fid);

writecomment(fid, '# set to true to use GPUs');
writebool(fid, 'GPU_MODE', params.GPU_MODE);
writeblank(fid);

writecomment(fid, '# ADIOS Options for I/Os');
writebool(fid, 'ADIOS_ENABLED', params.ADIOS_ENABLED);
writebool(fid, 'ADIOS_FOR_DATABASES', params.ADIOS_FOR_DATABASES);
writebool(fid, 'ADIOS_FOR_MESH', params.ADIOS_FOR_MESH);
writebool(fid, 'ADIOS_FOR_FORWARD_ARRAYS', ...
    params.ADIOS_FOR_FORWARD_ARRAYS);
writebool(fid, 'ADIOS_FOR_KERNELS', params.ADIOS_FOR_KERNELS);
writebool(fid, 'ADIOS_FOR_UNDO_ATTENUATION', ...
    params.ADIOS_FOR_UNDO_ATTENUATION);
writeblank(fid);

writecomment(fid, '# HDF5 Database I/O');
writecomment(fid, ['# (note the flags for HDF5 and ADIOS are mutually ' ...
    'exclusive, only one can be used)']);
writebool(fid, 'HDF5_ENABLED', params.HDF5_ENABLED);
writebool(fid, 'HDF5_FOR_MOVIES', params.HDF5_FOR_MOVIES);
writeint(fid, 'HDF5_IO_NODES', params.HDF5_IO_NODES, ...
    '# HDF5 IO server with number of IO dedicated procs');
writeblank(fid);

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