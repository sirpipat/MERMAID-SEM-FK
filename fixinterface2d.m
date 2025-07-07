function fixinterface2d(ddir3dlist)

for ii = 1:length(ddir3dlist)
    % load interface
    [itfs, ~] = loadinterfacefiles3d(fullfile(ddir3dlist{ii}, ...
        'DATA', 'meshfem3D_files', 'interfaces.dat'));
    
    % cross section of bathymetry profile at x=0
    tparams.X = itfs{1}.LAT_MIN + (0:itfs{1}.NETA-1)' * itfs{1}.SPACING_ETA;
    tparams.Z = indeks(itfs{1}.Z, ...
        sprintf('1:%d, ceil(%d/2)', itfs{1}.NETA, itfs{1}.NXI));
    
    % hydrophone depth
    depth = -1500;
    
    % load FKMODEL params
    inject = loadfkmodel(fullfile(ddir3dlist{ii}, 'DATA', 'FKMODEL'));
    inject.freq = inject.fmax;
    time_function_type = 3;
    
    solver.gpu_mode = true;
    
    % convert (X,Z) coordinates to those in SPECFEM2D
    % (0,0) in SPECFEM3D == (10000,9600) in SPECFEM2D 
    tparams.X = tparams.X + 10000;
    tparams.Z = tparams.Z +  9600;
    
    % determine directory for 2D runs
    words = split(ddir3dlist{ii}, filesep);
    ddir2d = fullfile(getenv('REMOTE2D'), sprintf('%s2D', words{end-1}), words{end});
    ddir2d = [ddir2d filesep];
    
    % at some point
    specfem2d_input_setup_response('from3d', ...
                                   'custom', ...
                                   tparams, ...
                                   -depth, ...
                                   'homogeneous', ...
                                   'homogeneous', ...
                                   inject.freq / pi, ... % Gaussian in 2D: exp(-(pi * t * f).^2)
                                   inject.theta, ...
                                   time_function_type, ...
                                   [], ...
                                   ddir2d, ...
                                   false, ...
                                   'devel', ...
                                   solver.gpu_mode);
                               
    % check the result
    plotoceanbottom3d(ddir3dlist{ii});
end
end