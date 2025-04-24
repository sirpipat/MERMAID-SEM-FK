function specfem3d_input_setup_routine(name)
% SPECFEM3D_INPUT_SETUP_ROUTINE(name)
%
% A script to create input files for several SPECFEM3D simulations. Please
% refer the source code for what it does in detail. I make this for my own
% purposes, and it is subjected to changes.
%
% INPUT:
% name          name to give for the input master directory
%               [Default: 'FEATURE_TEST']
%
% Last modified by sirawich-at-princeton.edu: 04/23/2025

defval('name', 'FEATURE_TEST')
thisday = datevec(datetime('today'));

testname = sprintf('%04d%02d%02d_%s', thisday(1), ...
    thisday(2), thisday(3), name);

masterdir3d = fullfile(getenv('REMOTE3D'), testname);
system(sprintf('mkdir %s', masterdir3d));
masterdir2d = fullfile(getenv('REMOTE2D'), testname);
system(sprintf('mkdir %s', masterdir2d));

% list of parameters
depth = -1500;
bottom = -5000;
theta = 0;
freq = 4;
fs = 10;
solver.gpu_mode = true;
solver.nsteps = 14000;

% FLAT
for ii = 1:length(theta)
    for jj = 1:length(freq)
        t_stf = (-20:0.01:100)';
        x_stf = exp(-(t_stf * freq(jj)).^2);
        inject.stf = [t_stf, x_stf];
        for kk = 1:length(fs)
            for ll = 1:length(bottom)
                ddir = sprintf('FLAT_BOTTOM-%04d_THETA-%02d_FREQ-%04.1f_FS-%02d', ...
                    -bottom(ll), theta(ii), freq(jj), fs(kk));
                ddir = fullfile(masterdir3d, ddir);

                topo.type = 'flat';
                topo.depth = bottom(ll);

                inject.type = 'stf';
                inject.theta = theta(ii);
                inject.freq = freq(jj);
                inject.fs = fs(kk);
                specfem3d_input_setup(ddir, depth, topo, inject, solver, true);
            end
        end
    end
end

% RANDOM FIELD
s2 = [100 10000 250000];
nu = [0.5 1 1.5];
rho = [500 1000 2000];
for mm = 1:length(freq)
    t_stf = (-20:0.01:100)';
    x_stf = exp(-(t_stf * freq(mm)).^2);
    inject.stf = [t_stf, x_stf];
    for ii = 1:length(s2)
        for jj = 1:length(nu)
            for kk = 1:length(rho)
                for ll = 1:length(theta)
                    for it = 1:100 % repeating iterations
                        % skip "problematic" random field param combination
                        if nu(jj) > 1 && rho(kk) > 1500
                            continue
                        end
                        if nu(jj) < Inf
                            ddir = sprintf('RANDOM_S2-%05d_NU-%3.1f_RHO-%04d_THETA-%02d_FREQ-%04.1f_IT-%03d', ...
                                s2(ii), nu(jj), rho(kk), theta(ll), freq(mm), it);
                        else
                            ddir = sprintf('RANDOM_S2-%05d_NU-inf_RHO-%04d_THETA-%02d_FREQ-%04.1f_IT-%03d', ...
                                s2(ii), rho(kk), theta(ll), freq(mm), it);
                        end

                        ddir = fullfile(masterdir3d, ddir);

                        topo.type = 'randomfield';
                        topo.depth = bottom(1);
                        topo.s2 = s2(ii);
                        topo.nu = nu(jj);
                        topo.rho = rho(kk);
                        topo.dydx = [500 500];
                        topo.NyNx = [57 57];
                        topo.blurs = Inf;
                        topo.taper = 1;

                        inject.type = 'stf';
                        inject.theta = theta(ll);
                        inject.freq = freq(mm);
                        inject.fs = 10;
                        try
                            specfem3d_input_setup(ddir, depth, topo, inject, solver, true);
                        catch
                            continue
                        end
                    end
                end
            end
        end
    end
end