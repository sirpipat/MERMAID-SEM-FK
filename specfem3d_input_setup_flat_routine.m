function specfem3d_input_setup_flat_routine
% SPECFEM3D_INPUT_SETUP_FLAT_ROUTINE
%
% A script to create input files for several SPECFEM3D simulations. Please
% refer the source code for what it does in detail. I make this for my own
% purposes, and it is subjected to changes.
%
% Last modified by sirawich-at-princeton.edu: 03/25/2025

thisday = datevec(datetime('today'));

testname = sprintf('ORIGIN_TIME_TEST_%04d%02d%02d', thisday(1), ...
    thisday(2), thisday(3));
masterdir = fullfile(getenv('REMOTE3D'), testname);
system(sprintf('mkdir %s', masterdir));

bottom = [4000 6000];
depth = 1500;
freq = 4;
baz = [225 270];
theta = [0 10 40];
fs = 10;
gpu_mode = true;
origin_time = [0 -2 -4 -6 -8];
nsteps = [10000 20000];

t_stf = (-100:0.01:100)';
x_stf_gaussian_upsidedown = -exp(-(t_stf / 0.25).^2);
w = exp(-0.09./(t_stf+0.1).^2) .* (1 + sign(t_stf)) / 2;
x_new = -w .* exp(-t_stf / 2) .* real(sin(2*pi/2*(t_stf).^1.5));

% default, upside-down gaussian, complicated stf
stf = {[], [t_stf, x_stf_gaussian_upsidedown], [t_stf, x_new]};

nn = 0;
for ii = 1:length(bottom)
    for jj = 1:length(baz)
        for kk = 1:length(theta)
            for ll = 1:length(stf)
                for mm = 1:length(origin_time)
                    for pp = 1:length(nsteps)
                        % skip origin_time "greater" than simulation time
                        if origin_time(mm) * (-20) > nsteps(pp) / 100
                            continue
                        end
                        nn = nn + 1;
                        ddir = sprintf(['FK_FLAT_BAZ-%03d_BOTTOM-%04d_' ...
                            'THETA-%02d_STF-%d_ORIGIN_TIME-%02d_' ...
                            'NSTEPS-%05d/'], baz(jj), bottom(ii), ...
                            theta(kk), ll, origin_time(mm), nsteps(pp));
                        ddir = fullfile(masterdir, ddir);
                        fprintf('nn = %3d, %s\n', nn, ddir);
                        specfem3d_input_setup_flat(ddir, bottom(ii), ...
                            depth, freq, baz(jj), theta(kk), fs, ...
                            gpu_mode, stf{ll}, origin_time(mm), ...
                            nsteps(pp));
                    end
                end
            end
        end
    end
end
end