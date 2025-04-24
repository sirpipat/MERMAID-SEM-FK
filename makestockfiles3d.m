function stockdir = makestockfiles3d
% MAKESTOCKFILES3D
%
% Create stock input files for SPECFEM3D runs.
%
% Last modified by sirawich-at-princeton.edu, 04/23/2025

%% create stockfiles folder
stockdir = fullfile(getenv('REMOTE3D'), 'stockfiles');
system(sprintf('mkdir -p %s', stockdir));

fkdir = fullfile(stockdir, 'fkmodelfiles');
system(sprintf('mkdir -p %s', fkdir));

stfdir = fullfile(stockdir, 'stffiles');
system(sprintf('mkdir -p %s', stfdir));

pardir = fullfile(stockdir, 'parfiles');
system(sprintf('mkdir -p %s', pardir));

stadir = fullfile(stockdir, 'stationfiles');
system(sprintf('mkdir -p %s', stadir));

cmtdir = fullfile(stockdir, 'cmtsolutionfiles');
system(sprintf('mkdir -p %s', cmtdir));

meshdir = fullfile(stockdir, 'meshparfiles');
system(sprintf('mkdir -p %s', meshdir));

interfdir = fullfile(stockdir, 'interfacefiles');
system(sprintf('mkdir -p %s', interfdir));

%% Layered model(s)

% Pipatprathanporn & Simons 2024 GJI
layers_2024(1) = struct('rho', 1020, 'vp', 1500, 'vs', 0, 'ztop', 0);
layers_2024(2) = struct('rho', 2500, 'vp', 3400, 'vs', 1963, 'ztop', 4000);

% crustal model drawn from CRUST1.0 entry #40270
layers(1) = struct('rho', 1020, 'vp', 1500, 'vs',    0, 'ztop',      0);
layers(2) = struct('rho', 1820, 'vp', 1750, 'vs',  340, 'ztop',  -4400);
layers(3) = struct('rho', 2550, 'vp', 5000, 'vs', 2700, 'ztop',  -4650);
layers(4) = struct('rho', 2850, 'vp', 6500, 'vs', 3700, 'ztop',  -5130);
layers(5) = struct('rho', 3050, 'vp', 7100, 'vs', 4050, 'ztop',  -6660);
layers(6) = struct('rho', 3330, 'vp', 8080, 'vs', 4490, 'ztop', -11370);

%% source-time function(s)
t_stf = (-100:0.01:100)';
freq = [0.1 0.5 1 2 4 10];

for ii = 1:length(freq)
    % Gaussian
    x_stf = exp(-(t_stf * freq(ii)).^2);
    fname = sprintf('stf_file_Gaussian_%s.txt', float2filestr(freq(ii)));
    writetimeseries(t_stf, x_stf, fullfile(stfdir, fname));
    
    % Ricker
    x_stf = (1 - 2 * freq(ii)^2 * t_stf.^2) .* exp(-(freq(ii) * t_stf).^2);
    fname = sprintf('stf_file_Ricker_%s.txt', float2filestr(freq(ii)));
    writetimeseries(t_stf, x_stf, fullfile(stfdir, fname));
end

%% FKMODEL(s)
theta = [0 10 40];

% Pipatprathanporn & Simons 2024 GJI
for ii = 1:length(freq)
    for jj = 1:length(theta)
        fkmodel.nlayers = 2;
        fkmodel.layers = num2cell(layers_2024);
        fkmodel.wave = 'p';
        fkmodel.baz = 180;
        fkmodel.theta = theta(jj);
        fkmodel.fmax = freq(ii);
        fkmodel.fs = 10;
        fkmodel.twindow = 100;
        fkmodel.amplitude = 1;
        fkmodel.origin_wavefront = [0 -10000 -10000];
        fkmodel.origin_time = 0;
        fkmodel.stf_type = 4;
        fkmodel.stf_file = 'stf_file.txt';
        
        fname = sprintf('FKMODEL_2024_THETA-%2d_FREQ-%s', theta(jj), ...
            float2filestr(freq(ii)));
        writefkmodel(fkmodel, fullfile(fkdir, fname));
    end
end

% Crust
end
